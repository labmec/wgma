/**
stepfiber.cpp

This target performs the modal analysis of a double-ridged waveguide.

It also illustrates how to import a .gmsh mesh and how to perform
directional mesh refinement.
***/

//wgma includes
#include "wganalysis.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include "slepcepshandler.hpp"
//pz includes
#include <MMeshType.h>                   //for MMeshType
#include <pzcmesh.h>                     //for TPZCompMesh
#include <pzgmesh.h>                     //for TPZGeoMesh
#include <pzlog.h>                       //for TPZLogger
#include <TPZElectromagneticConstants.h> //for pzelectromag::cZero
#include <TPZKrylovEigenSolver.h>        //for TPZKrylovEigenSolver
#include <TPZRefPatternTools.h>          //for TPZRefPatternTools
#include "TPZRefPatternDataBase.h"       //for TPZRefPatternDataBase



#include <TPZSimpleTimer.h>              //for TPZSimpleTimer


TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const int neigenpairs, const CSTATE target,
            TPZEigenSort sorting, bool usingSLEPC);

int main(int argc, char *argv[]) {
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif

  /***********************
   * setting the problem *
   ***********************/

  std::map<std::string,std::pair<CSTATE,CSTATE>> matmap;
  matmap["air"] = std::make_pair<CSTATE,CSTATE>(1.,1.);
  matmap["core"] = std::make_pair<CSTATE,CSTATE>(3.44*3.44,1.);
  matmap["substrate"] = std::make_pair<CSTATE,CSTATE>(3.40*3.40,1.);
  std::map<std::string,wgma::bc::type> bcmap;
  bcmap["bound"] = wgma::bc::type::PEC;
  // operational wavelength
  constexpr STATE lambda{1.15e-6};
  /*
    Given the small dimensions of the domain, scaling it can help in 
    achieving good precision. Using 1./k0 as a scale factor results in 
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  constexpr REAL scale{lambda/(2*M_PI)};

  /******************
   *  fem options   *
   ******************/
  // polynomial order to be used in the approximation
  constexpr int pOrder{2};

  /******************
   * solver options *
   ******************/
  
  //number of threads to use
  constexpr int nThreads{8};
  //number of genvalues to be computed
  constexpr int nEigenpairs{4};
  //whether to compute eigenvectors (instead of just eigenvalues)
  constexpr bool computeVectors{true};
  //how to sort the computed eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  /*
   The simulation uses a Krylov-based Arnoldi solver for solving the
   generalised EVP. A shift-and-inverse spectral transform is applied in 
   the system for better convergence of the eigenvalues.
   The target variable should be close to the desired eigenvalue (in this case,
   the effective index neff).*/
  constexpr CSTATE target = -12;

  //PML attenuation constant
  constexpr STATE alphaPML{0.04};

  /*********************
   * exporting options *
   *********************/

  //whether to print the geometric mesh in .txt and .vtk formats
  constexpr bool printGMesh{true};
  //whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  //prefix for exported files
  const std::string prefix{"ribwg"};
  //resolution of the .vtk file in which the solution will be exported
  constexpr int vtkRes{0};
  //if true, the real part of the electric fields is exported. otherwise, the magnitude
  constexpr bool printRealPart{true};

  /********************
   * advanced options *
   ********************/
  
  //reorder the equations in order to optimize bandwidth
  constexpr bool optimizeBandwidth{true};
  /*
    The equations corresponding to homogeneous dirichlet boundary condition(PEC)
    can be filtered out of the global system in order to achieve better conditioning.
   */
  constexpr bool filterBoundaryEqs{true};

  /*********
   * begin *
   *********/

  
  //scoped-timer 
  TPZSimpleTimer total("Total");
  

  bool usingSLEPC{true};
  //creates gmesh
  //filename
  const std::string filename{"meshes/ribWG.msh"};
  TPZVec<std::map<std::string,int>> gmshmats;
  auto gmesh = wgma::gmeshtools::ReadGmshMesh(filename, scale, gmshmats);

  /*
    Since we know that there is a singularity at the metallic corner,
    we will directionally refine all the elements with a node in this corner.
   */
  {
    TPZSimpleTimer refinement("Refining mesh");
    /*
      We initialise now the database of refinement patterns of 1 and 2 dimensions
      (lines, triangles and quadrilaterals)
     */
    gRefDBase.InitializeRefPatterns(1);
    gRefDBase.InitializeRefPatterns(2);
    //identifier of the corner where the singularity resides
    const auto cornermatid = gmshmats[0].at("corner");
    
    constexpr int nrefsteps{5};

    std::cout << "Refining..."<<std::endl;
    
    for(int iref = 0; iref < nrefsteps; iref++){
      std::set<int> matids = {cornermatid};
      const int nels = gmesh->NElements();
      for(int el = 0; el < nels; el++){
        auto *gel = gmesh->Element(el);
        if(gel && gel->NSubElements() == 0){
          TPZRefPatternTools::RefineDirectional(gel, matids);
        }
      }
    }
    std::cout << "\nFinished refining!"<<std::endl;
  }
  //print gmesh to .txt and .vtk format
  if(printGMesh)
  {
    //prefix for the gmesh files
    const std::string filename = prefix +"_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh,filename);
  }

  
  //setting up cmesh data
  TPZVec<int> volMatIdVec;
  TPZVec<CSTATE> urVec;
  TPZVec<CSTATE> erVec;
  TPZVec<wgma::pml::data> pmlDataVec;
  TPZVec<wgma::bc::data> bcDataVec;

  
  wgma::cmeshtools::SetupGmshMaterialData(gmshmats, matmap, bcmap,
                                          alphaPML, volMatIdVec,
                                          erVec, urVec,
                                          pmlDataVec, bcDataVec);
  /*
   The problem uses an H1 approximation space for the longitudinal component 
   and a HCurl approximation space for the transversal one. Therefore, three
  // computational meshes are generated. One for each space and a multiphysics mesh*/
  auto meshVec = wgma::cmeshtools::CreateCMesh(gmesh,pOrder,volMatIdVec,
                                               urVec, erVec, pmlDataVec,
                                               bcDataVec, lambda,scale);

  
  //WGAnalysis class is responsible for managing the modal analysis
  wgma::WGAnalysis analysis(meshVec,nThreads,optimizeBandwidth,filterBoundaryEqs);
  
  auto solver = SetupSolver(nEigenpairs, target, sortingRule, usingSLEPC);

  analysis.SetSolver(solver);
  analysis.Run(computeVectors);
  
  if (!computeVectors && !exportVtk) return 0;

  const std::string plotfile = prefix+"_field_";

  TPZSimpleTimer postProc("Post processing");
  
  analysis.PostProcess(plotfile, vtkRes, printRealPart);
  return 0;
}



TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const int neigenpairs, const CSTATE target,
            TPZEigenSort sorting, bool usingSLEPC)
{

#ifndef WGMA_USING_SLEPC
  if(usingSLEPC){
    std::cout<<"wgma was not configured with slepc. defaulting to: "
             <<"TPZKrylovSolver"<<std::endl;
    usingSLEPC = false;
  }
#endif

  TPZAutoPointer<TPZEigenSolver<CSTATE>> solver{nullptr};

  constexpr int krylovDim{-1};
  if (usingSLEPC){
    using namespace wgma::slepc;
    /*
      The following are suggested SLEPc settings.
      NOTE: -1 stands for PETSC_DECIDE
    */
    
    constexpr STATE eps_tol = -1;//PETSC_DECIDE
    constexpr int eps_max_its = -1;//PETSC_DECIDE
    constexpr EPSConv eps_conv_test = EPSConv::EPS_CONV_REL;
    constexpr EPSWhich eps_which = EPSWhich::EPS_TARGET_REAL;
    
    constexpr Precond pc = Precond::LU;
    constexpr KSPSolver linsolver = KSPSolver::PREONLY;
    constexpr STATE ksp_rtol = -1;//PETSC_DECIDE
    constexpr STATE ksp_atol = -1;//PETSC_DECIDE
    constexpr STATE ksp_dtol = -1;//PETSC_DECIDE
    constexpr STATE ksp_max_its = -1;//PETSC_DECIDE
    constexpr bool eps_true_residual = false;
    constexpr EPSProblemType eps_prob_type = EPSProblemType::EPS_GNHEP;//do NOT change
    constexpr EPSType eps_solver_type = EPSType::KRYLOVSCHUR;
    constexpr bool eps_krylov_locking = true;
    constexpr STATE eps_krylov_restart = 0.7;
    constexpr STATE eps_mpd = -1;//PETSC_DECIDE
    constexpr bool eps_verbosity = true;
    
    
    auto eps_solver = new EPSHandler<CSTATE>;
    eps_solver->SetType(eps_solver_type);
    eps_solver->SetProblemType(eps_prob_type);
    eps_solver->SetEPSDimensions(neigenpairs, krylovDim, eps_mpd);
    eps_solver->SetWhichEigenpairs(eps_which);
    eps_solver->SetTarget(target);
    eps_solver->SetTolerances(eps_tol,eps_max_its);
    eps_solver->SetConvergenceTest(eps_conv_test);
    eps_solver->SetKrylovOptions(eps_krylov_locking,eps_krylov_restart);
    eps_solver->SetVerbose(eps_verbosity);
    eps_solver->SetTrueResidual(eps_true_residual);
    
    eps_solver->SetLinearSolver(linsolver);
    eps_solver->SetLinearSolverTol(ksp_rtol,ksp_atol,ksp_dtol,ksp_max_its);
    eps_solver->SetPrecond(pc, 1e-16);

    solver = eps_solver;
  }else{
    auto krylov_solver = new TPZKrylovEigenSolver<CSTATE>;
    TPZSTShiftAndInvert<CSTATE> st;
    krylov_solver->SetSpectralTransform(st);
    krylov_solver->SetTarget(target);
    krylov_solver->SetNEigenpairs(neigenpairs);
    krylov_solver->SetAsGeneralised(true);
    krylov_solver->SetKrylovDim(krylovDim);
    
    solver = krylov_solver;
  }
  
  solver->SetEigenSorting(sorting);
  return solver;
}

