/**
ecf.cpp

This target performs the modal analysis of an Exposed Core Fiber (ECF).

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
#include <pzintel.h>
#include <pzmultiphysicselement.h>
#include <Electromagnetics/TPZWaveguideModalAnalysis.h>
#include <pzbuildmultiphysicsmesh.h>
#include <TPZSimpleTimer.h>              //for TPZSimpleTimer
#include <pzshapecube.h>

TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const int neigenpairs, const CSTATE target,
            TPZEigenSort sorting, bool usingSLEPC);

//Sets geometric info regarding all the circles in the mesh
TPZVec<wgma::gmeshtools::ArcData> SetUpArcData(const REAL scale);

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
  matmap["core"] = std::make_pair<CSTATE,CSTATE>(1.45*1.45,1.);
  std::map<std::string,wgma::bc::type> bcmap;
  bcmap["bound"] = wgma::bc::type::PEC;
  // operational wavelength
  constexpr STATE lambda{50e-6};
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
  constexpr int nEigenpairs{10};
  //whether to compute eigenvectors (instead of just eigenvalues)
  constexpr bool computeVectors{true};
  //how to sort the computed eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::RealAscending};
  /*
   The simulation uses a Krylov-based Arnoldi solver for solving the
   generalised EVP. A shift-and-inverse spectral transform is applied in 
   the system for better convergence of the eigenvalues.
   The target variable should be close to the desired eigenvalue (in this case,
   the effective index neff).*/
  constexpr CSTATE target = -2.12;

  //PML attenuation constant
  constexpr STATE alphaPML{0.02};

  /*********************
   * exporting options *
   *********************/

  //whether to print the geometric mesh in .txt and .vtk formats
  constexpr bool printGMesh{true};
  //whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  //whether to export the eigenvalues in .csv format
  constexpr bool exportCsv{true};
  //prefix for exported files
  const std::string prefix{"ecf"};
  //resolution of the .vtk file in which the solution will be exported
  constexpr int vtkRes{1};
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
  const std::string filename{"meshes/ecf.msh"};
  TPZVec<std::map<std::string,int>> gmshmats;
  auto gmesh = wgma::gmeshtools::ReadGmshMesh(filename, scale, gmshmats);

  /*
    In order to represent curved geometries correctly, we need to use
    TPZArc3D, which uses an exact mapping to represent circumference arcs.
    The neighbouring elements are also modified to take this deformation into account.
   */
  auto arcdata = SetUpArcData(scale);

  wgma::gmeshtools::SetExactArcRepresentation(gmesh, arcdata);
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

  if(exportCsv){
    const std::string csvfile = prefix+"_eigenvalues.csv";
    analysis.WriteToCsv(csvfile, lambda);
  }
  
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

  constexpr int krylovDim{50};
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

TPZVec<wgma::gmeshtools::ArcData> SetUpArcData(const REAL scale)
{
  TPZVec<wgma::gmeshtools::ArcData> arcdata(7);

  constexpr REAL um = 1e-6;
  const REAL r_core{14*um / scale};
  const REAL r_clad{(12.5*um) / scale}, r_clad_out{(100*um) / scale};

  //rotates a pt around the origin in the xy-plane
  auto RotatePt = [] (const REAL theta, const REAL x, const REAL y){
    const REAL xrot = cos(theta) * x - sin(theta) * y;
    const REAL yrot = sin(theta) * x + cos(theta) * y;
    return std::pair(xrot, yrot);
  };

  const REAL xini{1.05*r_clad+r_core}, yini{0};
  constexpr int ncircs{6};
  for(int ic = 0; ic < ncircs; ic++){
    const REAL theta = ic * M_PI/3;
    REAL xc, yc;
    std::tie(xc, yc) = RotatePt(theta, xini,yini);
    arcdata[ic].m_matid = 20 + ic+1;
    arcdata[ic].m_radius = r_clad;
    arcdata[ic].m_xc = xc;
    arcdata[ic].m_yc = yc;
    arcdata[ic].m_zc = 0;
  }

  arcdata[6].m_matid = 27;
  arcdata[6].m_radius = r_clad_out;
  arcdata[6].m_xc = 0;
  arcdata[6].m_yc = 0;
  arcdata[6].m_zc = 0;
  
  return arcdata;
}