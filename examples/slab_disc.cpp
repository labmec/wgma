/**
slab_disc.cpp

This target performs the modal analysis of a dielectric slab
and then the subsequent scattering analysis at a waveguide discontinuity.
***/

// wgma includes
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include <wganalysis.hpp>
#include <scattering.hpp>
#include <util.hpp>
#include <slepcepshandler.hpp>
#include <post/orthosol.hpp>
// pz includes
#include <MMeshType.h>      //for MMeshType
#include <TPZSimpleTimer.h> //for TPZSimpleTimer
#include <pzlog.h>          //for TPZLogger
#include <TPZVTKGenerator.h>

#include <regex>//for string search

TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const CSTATE target, const int neigenpairs,
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
  // operational wavelength

  // the meshes were designed in micrometers, so lambda has to follow
  constexpr STATE lambda{1.5};

  
  constexpr STATE ncore{1.5};
  constexpr STATE nclad{1.000};
  /*
    Given the small dimensions of the domain, scaling it can help in
    achieving good precision. Using 1./k0 as a scale factor results in
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  constexpr REAL scale{lambda / (2 * M_PI)};

  constexpr wgma::planarwg::mode mode{wgma::planarwg::mode::TE};

  /******************
   *  fem options   *
   ******************/
  // polynomial order to be used in the approximation
  constexpr int pOrder{3};

  /******************
   * solver options *
   ******************/

  // number of threads to use
  constexpr int nThreads{8};
  // how to sort eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  constexpr bool usingSLEPC {false};

  /*********************
   * exporting options *
   *********************/

  // whether to print the geometric mesh in .txt and .vtk formats
  constexpr bool printGMesh{true};
  // whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  // path for output files
  const std::string path {"res_slab_disc/"};
  // common prefix for both meshes and output files
  const std::string basisName{"slab_disc"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(prefix));


  constexpr int nEigenpairs{5};

  constexpr CSTATE target{ncore*ncore};
  
  // resolution of the .vtk file in which the solution will be exported
  constexpr int vtkRes{0};

  /********************
   * advanced options *
   ********************/

  // reorder the equations in order to optimize bandwidth
  constexpr bool optimizeBandwidth{false};
  /*
    The equations corresponding to homogeneous dirichlet boundary condition(PEC)
    can be filtered out of the global system in order to achieve better
    conditioning.
   */
  constexpr bool filterBoundaryEqs{true};

  /*********
   * begin *
   *********/




  /*************
   * geometry  *
   ************/
  // scoped-timer
  TPZSimpleTimer total("Total");

  // creates gmesh
  // file containing the .msh mesh
  const std::string meshfile{"meshes/"+basisName+".msh"};

  /**
     in order to exactly represent all the circles in the mesh,
     the python script that generated the mesh also generated this .csv file
  **/
  const std::string arcfile{"meshes/"+basisName+"_circdata.csv"};
  TPZVec<std::map<std::string, int>> gmshmats;
  constexpr bool verbosity_lvl{false};
  auto gmesh = wgma::gmeshtools::ReadGmshMesh(meshfile, scale,
                                              gmshmats,verbosity_lvl);

  // print wgma_gmesh to .txt and .vtk format
  if (printGMesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }


  /**************************
   * cmesh(modal analysis)  *
   **************************/

  auto modal_cmesh = [gmesh,&gmshmats,scale](){
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData modal_data;

    std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
    modal_mats["source_clad_left"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad, 1.);
    modal_mats["source_core_left"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    std::map<std::string, wgma::bc::type> modal_bcs;
    modal_bcs["source_left_bound"] = wgma::bc::type::PEC;
    //dimension of the modal analysis 
    constexpr int modal_dim{1};
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs, {},
                                            {}, modal_data, modal_dim);

    return wgma::wganalysis::CMeshWgma1D(gmesh,mode,pOrder+1,modal_data,
                                         lambda, scale);
  }();

  /*************************
   * solve(modal analysis) *
   *************************/

  std::streamsize ss = std::cout.precision();
  
  std::cout.precision(std::numeric_limits<STATE>::max_digits10);
  STATE rel_error{0};
  constexpr STATE tol = std::numeric_limits<STATE>::epsilon()*1000;

  constexpr bool computeVectors{true};
  auto solver = SetupSolver(target, nEigenpairs, sortingRule, usingSLEPC);

  wgma::wganalysis::WgmaPlanar
    modal_an(modal_cmesh, nThreads,
             optimizeBandwidth, filterBoundaryEqs);
  

  {
    TPZSimpleTimer modal_analysis("Modal analysis");

    solver->SetTarget(target);
    modal_an.SetSolver(*solver);
    modal_an.Assemble(TPZEigenAnalysis::Mat::A);
    modal_an.Assemble(TPZEigenAnalysis::Mat::B);
    modal_an.Solve(computeVectors);
    //load all obtained modes into the mesh
    modal_an.LoadAllSolutions();
    //leave empty for all valid matids
    std::set<int> matids {};
    auto ortho = wgma::post::OrthoSol(modal_cmesh, matids, nThreads);
    //orthogonalise the modes
    auto normsol = ortho.Orthogonalise();
    //let us set the orthogonalised modes
    modal_an.SetEigenvectors(normsol);
  }


  {
    TPZSimpleTimer postProc("Post processing");
    
    if(exportVtk){
      const std::string modal_file = prefix+"_modal";
      TPZVec<std::string> fvars = {
        "Field_real",
        "Field_imag",
        "Field_abs",
        // "Field_phase",
        "Deriv_real",
        "Deriv_imag",
        "Deriv_abs",
        // "Deriv_phase"
      };
      auto vtk = TPZVTKGenerator(modal_cmesh, fvars, modal_file, vtkRes);
      const auto nsol = modal_an.GetEigenvalues().size();
      for(auto isol = 0; isol < nsol; isol++){
        modal_an.LoadSolution(isol);
        vtk.Do();
      }
    }
  }  
  /*********************
   * cmesh(scattering) *
   *********************/

  //materials that will represent our source
  const std::string srcMat[] = {"source_core_left", "source_clad_left"};
  std::set<int> src_ids;
  for(const auto &mat : srcMat){
    src_ids.insert(gmshmats[1].at(mat));
  }
  

  auto scatt_cmesh = [gmesh,&gmshmats, src_ids](){
    
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData scatt_data;
    std::map<std::string, std::pair<CSTATE, CSTATE>> scatt_mats;
    scatt_mats["core_left"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    scatt_mats["core_right"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    scatt_mats["cladding_left"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad,1.);
    scatt_mats["cladding_right"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad,1.);
    std::map<std::string, wgma::bc::type> scatt_bcs;
    scatt_bcs["scatt_bound"] = wgma::bc::type::PEC;

    constexpr STATE alphaPML {1.5};
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, scatt_mats, scatt_bcs, alphaPML,
                                            alphaPML, scatt_data);
    

    return wgma::scattering::CMeshScattering2D(gmesh, mode, pOrder, scatt_data,src_ids,
                                               lambda,scale);
  }();
  /*********************
   * solve(scattering) *  
   *********************/  


  
  TPZSimpleTimer tscatt("Scattering");
  auto scatt_an = wgma::scattering::Analysis(scatt_cmesh, nThreads,
                                             optimizeBandwidth, filterBoundaryEqs);

  //in the first run we assemble the whole algebraic system. afterwards we can
  //just compute the rhs
  bool firstrun=true;

  //set up post processing
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_imag",
    "Field_abs"};
  const std::string scatt_file = prefix+"_scatt";
  auto vtk = TPZVTKGenerator(scatt_cmesh, fvars, scatt_file, vtkRes);

  //get id of source materials
  wgma::scattering::SourceWgma src;
  src.id = src_ids;
  src.modal_cmesh = modal_cmesh;
  //indexes of modal solutions to be used as sources
  TPZVec<int> src_index = {0,1};
  const int nsols = src_index.size();
    
  for(int isol = 0; isol < nsols; isol++){
    std::cout<<"running source "<<isol+1<<" out of "<<nsols<<std::endl;
    modal_an.LoadSolution(src_index[isol]);
    auto beta = std::sqrt(modal_an.GetEigenvalues()[isol]);
    wgma::scattering::LoadSource(scatt_cmesh, src);
    wgma::scattering::SetPropagationConstant(scatt_cmesh, beta);
      
    if(firstrun){
      firstrun=false;
      scatt_an.Assemble();
    }else{
      scatt_an.AssembleRhs(src_ids);
    }
    scatt_an.Solve();
    vtk.Do();
  }
}

#include <slepcepshandler.hpp>
#include <TPZKrylovEigenSolver.h>
//utility functions
TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const CSTATE target,const int neigenpairs,
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
  const int krylovDim = 2*neigenpairs;
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
    constexpr bool eps_verbosity = false;
    
    
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
    krylov_solver->SetKrylovDim(krylovDim);
    krylov_solver->SetNEigenpairs(neigenpairs);
    krylov_solver->SetAsGeneralised(true);
    
    solver = krylov_solver;
  }
  
  solver->SetEigenSorting(sorting);
  return solver;
}