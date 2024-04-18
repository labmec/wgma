/**
stepfiber.cpp

This target performs the modal analysis of a step-index optical fiber.

It also illustrates the usage of the SLEPc solver and how to create
complex curved structured meshes in NeoPZ.
***/

//wgma includes
#include "wganalysis.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include "slepcepshandler.hpp"
#include "util.hpp"
//pz includes
#include <pzcmesh.h>                     //for TPZCompMesh
#include <pzgmesh.h>                     //for TPZGeoMesh
#include <pzlog.h>                       //for TPZLogger
#include <TPZElectromagneticConstants.h> //for pzelectromag::cZero
#include <TPZKrylovEigenSolver.h>        //for TPZKrylovEigenSolver
#include <TPZVTKGenerator.h>
#include <TPZPardisoSolver.h>
#include <post/solutionnorm.hpp>
#include <TPZSimpleTimer.h>              //for TPZSimpleTimer

// Sets geometric info regarding all the circles in the mesh
TPZVec<wgma::gmeshtools::ArcData> SetUpArcData(std::string_view filename,
                                               const REAL scale);

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
  
  //refractive index of the fibers core
  constexpr STATE n_core{1.4457};
  //refractive index of the fibers cladding
  constexpr STATE n_clad{1.4378};
  //magnetic permeability of the core
  constexpr STATE ur_core{1};
  //magnetic permeability of the cladding
  constexpr STATE ur_clad{1};
  // operational wavelength
  constexpr STATE lambda{1.55};
  /*
    Given the small dimensions of the domain, scaling it can help in
    achieving good precision. Using 1./k0 as a scale factor results in
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  constexpr REAL scale{lambda / (2 * M_PI)};

  /******************
   *  fem options   *
   ******************/
  // polynomial order to be used in the modal analysis
  constexpr int pOrder2D{2};
  
  constexpr STATE modal_alphaPMLx{0.04};
  constexpr STATE modal_alphaPMLy{0.04};
  /******************
   * solver options *
   ******************/

  // number of threads to use
  const int nThreads = std::thread::hardware_concurrency();
  const int nThreadsDebug = std::thread::hardware_concurrency();
  // how to sort eigenvaluesn
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  bool usingSLEPC {false};
  constexpr int nEigenvalues{10};
  constexpr CSTATE target = -2.09;//n_core*n_core;
  constexpr bool symMat{false};

  constexpr bool computeVectors{true};
  /*********************
   * exporting options *
   *********************/

  // whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  // path for output files
  const std::string path {"res_stepfiber/"};
  // common prefix for both meshes and output files
  const std::string basisName{"stepfiber"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(prefix));

  
  // resolution of the .vtk file in which the solution will be exported
  constexpr int vtkRes{0};

  /********************
   * advanced options *
   ********************/

  // reorder the equations in order to optimize bandwidth
  constexpr bool optimizeBandwidth{true};
  /*
    The equations corresponding to homogeneous dirichlet boundary condition(PEC)
    can be filtered out of the global system in order to achieve better
    conditioning.
   */
  constexpr bool filterBoundaryEqs{true};
  /*
    Whether to use a non-linear representation for cylinders
   */
  constexpr bool cyl3D{true};

  //! Whether to use a direct solver
  constexpr bool direct{false};

  constexpr bool printGMesh{false};

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

  
  TPZVec<std::map<std::string, int>> gmshmats;
  constexpr bool verbosity_lvl{false};
  auto gmesh = wgma::gmeshtools::ReadGmshMesh(meshfile, scale,
                                              gmshmats,verbosity_lvl);

  if(cyl3D){
    /**
     in order to exactly represent all the circles in the mesh,
     the python script that generated the mesh also generated this .csv file
    **/
    const std::string arcfile{"meshes/"+basisName+"_circdata.csv"};
    auto arcdata = SetUpArcData(arcfile, scale);
    wgma::gmeshtools::SetExactArcRepresentation(gmesh, arcdata);
  }
  
  //print gmesh to .txt and .vtk format
  if(printGMesh)
  {
    //prefix for the gmesh files
    const auto filename = prefix+"_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh,filename);
  }

  //setting up cmesh data
  wgma::cmeshtools::PhysicalData modal_data;
  auto modal_cmesh = [gmesh,&gmshmats, &modal_data, scale](){
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData modal_data;

    std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
    modal_mats["core"] = std::make_pair<CSTATE, CSTATE>(n_core*n_core, 1.);
    modal_mats["cladding"] = std::make_pair<CSTATE, CSTATE>(n_clad*n_clad, 1.);
    std::map<std::string, wgma::bc::type> modal_bcs;
    modal_bcs["modal_bnd"] = wgma::bc::type::PEC;

    //dimension of the modal analysis 
    constexpr int modal_dim{2};
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                            {modal_alphaPMLx,modal_alphaPMLy},
                                            modal_data, modal_dim);
    return wgma::wganalysis::CMeshWgma2D(gmesh,pOrder2D, modal_data,lambda, scale);
  }();

  
  //WGAnalysis class is responsible for managing the modal analysis
  wgma::wganalysis::Wgma2D analysis(modal_cmesh,nThreads,optimizeBandwidth,filterBoundaryEqs);
  
  auto solver = SetupSolver(nEigenvalues,target, sortingRule, usingSLEPC);
  
  analysis.SetSolver(solver);
  analysis.Assemble();

  analysis.Solve(computeVectors);
  analysis.LoadAllSolutions();

  if (!computeVectors && !exportVtk) return 0;

  const std::string suffix = usingSLEPC ? "_slepc" : "_krylov";
  const std::string plotfile = path+"stepfiber_field" + suffix;

  
  {
    
    TPZSimpleTimer tpostprocess("Post processing");
    TPZVec<std::string> fvars = {
      "Ez_real",
      "Ez_abs",
      "Et_real",
      "Et_abs"};
    auto vtk = TPZVTKGenerator(modal_cmesh[0], fvars, plotfile, vtkRes);
    auto ev = analysis.GetEigenvalues();
    for (int isol = 0; isol < ev.size(); isol++) {
      auto currentKz = std::sqrt(-1.0*ev[isol]);
      std::cout<<"\rPost processing step "<<isol+1<<" out of "<<ev.size()
               <<"(kz = "<<currentKz<<")"<<std::endl;
      analysis.LoadSolution(isol);
      vtk.Do();
    }
  }
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
    
    constexpr PC pc = PC::LU;
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

// Sets geometric info regarding all the circles in the mesh
TPZVec<wgma::gmeshtools::ArcData> SetUpArcData(std::string_view filename,
                                               const REAL scale) {
  std::ifstream read(filename.data());
  if (!read) {
    std::cout << "Couldn't open the file " << filename << std::endl;
    DebugStop();
  }

  auto getNextLineAndSplitIntoTokens =
      [](std::istream &str) -> std::vector<std::string> {
    std::vector<std::string> result;
    std::string line;
    std::getline(str, line);

    std::stringstream lineStream(line);
    std::string cell;

    while (std::getline(lineStream, cell, ',')) {
      result.push_back(cell);
    }
    // This checks for a trailing comma with no data after it.
    if (!lineStream && cell.empty()) {
      // If there was a trailing comma then add an empty element.
      result.push_back("");
    }
    return result;
  };

  auto line = getNextLineAndSplitIntoTokens(read); // header
  line = getNextLineAndSplitIntoTokens(read);
  // we expect xc, yc, zc, r (in um), and matid
  TPZVec<wgma::gmeshtools::ArcData> arcs;
  const auto factor = 1./scale;
  while (line.size() == 5) {
    wgma::gmeshtools::ArcData arc;

    arc.m_xc = std::stod(line[0]) * factor;
    arc.m_yc = std::stod(line[1]) * factor;
    arc.m_zc = std::stod(line[2]) * factor;
    arc.m_radius = std::stod(line[3]) * factor;
    arc.m_matid = std::stoi(line[4]);
    const int narcs = arcs.size();
    arcs.Resize(narcs + 1);
    arcs[narcs] = arc;

    line = getNextLineAndSplitIntoTokens(read);
  }
  return arcs;
}