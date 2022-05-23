/**
pcwg.cpp

This target performs the modal analysis of a photonic crystal fiber
and then the subsequent scattering analysis.

Data taken from:
Yasuhide Tsuji and Masanori Koshiba,
"Finite Element Method Using Port Truncation by
 Perfectly Matched Layer Boundary Conditions for
 Optical Waveguide Discontinuity Problems,"
J. Lightwave Technol. 20, 463- (2002)

***/

// wgma includes
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include <wganalysis.hpp>
#include <scattering.hpp>
#include <util.hpp>
// pz includes
#include <MMeshType.h>      //for MMeshType
#include <TPZSimpleTimer.h> //for TPZSimpleTimer
#include <pzlog.h>          //for TPZLogger
// std includes
#include <regex>




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
  // operational wavelength

  constexpr STATE lambda{1.55e-6};
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
  constexpr int pOrder{1};

  /******************
   * solver options *
   ******************/

  // number of threads to use
  constexpr int nThreads{8};
  // how to sort eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  constexpr int nEigenpairs = 1;
  constexpr bool usingSLEPC {true};

  /*********************
   * exporting options *
   *********************/

  // whether to print the geometric mesh in .txt and .vtk formats
  constexpr bool printGMesh{true};
  // whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  // path for output files
  const std::string path {"res_planar_wg_periodic/"};
  // common prefix for both meshes and output files
  const std::string basisName{"periodic_test"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(prefix));


  //attenuation constant for the PML
  constexpr STATE alphaxPML {10.00};
  constexpr STATE alphayPML = alphaxPML/2;
  
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

  
  TPZVec<std::map<std::string, int>> gmshmats;
  std::map<int64_t,int64_t> periodic_els;
  constexpr bool verbosity_lvl{false};
  auto gmesh = wgma::gmeshtools::ReadPeriodicGmshMesh(meshfile, scale,
                                                      gmshmats,periodic_els,
                                                      verbosity_lvl);


  // print wgma_gmesh to .txt and .vtk format
  if (printGMesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }


  /**************************
   * cmesh(modal analysis)  *
   **************************/

  auto modal_cmesh = [&gmesh,&gmshmats,&periodic_els,alphaxPML, alphayPML](){
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData modal_data;

    std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
    modal_mats["subst_1"] = std::make_pair<CSTATE, CSTATE>(1., 1.);
    modal_mats["strip_1"] = std::make_pair<CSTATE, CSTATE>(1.55 * 1.55,1.);
    std::map<std::string, wgma::bc::type> modal_bcs;
    modal_bcs["gamma_1"] = wgma::bc::type::PERIODIC;
    modal_bcs["gamma_2"] = wgma::bc::type::PERIODIC;
    modal_bcs["modal_bound"] = wgma::bc::type::PEC;

    constexpr STATE alphaxPML {10.00};
    constexpr STATE alphayPML = alphaxPML/2;
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs, alphaxPML,
                                            alphayPML,modal_data);

    /*
     let us filter out the undesired PMLs*/
    
    std::vector<wgma::pml::data> modal_pmls;
    constexpr auto pattern{"modal"};
    auto rx = std::regex{ pattern,std::regex_constants::icase };
    for(auto & pml : modal_data.pmlvec){
      if(std::regex_search(pml.name, rx)){
        modal_pmls.push_back(pml);
        std::cout<<"inserting pml "<<pml.name<<" in modal analysis"<<std::endl;
      }
    }

    modal_data.pmlvec = modal_pmls;
    return wgma::wganalysis::CMeshWgmaPeriodic(gmesh,mode,pOrder,modal_data,
                                               periodic_els,lambda,scale);
  }();

  /*************************
   * solve(modal analysis) *
   *************************/
  wgma::wganalysis::WgmaPeriodic2D
    modal_an(modal_cmesh, nThreads,
             optimizeBandwidth, filterBoundaryEqs);

  //this is a non-linear problem on beta, we start with 0 as initial value
  CSTATE beta = 1.55;
  auto solver = SetupSolver(nEigenpairs, beta, sortingRule, usingSLEPC);

  modal_an.SetSolver(*solver);

  std::streamsize ss = std::cout.precision();
  
  std::cout.precision(std::numeric_limits<STATE>::max_digits10);
  STATE rel_error{0};
  constexpr STATE tol = std::numeric_limits<STATE>::epsilon() * 10000;

  bool computeVectors{false};
  int nit = 0;
  do{
    CSTATE old_beta = beta;
    std::cout<<std::fixed<<"beta: "<<beta.real()<<" "<<beta.imag()<<std::endl;
    modal_an.SetBeta(beta);
    modal_an.Run(computeVectors);
    beta = std::sqrt(modal_an.GetEigenvalues()[0]);
    const auto beta_abs = std::abs(beta);
    const auto old_beta_abs = std::abs(old_beta);
    rel_error = std::abs(beta_abs-old_beta_abs)/std::max(beta_abs,old_beta_abs);
    nit++;
  }while(rel_error > tol);
  std::setprecision(ss);
  std::cout<<"finished in "<<nit<<" iterations"<<std::endl;
  computeVectors = true;
  modal_an.Run(computeVectors);


  {
    TPZSimpleTimer postProc("Post processing");
  
    if(exportVtk){
      const std::string modal_file = prefix+"_modal";
      modal_an.PostProcess(modal_file, vtkRes);
    }
  }

  beta = modal_an.GetEigenvalues()[0];
  modal_an.LoadSolution(0);

  /*********************
   * cmesh(scattering) *
   *********************/

  
  auto scatt_cmesh = [&gmesh,&gmshmats, &modal_cmesh, alphaxPML, alphayPML](){
    //material that will represent our source
    constexpr auto srcMat{"gamma_1"};
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData scatt_data;
    std::map<std::string, std::pair<CSTATE, CSTATE>> scatt_mats;
    scatt_mats["subst_1"] = std::make_pair<CSTATE, CSTATE>(1., 1.);
    scatt_mats["strip_1"] = std::make_pair<CSTATE, CSTATE>(1.55 * 1.55, 1.);
    scatt_mats["subst_2"] = std::make_pair<CSTATE, CSTATE>(1., 1.);
    scatt_mats["strip_2"] = std::make_pair<CSTATE, CSTATE>(1.55 * 1.55, 1.);
    std::map<std::string, wgma::bc::type> scatt_bcs;
    scatt_bcs["modal_bound"] = wgma::bc::type::PEC;
    scatt_bcs["scatt_bound"] = wgma::bc::type::PEC;

    
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, scatt_mats, scatt_bcs,
                                            alphaxPML,alphayPML, scatt_data);
    //get id of source material
    const auto srcId = gmshmats[1].at(srcMat);
    wgma::scattering::SourceWgma src;
    src.id = {srcId};
    src.modal_cmesh = modal_cmesh;
    return wgma::scattering::CMeshScattering2D(gmesh, mode, pOrder, scatt_data,src,
                                               lambda,scale);
  }();

  wgma::scattering::SetPropagationConstant(scatt_cmesh, beta);
  /*********************
   * solve(scattering) *  
   *********************/  


  auto scatt_an = wgma::scattering::Analysis(scatt_cmesh, nThreads,
                                       optimizeBandwidth, filterBoundaryEqs);

  scatt_an.Run();
  //otherwise it will crash on destructor
  wgma::cmeshtools::RemovePeriodicity(modal_cmesh);

  const std::string scatt_file = prefix+"_scatt";
  scatt_an.PostProcess(scatt_file, vtkRes);
  
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
  const auto factor = 1e-6 / scale;
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

#include <slepcepshandler.hpp>
#include <TPZKrylovEigenSolver.h>
//utility functions
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