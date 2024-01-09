/**
meta_surf.cpp

This target performs the modal analysis of a two 1d lines
and use them to compute the waveguide port bc for the scenario of
light illuminating a meta surface
***/

// wgma includes
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include <TPZEigenAnalysis.h>
#include <wganalysis.hpp>
#include <scattering.hpp>
#include <util.hpp>
#include <slepcepshandler.hpp>
#include <post/solutionnorm.hpp>
#include <post/waveguideportbc.hpp>
#include <post/waveguidecoupling.hpp>
#include <post/orthowgsol.hpp>
// pz includes
#include <MMeshType.h>      //for MMeshType
#include <TPZSimpleTimer.h> //for TPZSimpleTimer
#include <pzlog.h>          //for TPZLogger
#include <TPZMaterial.h>
#include <TPZVTKGenerator.h>
#include <pzinterpolationspace.h>

#include <regex>//for string search
#include <thread>

//1i
using namespace std::complex_literals;

//!minimum shared sim data
struct SimData{
  //!.msh mesh
  std::string meshfile;
  //!wavelength
  STATE wavelength{0.741};
  //!geometric scaling (floating point precision)
  REAL scale{1};
  //!real part of refractive index of copper
  STATE n_copper{1};
  //!imag part of refractive index of copper
  STATE k_copper{1};
  //!real part of refractive index of photoresist
  STATE n_az{1};
  //!imag part of refractive index of photoresist
  STATE k_az{1};
  //!real part of refractive index of rib
  STATE n_rib{1};
  //!imag part of refractive index of rib
  STATE k_rib{1};
  //!refractive index of air 
  STATE n_air{1};
  //!mode currently analysed
  wgma::planarwg::mode mode{wgma::planarwg::mode::TM};
  //!polynomial order
  int porder{-1};
  //!number eigenpairs top
  int n_eigen_top{300};
  //!number eigenpairs bottom (0 for no modal analysis)
  int n_eigen_bot{300};
  //!number of modes of wgbc
  std::vector<int> nmodes = {1,10,15,100,200,250};
  //!target eigenvalue top
  CSTATE target_top{n_air*n_air*1.0001};
  //!target eigenvalue bottom
  CSTATE target_bot{n_copper*n_copper*1.0001};
  //!whether to filter dirichlet eqs
  bool filter_bnd_eqs{true};
  //!renumber equations
  bool optimize_bandwidth{true};
  //!output geometric mesh in .txt and .vtk files
  bool print_gmesh{false};
  //!post process modal fields
  bool export_vtk_modes{false};
  //!post process scatt fields
  bool export_vtk_scatt{false};
  /**
     @brief Computes reflection and exports to .csv file
     @note This will APPEND to a given csv file, be sure to
     erase it in between non-related runs
   */
  bool compute_reflection_norm{true};
  //!whether to compute coupling mat
  bool couplingmat{false};
  //!vtk resolution
  int vtk_res{0};
  //!initial count for vtk files
  int initial_count{0};
  //!number of threads
  int n_threads{(int)std::thread::hardware_concurrency()};
  //!prefix for both meshes and output files
  std::string prefix{""};
};

//! Reads sim data from file
SimData ReadSimData(const std::string &dataname);

//!needed data from modal analysis to create waveguide port bc
struct WgbcData{
  TPZAutoPointer<TPZCompMesh> cmesh;
  TPZFMatrix<CSTATE> wgbc_k;
  TPZVec<CSTATE> wgbc_f;
};

TPZAutoPointer<wgma::wganalysis::WgmaPlanar>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const std::map<int64_t,int64_t> &periodic_els,
  const CSTATE epsilon_mat,
  const SimData& simdata,
  const CSTATE target,
  int &nEigenpairs,
  const TPZEigenSort sortingRule,
  bool usingSLEPC,
  const std::string &name);

void SolveScattering(TPZAutoPointer<TPZGeoMesh>gmesh,
                     TPZAutoPointer<wgma::wganalysis::WgmaPlanar> src_an,
                     TPZAutoPointer<wgma::wganalysis::WgmaPlanar> match_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const std::map<int64_t,int64_t> &periodic_els,
                     const CSTATE epsilon_rib,
                     const CSTATE epsilon_copper,
                     const CSTATE epsilon_air,
                     const SimData& simdata);


SimData GetSimData(bool rib_copper)
{

  // path for output files
  const std::string path {"res_meta_surf_1d/"};
  // common prefix for both meshes and output files
  const std::string basisName{"meta_surf_1d"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  
  SimData data;
  data.meshfile="meshes/meta_surf_1d.msh";
  data.wavelength = rib_copper ? 0.746 : 0.741;
  /*
    Given the small dimensions of the domain, scaling it can help in
    achieving good precision. Using 1./k0 as a scale factor results in
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  data.scale = data.wavelength/(2*M_PI);
  data.mode = wgma::planarwg::mode::TM;
  data.n_copper = 0.1;
  data.k_copper = 7;
  data.n_az = 1.622;
  data.k_az = 0;
  data.n_air = 1;
  data.n_rib = rib_copper ? data.n_copper : data.n_az;
  data.k_rib = rib_copper ? data.k_copper : data.k_az;
  data.porder = 4;
  data.filter_bnd_eqs = true;
  data.print_gmesh=true;
  data.export_vtk_modes = false;
  data.compute_reflection_norm = true;
  data.export_vtk_scatt = true;
  data.couplingmat = false;
  data.n_eigen_top=300;
  data.n_eigen_bot=300;
  data.nmodes = {1,10,15,100,200,250};
  data.target_top=data.n_air*data.n_air*1.0001;
  data.target_bot=data.n_copper*data.n_copper*1.0001;
  data.vtk_res=0;
  data.prefix = prefix;
  return std::move(data);
}

int main(int argc, char *argv[]) {

#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif
  if(argc>2){
    PZError<<"Unexpected number of parameters. USAGE: ./meta_surf_1d param_file\n"
           <<"or ./meta_surf_1d to run with hardcoded params"<<std::endl;
    return -1;
  }
  
  SimData simdata= [argc, &argv](){
    constexpr bool rib_copper{true};
    if(argc==1){
      auto sd = GetSimData(rib_copper);
      sd.prefix+=rib_copper ? "_cu" : "_az";
      return sd;
    }
    const std::string dataname = argv[1];
    return ReadSimData(dataname);
  }();

  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(simdata.prefix));
  
  /******************
   * eigensolver options *
   ******************/

  // how to sort eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  constexpr bool usingSLEPC {true};

  /*********
   * begin *
   *********/




  /*************
   * geometry  *
   ************/
  // scoped-timer
  TPZSimpleTimer total("Total");

  TPZVec<std::map<std::string, int>> gmshmats;
  constexpr bool verbosity_lvl{false};
  std::map<int64_t,int64_t> periodic_els;
  auto gmesh = wgma::gmeshtools::ReadPeriodicGmshMesh(simdata.meshfile,
                                                      simdata.scale,
                                                      gmshmats,periodic_els,
                                                      verbosity_lvl);

  //first we rotate the mesh

  wgma::gmeshtools::RotateMesh(gmesh,{0,0,1},M_PI/2);
  // print wgma_gmesh to .txt and .vtk format
  if (simdata.print_gmesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = simdata.prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }

  const STATE &n_air = simdata.n_air;
  //complex permittivity
  const CSTATE epsilon_air = n_air*n_air;
  const STATE &n_copper = simdata.n_copper;
  const STATE &k_copper = simdata.k_copper;
  //complex permittivity
  const CSTATE epsilon_copper =
    n_copper*n_copper-k_copper*k_copper + 2i*n_copper*k_copper;

  const STATE &n_rib = simdata.n_rib;
  const STATE &k_rib = simdata.k_rib;
  //complex permittivity
  const CSTATE epsilon_rib =
    n_rib*n_rib-k_rib*k_rib + 2i*n_rib*k_rib;

  TPZAutoPointer<wgma::wganalysis::WgmaPlanar> modal_bottom_an, modal_top_an;
  
  {
    TPZSimpleTimer timer("Modal analysis",true);
    
    if(simdata.n_eigen_bot){
      //maybe there is no need to perform this
      /********************************
       * cmesh(modal analysis):bottom   *
       ********************************/
      modal_bottom_an = ComputeModalAnalysis(gmesh, gmshmats,
                                             periodic_els,
                                             epsilon_copper,simdata,
                                             simdata.target_bot, simdata.n_eigen_bot,
                                             sortingRule, usingSLEPC,
                                             "bottom");
    }
    /********************************
     * cmesh(modal analysis):top   *
     ********************************/
    modal_top_an = ComputeModalAnalysis(gmesh, gmshmats,
                                        periodic_els,
                                        epsilon_air,simdata,
                                        simdata.target_top, simdata.n_eigen_top,
                                        sortingRule, usingSLEPC,
                                        "top");
  }
  SolveScattering(gmesh, modal_top_an,
                  modal_bottom_an,gmshmats,
                  periodic_els,
                  epsilon_rib,
                  epsilon_copper,
                  epsilon_air,
                  simdata);
  //otherwise it will crash on destructor
  wgma::cmeshtools::RemovePeriodicity(modal_top_an->GetMesh());
  if(modal_bottom_an){wgma::cmeshtools::RemovePeriodicity(modal_bottom_an->GetMesh());}
  return 0;
}

#include <json_util.hpp>
#include <TPZKrylovEigenSolver.h>
#include <TPZPardisoSolver.h>

TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const CSTATE target, const int nEigen,
            TPZEigenSort sorting, bool &usingSLEPC);

void ComputeModes(wgma::wganalysis::WgmaPlanar &an,
                  const REAL scale,
                  const int n_threads);

void ComputeCouplingMat(wgma::wganalysis::WgmaPlanar &an,
                        std::string filename);

void PostProcessModes(wgma::wganalysis::WgmaPlanar &an,
                      std::string filename,
                      const int vtkres,
                      std::set<int> sols = {});
void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WgbcData& src_data,
                          WgbcData& match_data,
                          const int nmodes_src,
                          const int nmodes_match,
                          const SimData &simdata);

void ComputeWgbcCoeffs(wgma::wganalysis::WgmaPlanar& an,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       wgma::planarwg::mode mode,
                       const bool positive_z, const TPZVec<CSTATE> &coeff);

void AddWaveguidePortContribution(wgma::scattering::Analysis &scatt_an, 
                                  const int64_t indep_con_id,
                                  const int nm,
                                  const TPZFMatrix<CSTATE> &wgbc_k,
                                  const TPZVec<CSTATE> &wgbc_f);


//! Reads sim data from file
SimData ReadSimData(const std::string &dataname)
{
  using json = nlohmann::json;
  std::ifstream f(dataname);
  json data = json::parse(f);
  SimData sd;

  sd.meshfile = data["meshfile"];
  //!wavelength
  sd.wavelength=data["wavelength"];
  //!geometric scaling (floating point precision)
  sd.scale=data["scale"];
  //!real part of refractive index of copper
  sd.n_copper=data["n_copper"];
  //!imag part of refractive index of copper
  sd.k_copper=data["k_copper"];
  //!real part of refractive index of photoresist
  sd.n_az=data["n_az"];
  //!imag part of refractive index of photoresist
  sd.k_az=data["k_az"];
  //!real part of refractive index of rib
  sd.n_rib=data["n_rib"];
  //!imag part of refractive index of rib
  sd.k_rib=data["k_rib"];
  //!refractive index of air 
  sd.n_air=data["n_air"];
  //!mode currently analysed
  sd.mode=wgma::planarwg::string_to_mode(data["mode"]);
  //!polynomial order
  sd.porder=data["porder"];
  //!number eigenpairs top
  sd.n_eigen_top=data["n_eigen_top"];
  //!number eigenpairs bottom
  sd.n_eigen_bot=data["n_eigen_bot"];
  //!vector with number of modes of wgbc
  sd.nmodes = data["nmodes"].get<std::vector<int>>();
  //!target eigenvalue top
  data.at("target_top").get_to(sd.target_top);
  //!target eigenvalue bottom
  data.at("target_bot").get_to(sd.target_bot);
  //!whether to filter dirichlet eqs
  sd.filter_bnd_eqs=data["filter_bnd_eqs"];
  //!renumber equations
  sd.optimize_bandwidth=true;
  //!output geometric mesh in .txt and .vtk files
  sd.print_gmesh=data.value("print_gmesh",false);
  //!post process modal fields
  sd.export_vtk_modes=data.value("export_vtk_modes",false);
  //!export reflection norm
  sd.compute_reflection_norm = data.value("compute_reflection_norm",false);
  //!post process scatt fields
  sd.export_vtk_scatt=data.value("export_vtk_scatt",false);
  //!whether to compute coupling mat
  sd.couplingmat=data.value("couplingmat",false);
  //!vtk resolution
  sd.vtk_res=data.value("vtk_res",(int)0);
  //!vtk initial count
  sd.initial_count=data.value("initial_count",(int)0);
  //!prefix for both meshes and output files
  sd.prefix=data["prefix"];
  //!number of threads
  sd.n_threads = data.value("n_threads",(int)std::thread::hardware_concurrency());
  return sd;
}

TPZAutoPointer<wgma::wganalysis::WgmaPlanar>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const std::map<int64_t,int64_t> &periodic_els,
  const CSTATE epsilon_mat,
  const SimData& simdata,
  const CSTATE target,
  int &nEigenpairs,
  const TPZEigenSort sortingRule,
  bool usingSLEPC,
  const std::string &name)
{
  auto modal_cmesh = [gmesh,&gmshmats,
                      &simdata,
                      &periodic_els,
                      &name](CSTATE epsilon){
    // setting up cmesh data
    const auto &mode = simdata.mode;
    const auto &pOrder = simdata.porder;
    const auto &wavelength = simdata.wavelength;
    const auto &scale = simdata.scale;
    
    wgma::cmeshtools::PhysicalData modal_data;
    std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
    modal_mats["bnd_"+name] = std::pair<CSTATE, CSTATE>(epsilon, 1.);
    std::map<std::string, wgma::bc::type> modal_bcs;
    modal_bcs["bnd_wgma_"+name+"_0"] = wgma::bc::type::PERIODIC;
    modal_bcs["bnd_wgma_"+name+"_1"] = wgma::bc::type::PERIODIC;
    //dimension of the modal analysis 
    constexpr int modal_dim{1};
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                            {0}, modal_data, modal_dim);
    return wgma::wganalysis::CMeshWgma1DPeriodic(gmesh,mode,pOrder,modal_data,
                                                 periodic_els,
                                                 wavelength, scale);
  }(epsilon_mat);
  /******************************
   * solve(modal analysis left) *
   ******************************/

  TPZAutoPointer<wgma::wganalysis::WgmaPlanar>
    modal_an =
    new wgma::wganalysis::WgmaPlanar(modal_cmesh, simdata.n_threads,
                                     simdata.optimize_bandwidth,
                                     simdata.filter_bnd_eqs);

  //now we check if neigenpairs > neq
  const int neq = modal_an->StructMatrix()->EquationFilter().NActiveEquations();
  nEigenpairs = nEigenpairs >= neq ? neq : nEigenpairs;
  {
    auto solver = SetupSolver(target, nEigenpairs, sortingRule, usingSLEPC);
    modal_an->SetSolver(*solver);
  }

  const std::string modal_file{simdata.prefix+"_modal_"+name};
  ComputeModes(*modal_an, simdata.scale, simdata.n_threads);
  if(simdata.couplingmat){
    std::string couplingfile{simdata.prefix+"_coupling_"+name};
    ComputeCouplingMat(*modal_an,couplingfile);
  }
  if(simdata.export_vtk_modes){
    PostProcessModes(*modal_an, modal_file, simdata.vtk_res);
  }
  return modal_an;
}

//utility functions
TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const CSTATE target,const int neigenpairs,
            TPZEigenSort sorting, bool &usingSLEPC)
{

#ifndef WGMA_USING_SLEPC
  if(usingSLEPC){
    std::cout<<"wgma was not configured with slepc. defaulting to: "
             <<"TPZKrylovSolver"<<std::endl;
    usingSLEPC = false;
  }
#endif

  TPZAutoPointer<TPZEigenSolver<CSTATE>> solver{nullptr};
  const int krylovDim{5*neigenpairs};
  if (usingSLEPC){
    using namespace wgma::slepc;
    /*
      The following are suggested SLEPc settings.
      NOTE: -1 stands for PETSC_DECIDE
    */
    
    constexpr STATE eps_tol = 1e-18;//PETSC_DECIDE
    constexpr int eps_max_its = -1;//PETSC_DECIDE
    constexpr EPSConv eps_conv_test = EPSConv::EPS_CONV_REL;
    constexpr EPSWhich eps_which = EPSWhich::EPS_TARGET_REAL;
    
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


void ComputeModes(wgma::wganalysis::WgmaPlanar &an,
                  const REAL scale,
                  const int n_threads)
{
  
  TPZSimpleTimer analysis("Modal analysis");
  static constexpr bool computeVectors{true};
  
  an.Run(computeVectors);

  an.LoadAllSolutions();
  constexpr bool conj{false};
  TPZSimpleTimer timer("Ortho",true);
  constexpr STATE tol{1e-9};
  const int n_ortho = wgma::post::OrthoWgSol(an,tol,conj);
  std::cout<<"orthogonalised  "<<n_ortho<<" eigenvectors"<<std::endl;
  //now we normalise them in case we need to compute the reflective spectra
  {
    auto cmesh = an.GetMesh();
    //leave empty for all valid matids
    std::set<int> matids {};
    constexpr bool conj{true};
    wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(cmesh,matids,
                                                                conj,n_threads).Normalise();
    TPZFMatrix<CSTATE> &mesh_sol=cmesh->Solution();

    const int sz = mesh_sol.Rows() * mesh_sol.Cols();
    auto *sol_ptr = mesh_sol.Elem();
    for(int i = 0; i < sz; i++){
      *sol_ptr++= *sol_ptr*scale;
    }
    //we update analysis object
    an.SetEigenvectors(mesh_sol);
    an.LoadAllSolutions();
  }
  // constexpr bool ortho{true};
  // if(ortho){
  //   TPZSimpleTimer timer("Ortho",true);
  //   constexpr STATE tol{1e-9};
  //   const int n_ortho = wgma::post::OrthoWgSol(an,tol,conj);
  //   std::cout<<"orthogonalised  "<<n_ortho<<" eigenvectors"<<std::endl;
  // }else{
  //   auto cmesh = an.GetMesh();
  //   //leave empty for all valid matids
  //   std::set<int> matids {};
  //   wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(cmesh,matids,
  //                                                               conj,n_threads).Normalise();
  //   TPZFMatrix<CSTATE> &mesh_sol=cmesh->Solution();
  //   //we update analysis object
  //   an.SetEigenvectors(mesh_sol);
  //   an.LoadAllSolutions();
  // }
}
void ComputeCouplingMat(wgma::wganalysis::WgmaPlanar &an,
                        std::string filename)
{
  
  using namespace wgma::post;

  std::set<int> matids;
  constexpr bool conj{false};
  const int nthreads = std::thread::hardware_concurrency();
  WaveguideCoupling<SingleSpaceIntegrator> integrator(an.GetMesh(),
                                                       matids,
                                                       conj,
                                                       nthreads
                                                       );
  auto &an_beta = an.GetEigenvalues();
  const int nsol = an_beta.size();
  TPZVec<CSTATE> betavec(nsol,0);
  for(int i = 0; i < nsol; i++){
    betavec[i] = sqrt(-an_beta[i]);
  }
  integrator.SetBeta(betavec);
  
  integrator.ComputeCoupling();
  TPZFMatrix<CSTATE> couplingmat;
  integrator.GetCoupling(couplingmat);
  std::ofstream matfile(filename);
  couplingmat.Print("",matfile,ECSV);
}

void PostProcessModes(wgma::wganalysis::WgmaPlanar &an,
                      std::string filename,
                      const int vtkres,
                      std::set<int> sols){
  TPZSimpleTimer postProc("Post processing");
    
  const std::string file = filename;
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_abs",
  };
  auto cmesh = an.GetMesh();
  auto vtk = TPZVTKGenerator(cmesh, fvars, file, vtkres);
  const int nthreads = an.StructMatrix()->GetNumThreads();
  vtk.SetNThreads(nthreads);

  if(sols.size() == 0){
    const auto nsol = std::min((int64_t)20,an.GetEigenvectors().Cols());
    for(auto is = 0; is < nsol; is++){
      sols.insert(is);
    }
  }
  
  std::cout<<"Exporting "<<sols.size()<<" solutions"<<std::endl;
  for(auto isol : sols){
    an.LoadSolution(isol);
    vtk.Do();
  }
}

void SolveScattering(TPZAutoPointer<TPZGeoMesh>gmesh,
                     TPZAutoPointer<wgma::wganalysis::WgmaPlanar> src_an,
                     TPZAutoPointer<wgma::wganalysis::WgmaPlanar> match_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const std::map<int64_t,int64_t> &periodic_els,
                     const CSTATE epsilon_rib,
                     const CSTATE epsilon_copper,
                     const CSTATE epsilon_air,
                     const SimData& simdata)
{

  //maybe bottom is just PEC
  const bool has_bot = match_an;
  auto scatt_mesh = [gmesh,&gmshmats,&simdata, &periodic_els, &has_bot]
    (CSTATE e_rib, CSTATE e_copper, CSTATE e_air){
    const auto mode = simdata.mode;
    const auto &pOrder = simdata.porder;
    const auto &wavelength = simdata.wavelength;
    const auto &scale = simdata.scale;
    
    const std::string &prefix = simdata.prefix;

    
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData scatt_data;
    std::map<std::string, std::pair<CSTATE, CSTATE>> scatt_mats;
    scatt_mats["rib"] = std::pair<CSTATE, CSTATE>(e_rib, 1.);
    scatt_mats["copper"] = std::pair<CSTATE, CSTATE>(e_copper, 1.);
    scatt_mats["air"] = std::pair<CSTATE, CSTATE>(e_air,1.);
    std::map<std::string, wgma::bc::type> scatt_bcs;
    scatt_bcs["bnd_periodic_1"] = wgma::bc::type::PERIODIC;
    scatt_bcs["bnd_periodic_2"] = wgma::bc::type::PERIODIC;
    if(!has_bot){scatt_bcs["bnd_bottom"] = wgma::bc::type::PEC;}
    
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, scatt_mats, scatt_bcs,
                                            {0,0}, scatt_data);

    


    /*
      probe mats are regions of the domain in which we want to be able
      to evaluate our solution
      they are also used to ensure that the computational elements are created
      so every region that will be used in a waveguide port bc must be
      also inserted as a probe mat
    */
    {
      std::vector<std::string> probeMats;
      probeMats.push_back("bnd_top");
      if(has_bot){probeMats.push_back("bnd_bottom");}
      for(auto &mat : probeMats){
        const auto matdim = 1;
        const auto id = gmshmats[matdim].at(mat);
        scatt_data.probevec.push_back({id,matdim});
      }
    }
    //empty
    std::set<int> src_ids;
    return wgma::scattering::CMeshScattering2DPeriodic(gmesh, mode, pOrder,
                                                       scatt_data,periodic_els,
                                                       src_ids,wavelength,scale);
  }(epsilon_rib,epsilon_copper,epsilon_air);
  /*********************
   * solve(scattering) *  
   *********************/  
  TPZSimpleTimer tscatt("Scattering");


  const auto n_eigenpairs_top = src_an->GetEigenvalues().size();
  /*
    the source is written as a linear combination of the modes
    this vector contains the coefficients of such combination
   */
  TPZVec<CSTATE> src_coeffs(n_eigenpairs_top,0);
  src_coeffs[0] = 1;

  //index of the number of modes to be used to restrict the dofs on waveguide bcs
  auto nmodes = simdata.nmodes;
  
  if(has_bot){match_an->LoadAllSolutions();}

  /*
    this mesh will be used to compute the reflection
    so, first, we need to store the solution corresponding to our src
   */
  TPZFMatrix<CSTATE> src_sol;
  TPZAutoPointer<TPZCompMesh> ref_mesh = nullptr;

  if(simdata.compute_reflection_norm){
    ref_mesh = src_an->GetMesh()->Clone();
    TPZFMatrix<CSTATE> &sol = src_an->GetMesh()->Solution();
    const int64_t neq_ref_mesh = sol.Rows();
    src_sol.Redim(neq_ref_mesh,1);
    
    for(int i = 0; i < n_eigenpairs_top; i++){
      if(src_coeffs[i]!=0.){
        src_an->LoadSolution(i);
        TPZFMatrix<CSTATE> sol = src_an->GetMesh()->Solution();
        sol*=src_coeffs[i];
        src_sol += sol;
      }
    }
  }
  
  //now we solve varying the number of modes used in the wgbc
  src_an->LoadAllSolutions();
  
  const auto mode = simdata.mode;
  //compute wgbc coefficients
  WgbcData src_data;
  src_data.cmesh = src_an->GetMesh();
  ComputeWgbcCoeffs(*src_an, src_data.wgbc_k, src_data.wgbc_f, mode, false, src_coeffs);

  WgbcData match_data;
  if(has_bot){
    match_data.cmesh = match_an->GetMesh();
    ComputeWgbcCoeffs(*match_an, match_data.wgbc_k, match_data.wgbc_f,mode, true, {});
  }
  constexpr bool print_wgbc_mats{false};
  if(print_wgbc_mats){
    std::ofstream matfile{simdata.prefix+"_wgbc_k.csv"};
    src_data.wgbc_k.Print("",matfile,ECSV);
    const int nsol = src_data.wgbc_f.size();
    TPZFMatrix<CSTATE> wgbc_f(nsol,1,src_data.wgbc_f.begin(),nsol);
    std::ofstream vecfile{simdata.prefix+"_wgbc_f.csv"};
    wgbc_f.Print("",vecfile,ECSV);
  }

  //set up post processing
  const std::string scatt_file = simdata.prefix+"_scatt";
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_imag",
    "Field_abs"};
  TPZAutoPointer<TPZVTKGenerator> vtk = nullptr;
  if(simdata.export_vtk_scatt){
    vtk = new TPZVTKGenerator(scatt_mesh, fvars, scatt_file, simdata.vtk_res);
    vtk->SetNThreads(simdata.n_threads);
    const int count = simdata.initial_count;
    vtk->SetStep(count);
  }
    
  TPZFMatrix<CSTATE> ref_sol = src_sol;
  for(int im = 0; im < nmodes.size(); im++){
    const int nm = nmodes[im];
    if(!nm){continue;}
    if(im != 0){
      wgma::cmeshtools::SetPeriodic(scatt_mesh, periodic_els);
    }
    //now we check how many modes we have in each
    const int neq_src = src_an->StructMatrix()->EquationFilter().NActiveEquations();
    const int nmodes_src = neq_src > nm ? nm : neq_src;
    const int neq_match = has_bot ?
      (match_an->StructMatrix()->EquationFilter().NActiveEquations())
      :
      0;
    const int nmodes_match = neq_match > nm ? nm : neq_match;
    RestrictDofsAndSolve(scatt_mesh, src_data, match_data, nmodes_src,nmodes_match, simdata);

    
    if(simdata.compute_reflection_norm){
      //now we compute the reflection
      wgma::cmeshtools::ExtractSolFromMesh(ref_mesh, scatt_mesh, ref_sol);
      //subtract source
      ref_sol -= src_sol;
      ref_mesh->LoadSolution(ref_sol);
      wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator> sol_norm(ref_mesh);
      sol_norm.SetNThreads(simdata.n_threads);
      const STATE ref_norm = std::real(sol_norm.ComputeNorm()[0]);
      std::string outputfile = simdata.prefix+"_reflection.csv";
      std::ofstream ost;
      ost.open(outputfile, std::ios_base::app);
      ost << simdata.wavelength<<','<<ref_norm<<std::endl;
    }
    //plot
    if(simdata.export_vtk_scatt){vtk->Do();}
    //removing restrictions
    wgma::cmeshtools::RemovePeriodicity(scatt_mesh);
  }
  if(simdata.compute_reflection_norm){wgma::cmeshtools::RemovePeriodicity(ref_mesh);}
}

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WgbcData& src_data,
                          WgbcData& match_data,
                          const int nmodes_src,
                          const int nmodes_match,
                          const SimData &simdata)
{

  auto gmesh = scatt_mesh->Reference();
  

  TPZAutoPointer<TPZCompMesh> match_mesh = match_data.cmesh;

  if((match_mesh && nmodes_match == 0)||(!match_mesh && nmodes_match != 0)){
    DebugStop();
  }
   
  auto src_mesh = src_data.cmesh;

  /**
     no dirichlet connects must be restricted!!
   **/
  std::set<int64_t> boundConnects;
  wgma::cmeshtools::FindDirichletConnects(scatt_mesh, boundConnects);
  int64_t indep_con_id_match = -1;
  if(match_mesh){
    indep_con_id_match = wgma::cmeshtools::RestrictDofs(scatt_mesh, match_mesh,
                                                        nmodes_match, boundConnects);
  }
  const int64_t indep_con_id_src =
    wgma::cmeshtools::RestrictDofs(scatt_mesh, src_mesh, nmodes_src, boundConnects);

  constexpr bool sym{false};
  auto scatt_an = wgma::scattering::Analysis(scatt_mesh, simdata.n_threads,
                                             simdata.optimize_bandwidth,
                                             simdata.filter_bnd_eqs,
                                             sym);

  std::cout<<"nmodes on ingoing boundary: "<<nmodes_src<<std::endl;
  std::cout<<"nmodes on outgoing boundary: "<<nmodes_match<<std::endl;

  {
    TPZSimpleTimer timer("Assemble_"+std::to_string(nmodes_src),true);
    std::cout<<"Assembling...";
    scatt_an.Assemble();
    std::cout<<"\nAssembled!"<<std::endl;
  }
  {
    TPZSimpleTimer timer("Add contribution_"+std::to_string(nmodes_src),true);
  //now we must add the waveguide port terms
    if(match_mesh){
      AddWaveguidePortContribution(scatt_an, indep_con_id_match,
                                   nmodes_match, match_data.wgbc_k, match_data.wgbc_f);
    }
    AddWaveguidePortContribution(scatt_an, indep_con_id_src,
                                 nmodes_src, src_data.wgbc_k, src_data.wgbc_f);
  }

  {
    //get pardiso control
    auto *pardiso = scatt_an.GetSolver().GetPardisoControl();
    pardiso->SetMessageLevel(0);
    pardiso->ResetParam();
    constexpr auto sys_type = SymProp::NonSym;
    constexpr auto prop = TPZPardisoSolver<CSTATE>::MProperty::EIndefinite;
    pardiso->SetMatrixType(sys_type,prop);
  }
  
  TPZSimpleTimer tscatt("Solve",true);
  std::cout<<"Solving...";
  scatt_an.Solve();
  std::cout<<"\rSolved!"<<std::endl;
}

void ComputeWgbcCoeffs(wgma::wganalysis::WgmaPlanar& an,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       wgma::planarwg::mode mode,
                       const bool positive_z, const TPZVec<CSTATE> &coeff){
  auto mesh = an.GetMesh();
  
  wgma::post::WaveguidePortBC<wgma::post::SingleSpaceIntegrator> wgbc(mesh);
  TPZManVector<CSTATE,1000> betavec = an.GetEigenvalues();
  for(auto &b : betavec){b = std::sqrt(b);}
  wgbc.SetPositiveZ(positive_z);
  wgbc.SetAdjoint(false);
  if(coeff.size()){
    wgbc.SetSrcCoeff(coeff);
  }
  const bool is_te = mode == wgma::planarwg::mode::TE;
  wgbc.SetTE(is_te);
  wgbc.SetBeta(betavec);
  wgbc.ComputeContribution();
  wgbc.GetContribution(wgbc_k,wgbc_f);
}

void AddWaveguidePortContribution(wgma::scattering::Analysis &scatt_an, 
                                  const int64_t indep_con_id,
                                  const int nm,
                                  const TPZFMatrix<CSTATE> &wgbc_k,
                                  const TPZVec<CSTATE> &wgbc_f)
{
  auto mat = scatt_an.GetSolver().Matrix();
  TPZFMatrix<CSTATE>& fvec = scatt_an.Rhs();
  auto scatt_mesh = scatt_an.GetMesh();
  const auto &indep_con = scatt_mesh->ConnectVec()[indep_con_id];
  const auto &block = scatt_mesh->Block();
  const auto seqnum = indep_con.SequenceNumber();
  const auto pos_orig = block.Position(seqnum);
  const auto sz = block.Size(seqnum);
  if(sz!=nm){DebugStop();}

  int64_t pos_filt{-1};
  {
    TPZManVector<int64_t,1> posvec_orig(1,pos_orig), posvec_filt(1,pos_orig);
    scatt_an.StructMatrix()->EquationFilter().Filter(posvec_orig, posvec_filt);
    pos_filt = posvec_filt[0];
    //they are sequential
  }

  //we add to vec
  {
    CSTATE *ptr_f = &fvec.g(pos_orig,0);
    CSTATE *ptr_wgbc = wgbc_f.begin();
    for(auto imode = 0; imode < nm; imode++){
      *ptr_f++ += *ptr_wgbc++;
    }
  }
  //now we add to matrix
  for(auto imode = 0; imode < nm; imode++){
    const auto ipos = pos_filt+imode;
    for(auto jmode = 0; jmode < nm; jmode++){
      const auto jpos = pos_filt+jmode;
      const auto kval = mat->Get(ipos,jpos);
      mat->Put(ipos,jpos,kval+wgbc_k.Get(imode,jmode));
    }
  }
}
