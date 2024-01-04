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


//!minimum shared sim data
struct SimData{
  //!.msh mesh
  std::string meshfile;
  //!wavelength
  STATE lambda{0.741};
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
  //!Whether the rib is made of copper (instead of photoresist)
  bool rib_copper{false};
  //!mode currently analysed
  wgma::planarwg::mode mode{wgma::planarwg::mode::TM};
  //!polynomial order
  int porder{-1};
  //!whether to filter dirichlet eqs
  bool filterBoundEqs{true};
  //!renumber equations
  bool optimizeBandwidth{true};
  //!output geometric mesh in .txt and .vtk files
  bool printGmesh{false};
  //!post process fields
  bool exportVtk{false};
  //!whether to compute coupling mat
  bool couplingmat{false};
  //!vtk resolution
  int vtkRes{0};
  //!number of threads
  int nThreads{(int)std::thread::hardware_concurrency()};
  //!prefix for both meshes and output files
  std::string prefix{""};
};

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
  const int nEigenpairs,
  const TPZEigenSort sortingRule,
  bool usingSLEPC,
  const std::string &name);

// void SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
//                      wgma::wganalysis::WgmaPlanar &src_an,
//                      wgma::wganalysis::WgmaPlanar &match_an,
//                      const TPZVec<std::map<std::string, int>> &gmshmats,
//                      const SimData &simdata);


SimData GetSimData()
{

  // path for output files
  const std::string path {"res_meta_surf_1d/"};
  // common prefix for both meshes and output files
  const std::string basisName{"meta_surf_1d"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  
  SimData data;
  data.meshfile="meshes/meta_surf_1d.msh";
  data.lambda = 0.741;
  /*
    Given the small dimensions of the domain, scaling it can help in
    achieving good precision. Using 1./k0 as a scale factor results in
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  data.scale = data.lambda/(2*M_PI);
  data.mode = wgma::planarwg::mode::TM;
  data.rib_copper = false;
  data.n_copper = 0.1;
  data.k_copper = 7;
  data.n_az = 1.622;
  data.k_az = 0;
  data.n_air = 1;
  data.n_rib = data.rib_copper ? data.n_copper : data.n_az;
  data.k_rib = data.rib_copper ? data.k_copper : data.k_az;
  data.porder = 4;
  data.filterBoundEqs = true;
  data.printGmesh=true;
  data.exportVtk = true;
  data.couplingmat = true;
  data.vtkRes=0;
  data.prefix = prefix;
  return std::move(data);
}

int main(int argc, char *argv[]) {

#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif
  const SimData simdata=GetSimData();

  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(simdata.prefix));
  
  /******************
   * eigensolver options *
   ******************/

  // how to sort eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  constexpr bool usingSLEPC {true};
  constexpr int nEigenpairs_top{10};
  constexpr int nEigenpairs_bottom{10};
  const CSTATE target_top{simdata.n_air*simdata.n_air*1.0000001};
  const CSTATE target_bottom{simdata.n_copper*simdata.n_copper*1.00000001};

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
  if (simdata.printGmesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = simdata.prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }


  TPZAutoPointer<wgma::wganalysis::WgmaPlanar> modal_bottom_an, modal_top_an;
  
  {
    TPZSimpleTimer timer("Modal analysis",true);
    

    /********************************
   * cmesh(modal analysis):bottom   *
   ********************************/
    const STATE &n_copper = simdata.n_copper;
    const STATE &k_copper = simdata.k_copper;
    //complex permittivity
    const CSTATE epsilon_copper =
      n_copper*n_copper-k_copper*k_copper + 2i*n_copper*k_copper;
    modal_bottom_an = ComputeModalAnalysis(gmesh, gmshmats,
                                           periodic_els,
                                           epsilon_copper,simdata,
                                           target_bottom, nEigenpairs_bottom,
                                           sortingRule, usingSLEPC,
                                           "bottom");
  
    /********************************
     * cmesh(modal analysis):top   *
     ********************************/

    const STATE &n_air = simdata.n_air;
    //complex permittivity
    const CSTATE epsilon_air = n_air*n_air;
    modal_top_an = ComputeModalAnalysis(gmesh, gmshmats,
                                        periodic_els,
                                        epsilon_air,simdata,
                                        target_top, nEigenpairs_top,
                                        sortingRule, usingSLEPC,
                                        "top");
  }
  // SolveScattering(gmesh, modal_bottom_an,
  //                 modal_top_an,gmshmats,simdata);
  //otherwise it will crash on destructor
  wgma::cmeshtools::RemovePeriodicity(modal_top_an->GetMesh());
  wgma::cmeshtools::RemovePeriodicity(modal_bottom_an->GetMesh());
  return 0;
}

#include <TPZKrylovEigenSolver.h>

TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const CSTATE target, const int nEigen,
            TPZEigenSort sorting, bool &usingSLEPC);

void ComputeModes(wgma::wganalysis::WgmaPlanar &an,
                  const int nThreads);

void ComputeCouplingMat(wgma::wganalysis::WgmaPlanar &an,
                        std::string filename);

void PostProcessModes(wgma::wganalysis::WgmaPlanar &an,
                      std::string filename,
                      const int vtkres,
                      std::set<int> sols = {});

TPZAutoPointer<wgma::wganalysis::WgmaPlanar>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const std::map<int64_t,int64_t> &periodic_els,
  const CSTATE epsilon_mat,
  const SimData& simdata,
  const CSTATE target,
  const int nEigenpairs,
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
    const auto &lambda = simdata.lambda;
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
                                                 lambda, scale);
  }(epsilon_mat);
  /******************************
   * solve(modal analysis left) *
   ******************************/

  TPZAutoPointer<wgma::wganalysis::WgmaPlanar>
    modal_an =
    new wgma::wganalysis::WgmaPlanar(modal_cmesh, simdata.nThreads,
                                     simdata.optimizeBandwidth,
                                     simdata.filterBoundEqs);
  {
    auto solver = SetupSolver(target, nEigenpairs, sortingRule, usingSLEPC);
    modal_an->SetSolver(*solver);
  }

  const std::string modal_file{simdata.prefix+"_modal_"+name};
  ComputeModes(*modal_an, simdata.nThreads);
  if(simdata.couplingmat){
    std::string couplingfile{simdata.prefix+"_coupling_"+name};
    ComputeCouplingMat(*modal_an,couplingfile);
  }
  if(simdata.exportVtk){
    PostProcessModes(*modal_an, modal_file, simdata.vtkRes);
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
                  const int nThreads)
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
  //                                                               conj,nThreads).Normalise();
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