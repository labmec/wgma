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

void SolveScattering(TPZAutoPointer<TPZGeoMesh>gmesh,
                     TPZAutoPointer<wgma::wganalysis::WgmaPlanar> src_an,
                     TPZAutoPointer<wgma::wganalysis::WgmaPlanar> match_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const std::map<int64_t,int64_t> &periodic_els,
                     const CSTATE epsilon_rib,
                     const CSTATE epsilon_copper,
                     const CSTATE epsilon_air,
                     const SimData& simdata);


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
  data.rib_copper = true;
  data.lambda = data.rib_copper ? 0.746 : 0.741;
  /*
    Given the small dimensions of the domain, scaling it can help in
    achieving good precision. Using 1./k0 as a scale factor results in
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  data.scale = data.lambda/(2*M_PI);
  data.mode = wgma::planarwg::mode::TM;
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
  constexpr int nEigenpairs_top{300};
  constexpr int nEigenpairs_bottom{300};
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
    

    /********************************
   * cmesh(modal analysis):bottom   *
   ********************************/
    modal_bottom_an = ComputeModalAnalysis(gmesh, gmshmats,
                                           periodic_els,
                                           epsilon_copper,simdata,
                                           target_bottom, nEigenpairs_bottom,
                                           sortingRule, usingSLEPC,
                                           "bottom");
  
    /********************************
     * cmesh(modal analysis):top   *
     ********************************/
    modal_top_an = ComputeModalAnalysis(gmesh, gmshmats,
                                        periodic_els,
                                        epsilon_air,simdata,
                                        target_top, nEigenpairs_top,
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
  wgma::cmeshtools::RemovePeriodicity(modal_bottom_an->GetMesh());
  return 0;
}

#include <TPZKrylovEigenSolver.h>
#include <TPZPardisoSolver.h>

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
void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WgbcData& src_data,
                          WgbcData& match_data,
                          const int nmodes,
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
    
  auto scatt_mesh = [gmesh,&gmshmats,&simdata, &periodic_els]
    (CSTATE e_rib, CSTATE e_copper, CSTATE e_air){
    const auto mode = simdata.mode;
    const auto &pOrder = simdata.porder;
    const auto &lambda = simdata.lambda;
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
      probeMats.push_back("bnd_bottom");
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
                                                       src_ids,lambda,scale);
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
  TPZVec<int> nmodes = {1,10,15,100,200,250};
  //now we solve varying the number of modes used in the wgbc
  src_an->LoadAllSolutions();
  match_an->LoadAllSolutions();
  
  const auto mode = simdata.mode;
  //compute wgbc coefficients
  WgbcData src_data;
  src_data.cmesh = src_an->GetMesh();
  ComputeWgbcCoeffs(*src_an, src_data.wgbc_k, src_data.wgbc_f, mode, false, src_coeffs);

  WgbcData match_data;
  match_data.cmesh = match_an->GetMesh();
  ComputeWgbcCoeffs(*match_an, match_data.wgbc_k, match_data.wgbc_f,mode, true, {});

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
  const std::string suffix = simdata.rib_copper ? "_copper" : "_az";
  const std::string scatt_file = simdata.prefix+"_scatt"+suffix;
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_imag",
    "Field_abs"};
  auto vtk = TPZVTKGenerator(scatt_mesh, fvars, scatt_file, simdata.vtkRes);
  vtk.SetNThreads(simdata.nThreads);
    
  int count{0};
  for(int im = 0; im < nmodes.size(); im++){
    const int nm = nmodes[im];
    if(!nm){continue;}
    if(im != 0){
      wgma::cmeshtools::SetPeriodic(scatt_mesh, periodic_els);
    }
    RestrictDofsAndSolve(scatt_mesh, src_data, match_data, nm, simdata);
    //plot
    vtk.Do();
    //removing restrictions
    wgma::cmeshtools::RemovePeriodicity(scatt_mesh);
  }
}

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WgbcData& src_data,
                          WgbcData& match_data,
                          const int nmodes,
                          const SimData &simdata)
{

  auto gmesh = scatt_mesh->Reference();
  

  auto match_mesh = match_data.cmesh;
  auto src_mesh = src_data.cmesh;

  /**
     no dirichlet connects must be restricted!!
   **/
  std::set<int64_t> boundConnects;
  wgma::cmeshtools::FindDirichletConnects(scatt_mesh, boundConnects);



  const int64_t indep_con_id_match =
    wgma::cmeshtools::RestrictDofs(scatt_mesh, match_mesh, nmodes, boundConnects);
  const int64_t indep_con_id_src =
    wgma::cmeshtools::RestrictDofs(scatt_mesh, src_mesh, nmodes, boundConnects);

  constexpr bool sym{false};
  auto scatt_an = wgma::scattering::Analysis(scatt_mesh, simdata.nThreads,
                                             simdata.optimizeBandwidth,
                                             simdata.filterBoundEqs,
                                             sym);

  
  std::cout<<"nmodes on outgoing boundary: "<<nmodes<<std::endl;

  {
    TPZSimpleTimer timer("Assemble_"+std::to_string(nmodes),true);
    std::cout<<"Assembling...";
    scatt_an.Assemble();
    std::cout<<"\nAssembled!"<<std::endl;
  }
  {
    TPZSimpleTimer timer("Add contribution_"+std::to_string(nmodes),true);
  //now we must add the waveguide port terms
    AddWaveguidePortContribution(scatt_an, indep_con_id_match,
                                 nmodes, match_data.wgbc_k, match_data.wgbc_f);
    AddWaveguidePortContribution(scatt_an, indep_con_id_src,
                                 nmodes, src_data.wgbc_k, src_data.wgbc_f);
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
