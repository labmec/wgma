/**
slab_disc.cpp

This target performs the modal analysis of a dielectric slab
and then the subsequent scattering analysis at a waveguide discontinuity.
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
  STATE lambda{1.5};
  //!geometric scaling (floating point precision)
  REAL scale{1};
  //!reffractive index of core
  STATE ncore{1};
  //!reffractive index of clad
  STATE nclad{1};
  //!mode currently analysed
  wgma::planarwg::mode mode{wgma::planarwg::mode::TE};
  //!polynomial order
  int porder{-1};
  //!pml attenuation constant in x-direction
  CSTATE alphaPMLx{0};
  //!pml attenuation constant in y-direction
  CSTATE alphaPMLy{0};
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


TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const CSTATE target, const int neigenpairs,
            TPZEigenSort sorting, bool usingSLEPC);

void ComputeModes(wgma::wganalysis::WgmaPlanar &an,
                  TPZFMatrix<CSTATE> &sol_conj,
                  const bool orthogonalise,
                  const int nThreads);

void ComputeCouplingMat(wgma::wganalysis::WgmaPlanar &an,
                        std::string filename);
void ComputeCouplingMatTwoMeshes(wgma::wganalysis::WgmaPlanar &an_orig,
                                 const TPZFMatrix<CSTATE> &sol_conj,
                                 std::string filename);

void PostProcessModes(wgma::wganalysis::WgmaPlanar &an,
                      std::string filename,
                      const int vtkres);

void SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
                     wgma::wganalysis::WgmaPlanar &src_an,
                     TPZFMatrix<CSTATE> &src_adjoint_sol,
                     wgma::wganalysis::WgmaPlanar &match_an,
                     TPZFMatrix<CSTATE> &match_adjoint_sol,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const SimData &simdata);

void SolveWithPML(TPZAutoPointer<TPZCompMesh> scatt_cmesh,
                  wgma::wganalysis::WgmaPlanar& src_an,
                  const TPZVec<CSTATE> &src_coeffs,
                  const SimData &simdata);

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WgbcData& src_data,
                          WgbcData& match_data,
                          const TPZVec<CSTATE> &source_coeffs,
                          const int nmodes,
                          const SimData &simdata);

void ComputeWgbcCoeffs(wgma::wganalysis::WgmaPlanar& an,
                       const TPZFMatrix<CSTATE> &sol_conj,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       const bool positive_z, const TPZVec<CSTATE> &coeff);


SimData GetSimData()
{

  // path for output files
  const std::string path {"res_slab_disc/"};
  // common prefix for both meshes and output files
  const std::string basisName{"slab_disc"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  
  SimData data;
  data.meshfile="meshes/slab_disc.msh";
  data.lambda = 1.55;
  /*
    Given the small dimensions of the domain, scaling it can help in
    achieving good precision. Using 1./k0 as a scale factor results in
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  data.scale = data.lambda/(2*M_PI);
  data.mode = wgma::planarwg::mode::TE;
  data.ncore = 1.55;
  data.nclad = 1.00;
  data.alphaPMLx = {0.8, 0.0};
  data.alphaPMLy = {0.8, 0.0};
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
  constexpr int nEigenpairs_left{500};
  constexpr int nEigenpairs_right{500};
  const CSTATE target{simdata.ncore*simdata.ncore};

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
  auto gmesh = wgma::gmeshtools::ReadGmshMesh(simdata.meshfile, simdata.scale,
                                              gmshmats,
                                              verbosity_lvl);

  // print wgma_gmesh to .txt and .vtk format
  if (simdata.printGmesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = simdata.prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }


  /********************************
   * cmesh(modal analysis: left)  *
   ********************************/

  auto modal_l_cmesh = [gmesh,&gmshmats, &simdata](){
    // setting up cmesh data
    const auto &nclad = simdata.nclad;
    const auto &ncore = simdata.ncore;
    const auto &alphaPMLx = simdata.alphaPMLx;
    const auto &alphaPMLy = simdata.alphaPMLy;
    const auto &mode = simdata.mode;
    const auto &pOrder = simdata.porder;
    const auto &lambda = simdata.lambda;
    const auto &scale = simdata.scale;
    
    wgma::cmeshtools::PhysicalData modal_data;
    std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
    modal_mats["source_clad_left"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad, 1.);
    modal_mats["source_core_left"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    std::map<std::string, wgma::bc::type> modal_bcs;
    modal_bcs["source_left_bnd"] = wgma::bc::type::PEC;
    //dimension of the modal analysis 
    constexpr int modal_dim{1};
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                            {alphaPMLx,alphaPMLy}, modal_data, modal_dim);
    //we must now filter the 1D PMLs
    std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
    for(const auto &pml : modal_data.pmlvec){
      const std::string pattern{"source_clad_left"};
      const auto rx = std::regex{pattern, std::regex_constants::icase };
    
      const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
      if(found_pattern){pmlvec.push_back(pml);}
    }
    modal_data.pmlvec = pmlvec;
    return wgma::wganalysis::CMeshWgma1D(gmesh,mode,pOrder,modal_data,
                                         lambda, scale);
  }();
  /******************************
   * solve(modal analysis left) *
   ******************************/
  auto solver_left = SetupSolver(target, nEigenpairs_left, sortingRule, usingSLEPC);

  wgma::wganalysis::WgmaPlanar
    modal_l_an(modal_l_cmesh, simdata.nThreads,
             simdata.optimizeBandwidth, simdata.filterBoundEqs);
  modal_l_an.SetSolver(*solver_left);

  std::string modal_left_file{simdata.prefix+"_modal_left"};
  //no need to orthogonalise modes
  constexpr bool ortho{true};
  TPZFMatrix<CSTATE> sol_l_conj;
  modal_l_an.Assemble();
  constexpr bool export_mats{false};
  //now we export the matrices
  if(export_mats){
    auto mata = modal_l_an.GetSolver().MatrixA();
    auto matb = modal_l_an.GetSolver().MatrixB();
    std::ofstream mata_file(simdata.prefix+"_mata.csv");
    mata->Print("",mata_file,ECSV);
    std::ofstream matb_file(simdata.prefix+"_matb.csv");
    matb->Print("",matb_file,ECSV);
  }
  ComputeModes(modal_l_an, sol_l_conj, ortho, simdata.nThreads);
  if(export_mats) {
    auto &eqfilt = modal_l_an.StructMatrix()->EquationFilter();
    auto &ev = modal_l_an.GetEigenvectors();
    const int nrows = eqfilt.NActiveEquations();
    const int ncols = ev.Cols();
    TPZFMatrix<CSTATE> ev_gather(nrows,ncols,0.);
    eqfilt.Gather(ev, ev_gather);
    // auto matp = &ev_gather.g(0,0);
    // for(int i = 0; i < nrows*ncols; i++){
    //   *matp++ = std::conj(*matp);
    // }
    std::ofstream sol_file(simdata.prefix+"_ev.csv");
    ev_gather.Print("",sol_file,ECSV);
  }
  
  if(simdata.couplingmat){
    ComputeCouplingMat(modal_l_an, simdata.prefix+"_mat_src.csv");
    ComputeCouplingMatTwoMeshes(modal_l_an, sol_l_conj,
                                simdata.prefix+"_mat_src_conj.csv");
  }
  if(simdata.exportVtk){
    PostProcessModes(modal_l_an, modal_left_file, simdata.vtkRes);
  }
  

  /********************************
   * cmesh(modal analysis: right)  *
   ********************************/

  auto modal_r_cmesh = [gmesh,&gmshmats,&simdata](){
    const auto &nclad = simdata.nclad;
    const auto &ncore = simdata.ncore;
    const auto &alphaPMLx = simdata.alphaPMLx;
    const auto &alphaPMLy = simdata.alphaPMLy;
    const auto &mode = simdata.mode;
    const auto &pOrder = simdata.porder;
    const auto &lambda = simdata.lambda;
    const auto &scale = simdata.scale;
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData modal_data;

    std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
    modal_mats["source_clad_right"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad, 1.);
    modal_mats["source_core_right"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    std::map<std::string, wgma::bc::type> modal_bcs;
    modal_bcs["source_right_bnd"] = wgma::bc::type::PEC;
    //dimension of the modal analysis 
    constexpr int modal_dim{1};
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                            {alphaPMLx,alphaPMLy}, modal_data, modal_dim);

    //we must now filter the 1D PMLs
    std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
    for(const auto &pml : modal_data.pmlvec){
      const std::string pattern{"source_clad_right"};
      const auto rx = std::regex{pattern, std::regex_constants::icase };
    
      const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
      if(found_pattern){pmlvec.push_back(pml);}
    }
    modal_data.pmlvec = pmlvec;
    return wgma::wganalysis::CMeshWgma1D(gmesh,mode,pOrder,modal_data,
                                         lambda, scale);
  }();

  /******************************
   * solve(modal analysis: right) *
   ******************************/
  auto solver_right = SetupSolver(target, nEigenpairs_right, sortingRule, usingSLEPC);
  
  wgma::wganalysis::WgmaPlanar
    modal_r_an(modal_r_cmesh, simdata.nThreads,
             simdata.optimizeBandwidth, simdata.filterBoundEqs);
  modal_r_an.SetSolver(*solver_right);

  std::string modal_right_file{simdata.prefix+"_modal_right"};
  TPZFMatrix<CSTATE> sol_r_conj;
  modal_r_an.Assemble();
  ComputeModes(modal_r_an, sol_r_conj, ortho, simdata.nThreads);
  if(simdata.couplingmat){
    ComputeCouplingMat(modal_r_an, simdata.prefix+"_mat_match.txt");
  }
  if(simdata.exportVtk){
    std::string modal_left_file{simdata.prefix+"_modal_right"};
    PostProcessModes(modal_r_an, modal_right_file, simdata.vtkRes);
  }

  SolveScattering(gmesh, modal_l_an, sol_l_conj,
                  modal_r_an, sol_r_conj,
                  gmshmats,simdata);
  return 0;
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
  const int krylovDim = std::max(5,3*neigenpairs);
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
    constexpr STATE ksp_rtol = 1e-12;//PETSC_DECIDE
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

#include <TPZYSMPMatrix.h>
#include <post/orthowgsol.hpp>
void ComputeModes(wgma::wganalysis::WgmaPlanar &an,
                  TPZFMatrix<CSTATE> &sol_conj,
                  const bool orthogonalise,
                  const int nThreads)
{
  
  TPZSimpleTimer analysis("Modal analysis");
  static constexpr bool computeVectors{true};
  

  an.Solve(computeVectors);
  constexpr bool really_compute_adj{false};
  //currently we cannot orthogonalise solutions from adj problem
  if(really_compute_adj){
    /*we also need the solutions of the adjoint problem*/

    //our solver can either be a EPSHandler or KrylovEigenSolver
    TPZAutoPointer<TPZLinearEigenSolver<CSTATE>> solver_conj = [&an]()
      -> TPZLinearEigenSolver<CSTATE>*
      {
        auto clone = an.GetSolver().Clone();
        auto krylov = dynamic_cast<TPZKrylovEigenSolver<CSTATE>*>(clone);
        if(krylov){return krylov;}
        else{
          auto eps = dynamic_cast<wgma::slepc::EPSHandler<CSTATE>*>(clone);
          if(eps){return eps;}
        }
        DebugStop();
        return nullptr;
      }();

    /*
      modify matrices to represent adjoint problem
    */
    auto ConjMatrix = [](TPZAutoPointer<TPZMatrix<CSTATE>> mat){
      auto sparse_mat =
        TPZAutoPointerDynamicCast<TPZFYsmpMatrix<CSTATE>>(mat);
      if(!sparse_mat){
        DebugStop();
      }
      int64_t *ia{nullptr}, *ja{nullptr};
      CSTATE *avec{nullptr};
      sparse_mat->GetData(ia,ja,avec);
      const int nrows = sparse_mat->Rows();
      //number of non-zero entries in the sparse matrix
      const int last_pos=ia[nrows];
      for(int i = 0; i < last_pos; i++){
        *avec++ = std::conj(*avec);
      }
    };
    //copying the matrices
    solver_conj->SetMatrixA(solver_conj->MatrixA()->Clone());
    solver_conj->SetMatrixB(solver_conj->MatrixB()->Clone());

    ConjMatrix(solver_conj->MatrixA());
    ConjMatrix(solver_conj->MatrixB());
    /*
      now we need to know:
      n reduced equations (eq filter)
      n equations
      n eigenvalues
    */
    const auto nev = solver_conj->NEigenpairs();
    auto &eqfilt = an.StructMatrix()->EquationFilter();
    const auto neq = eqfilt.NEqExpand();
    const auto neqreduced = eqfilt.NActiveEquations();

    sol_conj.Redim(neq,nev);
    TPZVec<CSTATE> w_conj(nev, 0.);
    TPZFMatrix<CSTATE> sol_conj_gather(neqreduced,nev);
    solver_conj->Solve(w_conj, sol_conj_gather);
    //this is the solution that should be loaded in the computational mesh
    eqfilt.Scatter(sol_conj_gather,sol_conj);
  }
  
  //load all obtained modes into the mesh
  an.LoadAllSolutions();

  if(orthogonalise){
    TPZSimpleTimer timer("Ortho",true);
    constexpr STATE tol{1e-8};
    const int n_ortho = wgma::post::OrthoWgSol(an,tol);
    std::cout<<"orthogonalised  "<<n_ortho<<" eigenvectors"<<std::endl;
  }
  else{
    //we just normalise
    auto cmesh = an.GetMesh();
    constexpr bool conj{false};
    //leave empty for all valid matids
    std::set<int> matids {};
    wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(cmesh,matids,conj,nThreads).Normalise();
  
    TPZFMatrix<CSTATE> &mesh_sol=cmesh->Solution();
    //we update analysis object
    an.SetEigenvectors(mesh_sol);
  }

  if(!really_compute_adj){
    sol_conj = an.GetEigenvectors();
    const int nr = sol_conj.Rows();
    const int nc = sol_conj.Cols();
    CSTATE *ptr = sol_conj.Elem();
    const auto sz = nr*nc;
    for(int i = 0; i < sz; i++, ptr++){
      *ptr = std::conj(*ptr);
    }
  }
    
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
  integrator.ComputeCoupling();
  TPZFMatrix<CSTATE> couplingmat;
  integrator.GetCoupling(couplingmat);
  std::ofstream matfile(filename);
  couplingmat.Print("",matfile,ECSV);
}

#include <Electromagnetics/TPZScalarField.h>
#include <TPZNullMaterialCS.h>
#include <pzbuildmultiphysicsmesh.h>

void ComputeCouplingMatTwoMeshes(wgma::wganalysis::WgmaPlanar &an,
                                 const TPZFMatrix<CSTATE> &sol_conj,
                                 std::string filename)
{
  
  using namespace wgma::post;
  auto mesh = an.GetMesh();
  auto gmesh = mesh->Reference();
  const auto dim = mesh->Dimension();
  //now we transfer the solutions from adjoint problem to original problem
  TPZFMatrix<CSTATE>& sol_orig = mesh->Solution();
  const int neq = sol_orig.Rows();
  const int nsol = sol_orig.Cols();
  if(nsol != sol_conj.Cols() || sol_conj.Rows() != sol_orig.Rows()) {
    //how to deal with dependencies in original mesh sol?
    DebugStop();
  }
  sol_orig.Resize(neq,2*nsol);

  CSTATE *ptr_orig = &sol_orig.g(0, nsol);
  CSTATE *ptr_conj = &sol_conj.g(0, 0);
  for(int ipos = 0; ipos < nsol*neq; ipos++){
    *ptr_orig++ = *ptr_conj++;
  }

  std::set<int> matids;
  constexpr bool conj{true};
  const int nthreads = std::thread::hardware_concurrency();
  
  WaveguideCoupling<SingleSpaceIntegrator> integrator(mesh,
                                                      matids,
                                                      conj,
                                                      nthreads
                                                      );

  integrator.SetAdjoint(true);
  integrator.ComputeCoupling();
  TPZFMatrix<CSTATE> couplingmat;
  integrator.GetCoupling(couplingmat);
  std::ofstream matfile(filename);
  couplingmat.Print("",matfile,ECSV);
  //now we resize it again to its original size
  sol_orig.Resize(neq,nsol);
}

void PostProcessModes(wgma::wganalysis::WgmaPlanar &an,
                      std::string filename,
                      const int vtkres){
  TPZSimpleTimer postProc("Post processing");
    
  const std::string file = filename;
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

  auto cmesh = an.GetMesh();
  auto vtk = TPZVTKGenerator(cmesh, fvars, file, vtkres);
  const int64_t maxval{100};
  const auto nsol = std::min(maxval,an.GetEigenvectors().Cols());
  std::cout<<"Exporting "<<nsol<<" solutions"<<std::endl;
  for(auto isol = 0; isol < nsol ; isol++){
    an.LoadSolution(isol);
    vtk.Do();
  }
}


#include <materials/solutionprojection.hpp>
#include <TPZNullMaterial.h>

void ProjectSolIntoRestrictedMesh(wgma::wganalysis::WgmaPlanar &src_an,
                                  TPZAutoPointer<TPZCompMesh> scatt_mesh,
                                  const SimData &simdata,
                                  const TPZVec<int> &nmodes)
{
  std::map<int,STATE> error_proj;
  //load all solutions
  src_an.LoadAllSolutions();
  //1d mesh
  TPZAutoPointer<TPZCompMesh> src_mesh = src_an.GetMesh();
  //mesh used to project the solution into the restricted space
  TPZAutoPointer<TPZCompMesh> proj_mesh = src_mesh->Clone();
  //replace all materials for L2 projection
  TPZMaterialT<CSTATE> *last_mat{nullptr};
  for(auto [id,mat] : proj_mesh->MaterialVec()){
    auto bnd = dynamic_cast<TPZBndCond*>(mat);
    if(!bnd){
      constexpr int dim{1};
      constexpr int soldim{1};
      auto newmat = new wgma::materials::SolutionProjection<CSTATE>(id,dim,soldim);
      delete mat;
      proj_mesh->MaterialVec()[id]=newmat;
      last_mat = newmat;
    }
  }
  for(auto [id,mat] : proj_mesh->MaterialVec()){
    auto bnd = dynamic_cast<TPZBndCond*>(mat);
    if(bnd){
      TPZFMatrix<CSTATE> v1;
      TPZVec<CSTATE> v2;
      auto newmat = last_mat->CreateBC(last_mat,id,0,v1,v2);
      delete bnd;
      proj_mesh->MaterialVec()[id]=dynamic_cast<TPZMaterial*>(newmat);
    }
  }
  TPZAutoPointer<TPZCompMesh> error_mesh = proj_mesh->Clone();

  //setup vtk objects
  const std::string proj_file = simdata.prefix+"_proj";
  auto vtk =
    TPZVTKGenerator(proj_mesh, {"Solution"}, proj_file, simdata.vtkRes);
  const std::string error_file = simdata.prefix+"_proj_error";
  auto vtk_error =
    TPZVTKGenerator(error_mesh, {"Solution"}, error_file, simdata.vtkRes);

  //now we get reference sol
  TPZFMatrix<CSTATE> sol_pml = proj_mesh->Solution();
  sol_pml.Redim(sol_pml.Rows(),1);
  std::cout<<"neqs without restriction: "<<sol_pml.Rows()<<std::endl;
  wgma::cmeshtools::ExtractSolFromMesh(proj_mesh, scatt_mesh, sol_pml);
  //just to set size
  auto sol_proj = sol_pml;



  {
    auto eqfilt = src_an.StructMatrix()->EquationFilter();
    auto neqcondense = eqfilt.NActiveEquations();
    std::string sol2dfilename = simdata.prefix+"_sol.csv";
    std::ofstream sol2dfile(sol2dfilename);
    TPZFMatrix<CSTATE> sol2d(neqcondense,1);
    eqfilt.Gather(sol_proj, sol2d);
    std::cout<<"sol2d dimensions: "<<sol2d.Rows()<<","<<sol2d.Cols()<<std::endl;
    sol2d.Print("",sol2dfile,ECSV);

    // std::string vmatfilename = simdata.prefix+"_ev.csv";
    // std::ofstream vmatfile(vmatfilename);
    // std::string vmatname = "vmat=";
    // auto &eigenvectors = src_an.GetEigenvectors();
    // TPZFMatrix<CSTATE> vmat(neqcondense,eigenvectors.Cols());
    // eqfilt.Gather(eigenvectors, vmat);
    // std::cout<<"vmat dimensions: "<<vmat.Rows()<<","<<vmat.Cols()<<std::endl;
    // vmat.Print("",vmatfile,ECSV);
  }
  
  /*
    dirichlet boundary connects should not be restricted, otherwise
    this will result in all the equations on the same dependency
    being removed as well
   */
  std::set<int64_t> bound_connects;
  wgma::cmeshtools::FindDirichletConnects(proj_mesh, bound_connects);

  auto error_an = wgma::scattering::Analysis(error_mesh, simdata.nThreads,
                                             false,
                                             simdata.filterBoundEqs,false);

  std::ofstream s_cmesh_file{simdata.prefix+"_scatt_mesh_.txt"};
  scatt_mesh->Print(s_cmesh_file);

  constexpr bool export_mats{false};
  TPZFMatrix<CSTATE> couplingmat;
  {
    using namespace wgma::post;
    std::set<int> matids;
    constexpr bool conj{true};
    const int nthreads = std::thread::hardware_concurrency();
    WaveguideCoupling<SingleSpaceIntegrator> integrator(src_an.GetMesh(),
                                                        matids,
                                                        conj,
                                                        nthreads
                                                        );
    TPZVec<CSTATE> betavec = src_an.GetEigenvalues();
    for(auto &b : betavec){b = sqrt(b);}
    integrator.SetBeta(betavec);
    //we want just std::conj(et_i) * et_j
    integrator.SetMu(false);
    integrator.ComputeCoupling();
    integrator.GetCoupling(couplingmat);
  }

  if(export_mats){
    std::ofstream couplfile{simdata.prefix+"_proj_coupl.csv"};
    couplingmat.Print("",couplfile,ECSV);
    couplfile.close();
  }
  for(int im = 0; im < nmodes.size(); im++){
    const int nm = nmodes[im];
    //restrict dofs in proj mesh
    if(nm){
      wgma::cmeshtools::RestrictDofs(proj_mesh, src_mesh, nm, bound_connects);
      TPZFMatrix<CSTATE> &sol = proj_mesh->Solution();
      std::cout<<"neqs after restriction: "<<sol.Rows()<<std::endl;
    }
    /*The TPZLinearAnalysis class manages the creation of the algebric
     * problem and the matrix inversion*/
    auto proj_an = wgma::scattering::Analysis(proj_mesh, simdata.nThreads,
                                              false,
                                              simdata.filterBoundEqs,false);
    std::ofstream cmesh_file{simdata.prefix+"_proj_mesh_"+std::to_string(nm)+".txt"};
    proj_mesh->Print(cmesh_file);
    //now we load desired solution into proj_mesh so we can project it
    {
      TPZFMatrix<CSTATE> &sol_reference = proj_mesh->Solution();
      wgma::cmeshtools::ExtractSolFromMesh(proj_mesh, scatt_mesh, sol_reference);
    }
        
    proj_an.Assemble();    

    if(export_mats){
      std::ofstream projfile{simdata.prefix+"_proj_mat_"+std::to_string(nm)+".csv"};
      auto mat = proj_an.GetSolver().Matrix();
      mat->Print("",projfile,ECSV);
      projfile.close();
    }

    // if(nm){
    //   //this works if coupling mat is computed with ur=1
    //   auto mat = proj_an.GetSolver().Matrix();
    //   for(int i = 0; i < nm; i++){
    //     for(int j = 0; j < nm; j++){
    //       mat->Put(i,j,couplingmat.Get(i,j));
    //     }
    //   }
    // }
    
    proj_an.Solve();

    
    //plot
    vtk.Do();
    //now we compute the error
    wgma::cmeshtools::ExtractSolFromMesh(error_mesh, proj_mesh, sol_proj);
    
    sol_proj -= sol_pml;

    error_mesh->LoadSolution(sol_proj);
    vtk_error.Do();
    auto normsol =
      wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
    normsol.SetNThreads(std::thread::hardware_concurrency());
    const auto norm = std::real(normsol.ComputeNorm()[0]);
    std::cout<<"nmodes "<<nm<<" error "<<norm<<std::endl;
    error_proj.insert({nm,norm});
    if(nm){
      wgma::cmeshtools::RemovePeriodicity(proj_mesh);
      proj_mesh->InitializeBlock();
    }
  }

  std::cout<<"**********PROJECTION**********"<<std::endl;
  std::cout<<"nmodes  norm error: "<<std::endl;
  for(auto [nm,error] : error_proj){
    std::cout<<nm<<","<<error<<std::endl;
  }
  std::cout<<"******************************"<<std::endl;
}

void SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
                     wgma::wganalysis::WgmaPlanar &src_an,
                     TPZFMatrix<CSTATE> &src_adjoint_sol,
                     wgma::wganalysis::WgmaPlanar &match_an,
                     TPZFMatrix<CSTATE> &match_adjoint_sol,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const SimData &simdata)
{
  
  auto CreateScattMesh = [gmesh,&gmshmats,&simdata]
    (bool extendDomains){

    const CSTATE &ncore = simdata.ncore;
    const CSTATE &nclad = simdata.nclad;
  
    const CSTATE alphaPMLx = simdata.alphaPMLx;
    const CSTATE alphaPMLy = simdata.alphaPMLy;

    const auto &mode = simdata.mode;
    const auto &pOrder = simdata.porder;
    const auto &lambda = simdata.lambda;
    const auto &scale = simdata.scale;
    
    const std::string &prefix = simdata.prefix;

    
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData scatt_data;
    std::map<std::string, std::pair<CSTATE, CSTATE>> scatt_mats;
    scatt_mats["core_left"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    scatt_mats["core_right"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    scatt_mats["cladding_left"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad,1.);
    scatt_mats["cladding_right"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad,1.);
    std::map<std::string, wgma::bc::type> scatt_bcs;
    scatt_bcs["scatt_bnd_mid"] = wgma::bc::type::PEC;
    if(extendDomains){
      scatt_bcs["scatt_bnd_right"] = wgma::bc::type::PEC;
      scatt_bcs["scatt_bnd_left"] = wgma::bc::type::PEC;
    }

    
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, scatt_mats, scatt_bcs,
                                            {alphaPMLx,alphaPMLy}, scatt_data);


    //materials that will represent our source
    std::set<int> src_ids;
    /*
      if the domain is extended, we will create a current source at
      source_core_left and source_clad_left.
      otherwise, the source is modelled as a waveguide port boundary
      condition and there is no need for source mats
     */
    if (extendDomains){
      const std::string srcMat[] = {"source_core_left", "source_clad_left"};
      for(const auto &mat : srcMat){
        src_ids.insert(gmshmats[1].at(mat));
      }
      //now we add the 1d pml mats to source mats
      for(const auto &pml : scatt_data.pmlvec){
        const std::string pattern{"source_clad_left"};
        const auto rx = std::regex{pattern, std::regex_constants::icase };
    
        const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
        if(found_pattern){
          const auto matdim = 1;
          const auto id = gmshmats[matdim].at(*pml->names.begin());
          src_ids.insert(id);
        }
      }
    }
    


    /*
      probe mats are regions of the domain in which we want to be able
      to evaluate our solution
      they are also used to ensure that the computational elements are created
      so every region that will be used in a waveguide port bc must be
      also inserted as a probe mat
    */
    std::vector<std::string> probeMats;
    {
      probeMats.push_back("source_clad_right");
      probeMats.push_back("source_core_right");
      probeMats.push_back("source_clad_left");
      probeMats.push_back("source_core_left");
      probeMats.push_back("eval_clad_left");
      probeMats.push_back("eval_core_left");
      //now for 1d pml mats
      for(const auto &pml : scatt_data.pmlvec){
        const std::string pattern_left{"source_clad_left"};
        const auto rx_left =
          std::regex{pattern_left, std::regex_constants::icase };
        const std::string pattern_right{"source_clad_right"};
        const auto rx_right =
          std::regex{pattern_right, std::regex_constants::icase };
        const std::string pattern_eval{"eval_clad_left"};
        const auto rx_eval =
          std::regex{pattern_eval, std::regex_constants::icase };
        
        const bool found_pattern =
          std::regex_search(*(pml->names.begin()), rx_left) ||
          std::regex_search(*(pml->names.begin()), rx_right) ||
          std::regex_search(*(pml->names.begin()), rx_eval);

        if(found_pattern){
          probeMats.push_back(*pml->names.begin());
        }
      }
      for(auto &mat : probeMats){
        const auto matdim = 1;
        const auto id = gmshmats[matdim].at(mat);
        if(src_ids.count(id)==0){//just to avoid inserting materials twice when extending domains
          scatt_data.probevec.push_back({id,matdim});
        }
      }
    }
    
    

    
    
    /*when using waveguide port bc we need to filter out some PML regions*/
    if (!extendDomains){
      std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
      for(const auto &pml : scatt_data.pmlvec){
        const std::string pattern_left{"xm"};
        const auto rx_left = std::regex{pattern_left, std::regex_constants::icase };
        const std::string pattern_right{"xp"};
        const auto rx_right = std::regex{pattern_right, std::regex_constants::icase };
    
        const bool found_pattern =
          std::regex_search(*(pml->names.begin()), rx_left) ||
          std::regex_search(*(pml->names.begin()), rx_right);
        
        if(found_pattern){//we discard these regions
          continue;
        }
        pmlvec.push_back(pml);
      }
      scatt_data.pmlvec = pmlvec;
    }

    return wgma::scattering::CMeshScattering2D(gmesh, mode, pOrder, scatt_data,src_ids,
                                               lambda,scale);
  };
  /*********************
   * solve(scattering) *  
   *********************/  
  TPZSimpleTimer tscatt("Scattering");


  const auto n_eigenpairs_left = src_an.GetEigenvalues().size();
  /*
    the source is written as a linear combination of the modes
    this vector contains the coefficients of such combination
   */
  TPZVec<CSTATE> src_coeffs(n_eigenpairs_left,0);
  src_coeffs[0] = 1;

  //set up post processing
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_imag",
    "Field_abs"};

  
  auto scatt_mesh_pml = CreateScattMesh(true);
  //solve using PML as a reference solution
  {
    const std::string suffix = "pml";
    const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
    auto vtk = TPZVTKGenerator(scatt_mesh_pml, fvars, scatt_file, simdata.vtkRes+2);
    SolveWithPML(scatt_mesh_pml,src_an,src_coeffs,simdata);
    vtk.Do();
  }
  //now we solve varying the number of modes used in the wgbc
  
  //index of the number of modes to be used to restrict the dofs on waveguide bcs
  TPZVec<int> nmodes = {0,1,2,5,10,15,20,100,300};
  src_an.LoadAllSolutions();
  match_an.LoadAllSolutions();

  //as an initial test, one could just simply project the solution and check
  //the results
  constexpr bool test_proj{false};
  if(test_proj){
    ProjectSolIntoRestrictedMesh(src_an, scatt_mesh_pml, simdata, nmodes);
  }
  
  auto scatt_mesh_wgbc = CreateScattMesh(false);
  const std::string suffix = "wgbc";
  const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
  auto vtk = TPZVTKGenerator(scatt_mesh_wgbc, fvars, scatt_file, simdata.vtkRes);

  //compute wgbc coefficients
  WgbcData src_data;
  src_data.cmesh = src_an.GetMesh();
  ComputeWgbcCoeffs(src_an, src_adjoint_sol,
                    src_data.wgbc_k, src_data.wgbc_f, false, src_coeffs);

  WgbcData match_data;
  match_data.cmesh = match_an.GetMesh();
  ComputeWgbcCoeffs(match_an, match_adjoint_sol,
                    match_data.wgbc_k, match_data.wgbc_f,true, {});
    
    
  //here we will store the error between pml approx and wgbc approx
  /*
  TPZAutoPointer<TPZCompMesh> error_mesh = [gmesh,&gmshmats, &simdata](){
    // setting up cmesh data
    const auto &nclad = simdata.nclad;
    const auto &ncore = simdata.ncore;
    const auto &alphaPMLx = simdata.alphaPMLx;
    const auto &alphaPMLy = simdata.alphaPMLy;
    const auto &mode = simdata.mode;
    const auto &pOrder = simdata.porder;
    const auto &lambda = simdata.lambda;
    const auto &scale = simdata.scale;
    
    wgma::cmeshtools::PhysicalData modal_data;
    std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
    modal_mats["eval_clad_left"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad, 1.);
    modal_mats["eval_core_left"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    std::map<std::string, wgma::bc::type> modal_bcs;
    modal_bcs["eval_left_bnd"] = wgma::bc::type::PEC;
    //dimension of the modal analysis 
    constexpr int modal_dim{1};
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                            {alphaPMLx,alphaPMLy}, modal_data, modal_dim);
    //we must now filter the 1D PMLs
    std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
    for(const auto &pml : modal_data.pmlvec){
      const std::string pattern{"eval_clad_left"};
      const auto rx = std::regex{pattern, std::regex_constants::icase };
    
      const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
      if(found_pattern){pmlvec.push_back(pml);}
    }
    modal_data.pmlvec = pmlvec;
    return wgma::wganalysis::CMeshWgma1D(gmesh,mode,pOrder,modal_data,
                                         lambda, scale);
  }();
  */
  TPZAutoPointer<TPZCompMesh> error_mesh = src_an.GetMesh()->Clone();
  
  TPZFMatrix<CSTATE> sol_pml = error_mesh->Solution();
  {
    wgma::cmeshtools::ExtractSolFromMesh(error_mesh, scatt_mesh_pml, sol_pml);
  }
  //just to get the same size, we will zero it later
  auto sol_wgbc = sol_pml;
  const std::string error_file = simdata.prefix+"_error_";
  auto vtk_error = TPZVTKGenerator(error_mesh, fvars, error_file, simdata.vtkRes);
  std::map<int,STATE> error_res;
  for(int im = 0; im < nmodes.size(); im++){
    const int nm = nmodes[im];
    if(!nm){continue;}
    RestrictDofsAndSolve(scatt_mesh_wgbc, src_data, match_data,
                         src_coeffs, nm,simdata);
    //plot
    vtk.Do();
    //now we compute the error
    wgma::cmeshtools::ExtractSolFromMesh(error_mesh, scatt_mesh_wgbc, sol_wgbc);
    sol_wgbc -= sol_pml;
    error_mesh->LoadSolution(sol_wgbc);
    vtk_error.Do();
    auto normsol =
      wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
    normsol.SetNThreads(std::thread::hardware_concurrency());
    const auto norm = std::real(normsol.ComputeNorm()[0]);
    std::cout<<"nmodes "<<nm<<" error "<<norm<<std::endl;
    error_res.insert({nm,norm});
      
    //removing restrictions
    wgma::cmeshtools::RemovePeriodicity(scatt_mesh_wgbc);
    scatt_mesh_wgbc->ComputeNodElCon();
    scatt_mesh_wgbc->CleanUpUnconnectedNodes();
  }
  std::cout<<"++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  for(auto [nm,error] : error_res){
    std::cout<<"nmodes "<<nm<<" norm error: "<<error<<std::endl;
  }
  std::cout<<"++++++++++++++++++++++++++++++++++++++++"<<std::endl;
}

#include <TPZPardisoSolver.h>
#include <TPZYSMPPardiso.h>

void SolveWithPML(TPZAutoPointer<TPZCompMesh> scatt_cmesh,
                  wgma::wganalysis::WgmaPlanar& src_an,
                  const TPZVec<CSTATE> &src_coeffs,
                  const SimData &simdata)
{
  constexpr bool sym{false};
  auto scatt_an = wgma::scattering::Analysis(scatt_cmesh, simdata.nThreads,
                                             simdata.optimizeBandwidth,
                                             simdata.filterBoundEqs,
                                             sym);
    
  auto src_mesh = src_an.GetMesh();
  //get id of source materials
  wgma::scattering::SourceWgma src;
  for(auto [id,mat] : src_mesh->MaterialVec()){
    src.id.insert(id);
  }
  src.modal_cmesh = src_mesh;

  TPZFMatrix<CSTATE> sol(scatt_cmesh->NEquations(),1);
  //we will add the components of the solution one by one
  const int nsol = src_coeffs.size();
  bool first_assemble{true};
  for(int isol = 0; isol < nsol; isol++){
    if(IsZero(src_coeffs[isol])){continue;}
    src_an.LoadSolution(isol);
    auto beta = std::sqrt(src_an.GetEigenvalues()[isol]);
    wgma::scattering::LoadSource1D(scatt_cmesh, src, src_coeffs[isol]);
    wgma::scattering::SetPropagationConstant(scatt_cmesh, beta);
    if(first_assemble){
      first_assemble = false;
      scatt_an.Assemble();
    }else{
      scatt_an.AssembleRhs(src.id);
    }
    scatt_an.Solve();
    scatt_an.LoadSolution();
    TPZFMatrix<CSTATE> &curr_sol = scatt_an.Solution();
    sol+=curr_sol;
  }
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


void ComputeWgbcCoeffs(wgma::wganalysis::WgmaPlanar& an,
                       const TPZFMatrix<CSTATE> &sol_conj,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       const bool positive_z, const TPZVec<CSTATE> &coeff){
  auto mesh = an.GetMesh();

  TPZFMatrix<CSTATE>& sol_orig = mesh->Solution();
  const int neq = sol_orig.Rows();
  const int nsol = sol_orig.Cols();

  constexpr bool adj{false};
  if(adj){
    
    if(nsol != sol_conj.Cols() || sol_conj.Rows() != sol_orig.Rows()) {
      //how to deal with dependencies in original mesh sol?
      DebugStop();
    }
    sol_orig.Resize(neq,2*nsol);

    CSTATE *ptr_orig = &sol_orig.g(0, nsol);
    CSTATE *ptr_conj = &sol_conj.g(0, 0);
    for(int ipos = 0; ipos < nsol*neq; ipos++){
      *ptr_orig++ = *ptr_conj++;
    }
  }
  
  wgma::post::WaveguidePortBC<wgma::post::SingleSpaceIntegrator> wgbc(mesh);
  TPZManVector<CSTATE,1000> betavec = an.GetEigenvalues();
  for(auto &b : betavec){b = std::sqrt(b);}
  wgbc.SetPositiveZ(positive_z);
  wgbc.SetAdjoint(adj);
  if(coeff.size()){
    wgbc.SetSrcCoeff(coeff);
  }
  wgbc.SetBeta(betavec);
  wgbc.ComputeContribution();
  wgbc.GetContribution(wgbc_k,wgbc_f);
  //now we resize it again to its original size
  if(adj){sol_orig.Resize(neq,nsol);}
}
void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WgbcData& src_data,
                          WgbcData& match_data,
                          const TPZVec<CSTATE> &source_coeffs,
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
    
  scatt_an.Assemble();
  //now we must add the waveguide port terms
  AddWaveguidePortContribution(scatt_an, indep_con_id_match,
                               nmodes, match_data.wgbc_k, match_data.wgbc_f);
  AddWaveguidePortContribution(scatt_an, indep_con_id_src,
                               nmodes, src_data.wgbc_k, src_data.wgbc_f);
  //get pardiso control
  auto *pardiso = scatt_an.GetSolver().GetPardisoControl();
  pardiso->SetMessageLevel(0);
  
  pardiso->ResetParam();
  constexpr auto sys_type = SymProp::Sym;
  constexpr auto prop = TPZPardisoSolver<CSTATE>::MProperty::EIndefinite;
  pardiso->SetMatrixType(sys_type,prop);
  TPZSimpleTimer tscatt("Solve",true);
  scatt_an.Solve();
}

