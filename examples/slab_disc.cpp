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
#include <post/orthosol.hpp>
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
                  const bool orthogonalise,
                  const int nThreads);

void ComputeCouplingMat(wgma::wganalysis::WgmaPlanar &an,
                        std::string filename,
                        std::string matname);

void PostProcessModes(wgma::wganalysis::WgmaPlanar &an,
                      std::string filename,
                      const int vtkres);

void SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
                     wgma::wganalysis::WgmaPlanar &src_an,
                     wgma::wganalysis::WgmaPlanar &match_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const SimData &simdata);

void SolveWithPML(TPZAutoPointer<TPZCompMesh> scatt_cmesh,
                  wgma::wganalysis::WgmaPlanar& src_an,
                  const TPZVec<CSTATE> &src_coeffs,
                  const SimData &simdata);

int64_t RestrictDofs(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                     TPZAutoPointer<TPZCompMesh> match_mesh,
                     const int nm,
                     const std::set<int64_t> bound_connects);

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WgbcData& src_data,
                          WgbcData& match_data,
                          const TPZVec<CSTATE> &source_coeffs,
                          const int nmodes,
                          const SimData &simdata);

void ComputeWgbcCoeffs(wgma::wganalysis::WgmaPlanar& an,
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
  data.alphaPMLx = 0.75;
  data.alphaPMLy = 0.75;
  data.porder = 2;
  data.filterBoundEqs = true;
  data.printGmesh=true;
  data.exportVtk = true;
  data.couplingmat = false;
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
  constexpr int nEigenpairs_left{50};
  constexpr int nEigenpairs_right{50};
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
  constexpr bool ortho{false};
  ComputeModes(modal_l_an, ortho, simdata.nThreads);
  if(simdata.couplingmat){
    ComputeCouplingMat(modal_l_an, simdata.prefix+"_mat_src.txt", "src");
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
  
  ComputeModes(modal_r_an, ortho, simdata.nThreads);
  if(simdata.couplingmat){
    ComputeCouplingMat(modal_r_an, simdata.prefix+"_mat_match.txt", "match");
  }
  if(simdata.exportVtk){
    std::string modal_left_file{simdata.prefix+"_modal_right"};
    PostProcessModes(modal_r_an, modal_right_file, simdata.vtkRes);
  }

  SolveScattering(gmesh, modal_l_an, modal_r_an,
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

void ComputeModes(wgma::wganalysis::WgmaPlanar &an,
                  const bool orthogonalise,
                  const int nThreads)
{
  
  TPZSimpleTimer analysis("Modal analysis");
  an.Assemble();
  static constexpr bool computeVectors{true};
  an.Solve(computeVectors);
  //load all obtained modes into the mesh
  an.LoadAllSolutions();

  auto cmesh = an.GetMesh();
  constexpr bool conj{false};
  //leave empty for all valid matids
  std::set<int> matids {};
  wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(cmesh,matids,conj,nThreads).Normalise();
  
  TPZFMatrix<CSTATE> &mesh_sol=cmesh->Solution();
  //we update analysis object
  an.SetEigenvectors(mesh_sol);
  if(orthogonalise){
    auto ortho = wgma::post::OrthoSol<wgma::post::SingleSpaceIntegrator>(cmesh, matids, conj,nThreads);
    //orthogonalise the modes
    auto normsol = ortho.Orthogonalise();
    //let us set the orthogonalised modes
    an.SetEigenvectors(normsol);
    an.LoadAllSolutions();
  }
  
}

void ComputeCouplingMat(wgma::wganalysis::WgmaPlanar &an,
                        std::string filename,
                        std::string matname)
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
  TPZVec<CSTATE> betavec = an.GetEigenvalues();
  for(auto &b : betavec){b = sqrt(b);}
  integrator.SetBeta(betavec);
  integrator.ComputeCoupling();
  TPZFMatrix<CSTATE> couplingmat;
  integrator.GetCoupling(couplingmat);
  std::ofstream matfile(filename);
  std::string name = matname+"=";
  couplingmat.Print(name.c_str(),matfile,EMathematicaInput);
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




#include <TPZParallelUtils.h>

void Extract1DSolFrom2DMesh(TPZAutoPointer<TPZCompMesh> mesh_1d,
                            TPZAutoPointer<TPZCompMesh> mesh_2d,
                            TPZFMatrix<CSTATE> &sol_1d)
{
  sol_1d.Zero();
  mesh_2d->LoadReferences();
  const TPZFMatrix<CSTATE> &sol_2d = mesh_2d->Solution();
  const auto &block_1d = mesh_1d->Block();
  const auto &block_2d = mesh_2d->Block();
  const int nel = mesh_1d->NElements();
  pzutils::ParallelFor(0,nel,[mesh_1d,&block_1d,&block_2d,
                              &sol_1d,&sol_2d](int iel){
    auto el_1d = mesh_1d->Element(iel);
    if(el_1d->Dimension() < mesh_1d->Dimension()){return;}
    //geometric mesh has 2d mesh as reference
    auto el_2d = el_1d->Reference()->Reference();
    if(!el_2d){DebugStop();}
    const auto ncon = el_2d->NConnects();
    if(ncon != el_1d->NConnects()){
      DebugStop();
    }
    for(auto icon = 0; icon < ncon; icon++){
      auto &con_2d = el_2d->Connect(icon);
      auto &con_1d = el_1d->Connect(icon);
        
      const int64_t dfseq_1d = con_1d.SequenceNumber();
      const int64_t pos_1d = block_1d.Position(dfseq_1d);

      const int64_t dfseq_2d = con_2d.SequenceNumber();
      const int64_t pos_2d = block_2d.Position(dfseq_2d);
        
      const int nshape = block_1d.Size(dfseq_1d);
      if(nshape!= block_2d.Size(dfseq_2d)){
        DebugStop();
      }
      for(int is = 0; is < nshape; is++){
        sol_1d.Put(pos_1d+is,0,sol_2d.Get(pos_2d+is,0));
      }
        
    }
  });
}

#include <materials/solutionprojection.hpp>
#include <TPZNullMaterial.h>

void ProjectSolIntoRestrictedMesh(wgma::wganalysis::WgmaPlanar &src_an,
                                  TPZAutoPointer<TPZCompMesh> scatt_mesh,
                                  const SimData &simdata,
                                  const TPZVec<int> &nmodes)
{
  std::map<int,STATE> error_proj;
  //1d mesh
  TPZAutoPointer<TPZCompMesh> src_mesh = src_an.GetMesh();
  //mesh used to project the solution into the restricted space
  TPZAutoPointer<TPZCompMesh> proj_mesh = src_mesh->Clone();
  //replace all materials for L2 projection
  for(auto [id,mat] : proj_mesh->MaterialVec()){
    auto bnd = dynamic_cast<TPZBndCond*>(mat);
    if(!bnd){
      auto newmat = new wgma::materials::SolutionProjection<CSTATE>(id,1);
      delete mat;
      proj_mesh->MaterialVec()[id]=newmat;
    }else{
      auto newmat = new TPZNullMaterial<CSTATE>(id,0);
      delete mat;
      proj_mesh->MaterialVec()[id]=newmat;
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
  std::cout<<"neqs without restriction: "<<sol_pml.Rows()<<std::endl;
  Extract1DSolFrom2DMesh(proj_mesh, scatt_mesh, sol_pml);
  //just to set size
  auto sol_proj = sol_pml;


  /*
    dirichlet boundary connects should not be restricted, otherwise
    this will result in all the equations on the same dependency
    being removed as well
   */
  std::set<int64_t> bound_connects;
  wgma::cmeshtools::FindDirichletConnects(proj_mesh, bound_connects);

  TPZVec<CSTATE> mode_norms;
  constexpr bool check_diag_mat{false};
  if(check_diag_mat){
    std::set<int> matids = {};
    constexpr bool conj{true};
    wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator> normcalc(src_mesh,matids,conj,simdata.nThreads);
    mode_norms = normcalc.ComputeNorm();
  }
  
  for(int im = 0; im < nmodes.size(); im++){
    const int nm = nmodes[im];
    //restrict dofs in proj mesh
    if(nm){
      RestrictDofs(proj_mesh, src_mesh, nm, bound_connects);
      TPZFMatrix<CSTATE> &sol = proj_mesh->Solution();
      std::cout<<"neqs after restriction: "<<sol.Rows()<<std::endl;
    }
    /*The TPZLinearAnalysis class manages the creation of the algebric
     * problem and the matrix inversion*/
    auto proj_an = wgma::scattering::Analysis(proj_mesh, simdata.nThreads,
                                              true,
                                              simdata.filterBoundEqs,false);
    

    //now we load desired solution into proj_mesh so we can project it
    {
      TPZFMatrix<CSTATE> &sol_reference = proj_mesh->Solution();
      Extract1DSolFrom2DMesh(proj_mesh, scatt_mesh, sol_reference);
    }
        
    proj_an.Run();

    /*now we get the matrix and check the diagonal:
     diag values should be phi_i*conj(phi_i)*/
    if(check_diag_mat){
      auto mat = proj_an.GetSolver().Matrix();

      for(int imode = 0; imode < nm; imode++){
        const auto matval = sqrt(mat->Get(imode,imode));
        const auto modenorm = mode_norms[imode];
        if(!IsZero(matval-modenorm)){
          std::cout<<"matrix error:\t";
          std::cout<<"nm "<<nm<<" imode "<<imode<<std::endl;
          std::cout<<"\tmatval : "<<matval<<"\tmodenorm: "<<modenorm<<std::endl;
        }
      }
    }
    
    //plot
    vtk.Do();
    //now we compute the error
    Extract1DSolFrom2DMesh(error_mesh, proj_mesh, sol_proj);
    // if(nm && nm<30){
    //   std::cout<<"nm : "<<nm<<std::endl;
    //   auto mat = proj_an.GetSolver().Matrix();
    //   std::cout<<"mat: "<<std::endl;
    //   for (int i = 0; i < mat->Rows(); i++)
    //   {
    //     for (int j = 0; j < mat->Cols(); j++)
    //     {
    //       std::cout<<'\t'<<mat->Get(i,j);
    //     }
    //     std::cout<<std::endl;
    //   }
    //   TPZFMatrix<CSTATE> &an_rhs = proj_an.Rhs();
    //   std::cout<<"rhs: "<<std::endl;
    //   for (int i = 0; i < an_rhs.Rows(); i++)
    //   {
    //     std::cout<<'\t'<<an_rhs.Get(i,0)<<std::endl;
    //   }
    //   TPZFMatrix<CSTATE> &an_sol = proj_an.Solution();
    //   std::cout<<"solution: "<<std::endl;
    //   for (int i = 0; i < an_sol.Rows(); i++)
    //   {
    //     std::cout<<'\t'<<an_sol.Get(i,0)<<std::endl;
    //   }
    // }
    
    sol_proj -= sol_pml;

    error_mesh->LoadSolution(sol_proj);
    vtk_error.Do();
    auto normsol =
      wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
    normsol.SetNThreads(std::thread::hardware_concurrency());
    const auto norm = std::real(normsol.ComputeNorm()[0]);
    std::cout<<"nmodes "<<nm<<" error "<<norm<<std::endl;
    error_proj.insert({nm,norm});
    wgma::cmeshtools::RemovePeriodicity(proj_mesh);
    proj_mesh->ComputeNodElCon();
    proj_mesh->CleanUpUnconnectedNodes();
  }

  std::cout<<"**********PROJECTION**********"<<std::endl;
  for(auto [nm,error] : error_proj){
    std::cout<<"nmodes "<<nm<<" norm error: "<<error<<std::endl;
  }
  std::cout<<"******************************"<<std::endl;
}

void SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
                     wgma::wganalysis::WgmaPlanar &src_an,
                     wgma::wganalysis::WgmaPlanar &match_an,
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
      //now for 1d pml mats
      for(const auto &pml : scatt_data.pmlvec){
        const std::string pattern_left{"source_clad_left"};
        const auto rx_left =
          std::regex{pattern_left, std::regex_constants::icase };
        const std::string pattern_right{"source_clad_right"};
        const auto rx_right =
          std::regex{pattern_right, std::regex_constants::icase };
        
        const bool found_pattern =
          std::regex_search(*(pml->names.begin()), rx_left) ||
          std::regex_search(*(pml->names.begin()), rx_right);

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
  TPZVec<int> nmodes = {0,1,2,5,10,15,20,30};//,50,100,200,400};
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
  ComputeWgbcCoeffs(src_an, src_data.wgbc_k, src_data.wgbc_f, false, src_coeffs);

  WgbcData match_data;
  match_data.cmesh = match_an.GetMesh();
  ComputeWgbcCoeffs(match_an, match_data.wgbc_k, match_data.wgbc_f,true, {});
    
    
  //here we will store the error between pml approx and wgbc approx
  TPZAutoPointer<TPZCompMesh> error_mesh = src_an.GetMesh()->Clone();
  TPZFMatrix<CSTATE> sol_pml = error_mesh->Solution();
  {
    Extract1DSolFrom2DMesh(error_mesh, scatt_mesh_pml, sol_pml);
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
    Extract1DSolFrom2DMesh(error_mesh, scatt_mesh_wgbc, sol_wgbc);
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

void SolveWithPML(TPZAutoPointer<TPZCompMesh> scatt_cmesh,
                  wgma::wganalysis::WgmaPlanar& src_an,
                  const TPZVec<CSTATE> &src_coeffs,
                  const SimData &simdata)
{
  auto scatt_an = wgma::scattering::Analysis(scatt_cmesh, simdata.nThreads,
                                             simdata.optimizeBandwidth,
                                             simdata.filterBoundEqs);
    
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

int64_t RestrictDofs(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                     TPZAutoPointer<TPZCompMesh> match_mesh,
                     const int nm,
                     const std::set<int64_t> bound_connects)
{
  //so we can access the dofs
  TPZFMatrix<CSTATE> &modal_sol = match_mesh->Solution();

  const int nmodes_total = modal_sol.Cols();
  if(nm > nmodes_total){
    std::cout<<"nm "<<nm<<" bigger than n of computed modes: "
             <<nmodes_total<<std::endl;
    DebugStop();
  }
  
  std::set<int64_t> dependent_connects;
  //this "dummy" connect will impose the restriction over the 1D line

  //TPZCompMesh::AllocateNewConnect(int nshape, int nstate, int order)
  const auto indep_con_id = scatt_mesh->AllocateNewConnect(nm,1,1);
  auto &indep_con = scatt_mesh->ConnectVec()[indep_con_id];
  //the geometric mesh will point to the scattering mesh
  scatt_mesh->LoadReferences();

      
  const auto &modal_block = match_mesh->Block();
  for(auto modal_el : match_mesh->ElementVec()){
    //select only 1d elements, skip boundary els
    if(modal_el->Dimension() < match_mesh->Dimension()){continue;}
    auto scatt_el = modal_el->Reference()->Reference();
    if(!scatt_el){DebugStop();}
    const auto ncon = scatt_el->NConnects();
    if(ncon != modal_el->NConnects()){
      DebugStop();
    }
    for(auto icon = 0; icon < ncon; icon++){
      auto &scatt_con = scatt_el->Connect(icon);
      auto &modal_con = modal_el->Connect(icon);
        
      const int64_t modal_dfseq = modal_con.SequenceNumber();
      const int64_t modal_pos = modal_block.Position(modal_dfseq);
      const int nshape = modal_block.Size(modal_dfseq);
      const auto dep_con_id = scatt_el->ConnectIndex(icon);
      if(dependent_connects.count(dep_con_id) == 0 &&
         bound_connects.count(dep_con_id) == 0 ){
        dependent_connects.insert(dep_con_id);
        //TPZDepend<TVar> *AddDependency(int64_t myindex, int64_t dependindex,TPZFMatrix<TVar> &depmat,int64_t ipos,int64_t jpos, int isize, int jsize);
        scatt_con.AddDependency(dep_con_id, indep_con_id, modal_sol, modal_pos, 0, nshape, nm);
      }
    }
  }
  scatt_mesh->ComputeNodElCon();
  scatt_mesh->CleanUpUnconnectedNodes();
  return indep_con_id;
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
  const auto pos = block.Position(seqnum);
  const auto sz = block.Size(seqnum);
  if(sz!=nm){DebugStop();}
  TPZManVector<int64_t,300> posvec_orig(sz,-1), posvec_filt(sz,-1);
  for(int i = 0; i < sz; i++){posvec_orig[i] = pos+i;};
  posvec_filt = posvec_orig;
  scatt_an.StructMatrix()->EquationFilter().Filter(posvec_orig, posvec_filt);
  for(auto imode = 0; imode < nm; imode++){
    const auto pos_f = posvec_orig[imode];
    const auto fval = fvec.Get(pos_f,0);
    fvec.Put(pos_f,0, fval+wgbc_f[imode]);
    for(auto jmode = 0; jmode < nm; jmode++){
      const auto ipos = posvec_filt[imode];
      const auto jpos = posvec_filt[jmode];
      const auto kval = mat->Get(ipos,jpos);
      mat->Put(ipos,jpos,kval+wgbc_k.Get(imode,jmode));
    }
  }
}


void ComputeWgbcCoeffs(wgma::wganalysis::WgmaPlanar& an,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       const bool positive_z, const TPZVec<CSTATE> &coeff){
  auto mesh = an.GetMesh();
  wgma::post::WaveguidePortBC<wgma::post::SingleSpaceIntegrator> wgbc(mesh);
  TPZManVector<CSTATE,1000> betavec = an.GetEigenvalues();
  for(auto &b : betavec){b = std::sqrt(b);}
  wgbc.SetPositiveZ(positive_z);
  if(coeff.size()){
    wgbc.SetSrcCoeff(coeff);
  }
  wgbc.SetBeta(betavec);
  wgbc.ComputeContribution();
  wgbc.GetContribution(wgbc_k,wgbc_f);
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
    RestrictDofs(scatt_mesh, match_mesh, nmodes, boundConnects);
  const int64_t indep_con_id_src =
    RestrictDofs(scatt_mesh, src_mesh, nmodes, boundConnects);

  auto scatt_an = wgma::scattering::Analysis(scatt_mesh, simdata.nThreads,
                                             simdata.optimizeBandwidth,
                                             simdata.filterBoundEqs);

  
  std::cout<<"nmodes on outgoing boundary: "<<nmodes<<std::endl;
    
  scatt_an.Assemble();
  //now we must add the waveguide port terms
  AddWaveguidePortContribution(scatt_an, indep_con_id_match,
                               nmodes, match_data.wgbc_k, match_data.wgbc_f);
  AddWaveguidePortContribution(scatt_an, indep_con_id_src,
                               nmodes, src_data.wgbc_k, src_data.wgbc_f);
  scatt_an.Solve();
}

