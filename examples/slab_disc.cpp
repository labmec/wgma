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
#include <post/waveguideportbc.hpp>
// pz includes
#include <MMeshType.h>      //for MMeshType
#include <TPZSimpleTimer.h> //for TPZSimpleTimer
#include <pzlog.h>          //for TPZLogger
#include <TPZMaterial.h>
#include <TPZVTKGenerator.h>
#include <pzinterpolationspace.h>

#include <regex>//for string search


constexpr bool optimizeBandwidth{true};
constexpr bool filterBoundaryEqs{true};
constexpr int nThreads{2};
constexpr int nThreadsDebug{0};
constexpr int vtkRes{0};
constexpr bool extendRightDomain{false};
constexpr bool extendLeftDomain{false};

TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const CSTATE target, const int neigenpairs,
            TPZEigenSort sorting, bool usingSLEPC);

void ComputeModes(wgma::wganalysis::WgmaPlanar &an,
                  TPZAutoPointer<TPZCompMesh> mesh,
                  const bool orthogonalise,
                  const int nThreads);

void PostProcessModes(wgma::wganalysis::WgmaPlanar &an,
                      TPZAutoPointer<TPZCompMesh> cmesh,
                      std::string filename,
                      const int vtkres);

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          wgma::wganalysis::WgmaPlanar& src_an,
                          wgma::wganalysis::WgmaPlanar& match_an,
                          const TPZVec<CSTATE> &source_coeffs,
                          const TPZVec<int> &nmodes,
                          const std::string &prefix);

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
  constexpr STATE lambda{1.55};

  
  constexpr STATE ncore{1.55};
  constexpr STATE nclad{1.000};

  constexpr STATE alphaPMLx {0.75};
  constexpr STATE alphaPMLy {0.75};
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

  // how to sort eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  constexpr bool usingSLEPC {true};

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


  constexpr int nEigenpairs_left{5};
  constexpr int nEigenpairs_right{5};

  constexpr CSTATE target{ncore*ncore};

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
  std::map<int64_t,int64_t> periodic_els;
  auto gmesh = wgma::gmeshtools::ReadPeriodicGmshMesh(meshfile, scale,gmshmats,
                                                      periodic_els,verbosity_lvl);

  // print wgma_gmesh to .txt and .vtk format
  if (printGMesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }


  /********************************
   * cmesh(modal analysis: left)  *
   ********************************/

  auto modal_l_cmesh = [gmesh,&gmshmats,scale](){
    // setting up cmesh data
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
    modal_l_an(modal_l_cmesh, nThreads,
             optimizeBandwidth, filterBoundaryEqs);
  modal_l_an.SetSolver(*solver_left);

  std::string modal_left_file{prefix+"_modal_left"};
  //no need to orthogonalise modes
  constexpr bool ortho{false};
  ComputeModes(modal_l_an, modal_l_cmesh, ortho, nThreads);
  if(exportVtk){
    PostProcessModes(modal_l_an, modal_l_cmesh, modal_left_file, vtkRes);
  }

  /********************************
   * cmesh(modal analysis: right)  *
   ********************************/

  auto modal_r_cmesh = [gmesh,&gmshmats,scale](){
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
    modal_r_an(modal_r_cmesh, nThreads,
             optimizeBandwidth, filterBoundaryEqs);
  modal_r_an.SetSolver(*solver_right);

  std::string modal_right_file{prefix+"_modal_right"};
  
  ComputeModes(modal_r_an, modal_r_cmesh, ortho, nThreads);
  if(exportVtk){
    std::string modal_left_file{prefix+"_modal_right"};
    PostProcessModes(modal_r_an, modal_r_cmesh, modal_right_file, vtkRes);
  }
  //now we need to project the scattering solution over the computed modes

  
  
  auto scatt_cmesh = [gmesh,&gmshmats](){
    
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData scatt_data;
    std::map<std::string, std::pair<CSTATE, CSTATE>> scatt_mats;
    scatt_mats["core_left"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    scatt_mats["core_right"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    scatt_mats["cladding_left"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad,1.);
    scatt_mats["cladding_right"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad,1.);
    std::map<std::string, wgma::bc::type> scatt_bcs;
    scatt_bcs["scatt_bnd_mid"] = wgma::bc::type::PEC;
    if constexpr(extendRightDomain){
      scatt_bcs["scatt_bnd_right"] = wgma::bc::type::PEC;
    }
    if constexpr(extendLeftDomain){
      scatt_bcs["scatt_bnd_left"] = wgma::bc::type::PEC;
    }

    
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, scatt_mats, scatt_bcs,
                                            {alphaPMLx,alphaPMLy}, scatt_data);


    //materials that will represent our source
    std::set<int> src_ids;
    if constexpr(extendLeftDomain){
      const std::string srcMat[] = {"source_core_left", "source_clad_left"};
      for(const auto &mat : srcMat){
        src_ids.insert(gmshmats[1].at(mat));
      }
      //now we add the pml mats to source mats

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
    


    //materials in which we would like to evaluate the solution
    std::vector<std::string> probeMats = {"source_clad_right", "source_core_right"};
    if(!extendLeftDomain){
      probeMats.push_back("source_clad_left");
      probeMats.push_back("source_core_left");
    }
    
    for(auto &mat : probeMats){
      const auto matdim = 1;
      const auto id = gmshmats[matdim].at(mat);
      scatt_data.probevec.push_back({id,matdim});
    }

    for(const auto &pml : scatt_data.pmlvec){
      const std::string pattern{"source_clad_right"};
      const auto rx = std::regex{pattern, std::regex_constants::icase };
    
      bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
      if constexpr(!extendLeftDomain){
        const std::string pattern{"source_clad_left"};
        const auto rx = std::regex{pattern, std::regex_constants::icase };
    
        found_pattern = found_pattern || std::regex_search(*(pml->names.begin()), rx);
      }
      if(found_pattern){
        const auto matdim = 1;
        const auto id = gmshmats[matdim].at(*pml->names.begin());
        std::cout<<"id "<<id<<" name "<<*pml->names.begin()<<std::endl;
        scatt_data.probevec.push_back({id,matdim});
      }
    }
    
    
    if constexpr(!extendRightDomain || !extendLeftDomain){
      std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
      for(const auto &pml : scatt_data.pmlvec){
        const std::string pattern_left{"xm"};
        const auto rx_left = std::regex{pattern_left, std::regex_constants::icase };
    
        const bool found_left = std::regex_search(*(pml->names.begin()), rx_left);

        const std::string pattern_right{"xp"};
        const auto rx_right = std::regex{pattern_right, std::regex_constants::icase };
    
        const bool found_right = std::regex_search(*(pml->names.begin()), rx_right);
        
        if((found_right && !extendRightDomain) ||(found_left && !extendLeftDomain)){
          continue;
        }
        pmlvec.push_back(pml);
      }

      
      scatt_data.pmlvec = pmlvec;
    }
    
    
    return wgma::scattering::CMeshScattering2D(gmesh, mode, pOrder, scatt_data,src_ids,
                                               lambda,scale);
  }();
  /*********************
   * solve(scattering) *  
   *********************/  
  TPZSimpleTimer tscatt("Scattering");


  //index of the mode to be used as a source (left wg)
  TPZVec<CSTATE> src_index = {0.,0.,1.,0.,0.};
  //index of the number of modes to be used to restrict the dofs of the scatt mesh(right wg)
  TPZVec<int> nmodes = {5};

  
  RestrictDofsAndSolve(scatt_cmesh, modal_l_an, modal_r_an,
                       src_index, nmodes, prefix);
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
                  TPZAutoPointer<TPZCompMesh> cmesh,
                  const bool orthogonalise,
                  const int nThreads)
{
  
  TPZSimpleTimer analysis("Modal analysis");
  an.Assemble();
  static constexpr bool computeVectors{true};
  an.Solve(computeVectors);
  //load all obtained modes into the mesh
  an.LoadAllSolutions();
  if(orthogonalise){
    //leave empty for all valid matids
    std::set<int> matids {};
    auto ortho = wgma::post::OrthoSol<wgma::post::SingleSpaceIntegrator>(cmesh, matids, nThreads);
    //orthogonalise the modes
    auto normsol = ortho.Orthogonalise();
    //let us set the orthogonalised modes
    an.SetEigenvectors(normsol);
    an.LoadAllSolutions();
  }
  
}

void PostProcessModes(wgma::wganalysis::WgmaPlanar &an,
                      TPZAutoPointer<TPZCompMesh> cmesh,
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
  auto vtk = TPZVTKGenerator(cmesh, fvars, file, vtkres);
  const auto nsol = an.GetEigenvectors().Cols();
  std::cout<<"Exporting "<<nsol<<" solutions"<<std::endl;
  for(auto isol = 0; isol < nsol; isol++){
    an.LoadSolution(isol);
    vtk.Do();
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
                                  const TPZVec<CSTATE> &wgbc_k,
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
    const auto pos_k = posvec_filt[imode];
    const auto kval = mat->Get(pos_k,pos_k);
    mat->Put(pos_k,pos_k, kval+wgbc_k[imode]);
    const auto pos_f = posvec_orig[imode];
    const auto fval = fvec.Get(pos_f,0);
    fvec.Put(pos_f,0, fval+wgbc_f[imode]);
  }
}

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          wgma::wganalysis::WgmaPlanar& src_an,
                          wgma::wganalysis::WgmaPlanar& match_an,
                          const TPZVec<CSTATE> &source_coeffs,
                          const TPZVec<int> &nmodes,
                          const std::string &prefix)
{
  TPZFMatrix<CSTATE> ev = match_an.GetEigenvectors();
  const auto evr = ev.Rows();
  const auto evc = ev.Cols();

  auto gmesh = scatt_mesh->Reference();
  
  //let us load all the modes into the comp mesh
  match_an.LoadAllSolutions();
  src_an.LoadAllSolutions();

  auto match_mesh = match_an.GetMesh();
  auto src_mesh = src_an.GetMesh();
  //first we compute the waveguide port bc values for the match mesh
  TPZVec<CSTATE> wgbc_k_right, wgbc_f_right;
  {
    wgma::post::WaveguidePortBC<wgma::post::SingleSpaceIntegrator> wgbc(match_mesh);
    TPZVec<CSTATE> betavec = match_an.GetEigenvalues();
    for(auto &b : betavec){b = std::sqrt(b);}
    wgbc.SetPositiveZ(true);
    wgbc.SetBeta(betavec);
    wgbc.ComputeContribution();
    wgbc.GetContribution(wgbc_k_right,wgbc_f_right);
  }

  TPZVec<CSTATE> wgbc_k_left,wgbc_f_left;
  //now we compute the waveguide port bc values for the source mesh
  if(!extendLeftDomain){
    wgma::post::WaveguidePortBC<wgma::post::SingleSpaceIntegrator> wgbc(src_mesh);
    TPZVec<CSTATE> betavec = src_an.GetEigenvalues();
    for(auto &b : betavec){b = std::sqrt(b);}
    wgbc.SetSrcCoeff(source_coeffs);
    wgbc.SetPositiveZ(false);
    wgbc.SetBeta(betavec);
    wgbc.ComputeContribution();
    wgbc.GetContribution(wgbc_k_left,wgbc_f_left);
  }
  

  //set up post processing
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_imag",
    "Field_abs"};
  const std::string scatt_file = prefix+"_scatt";
  auto vtk = TPZVTKGenerator(scatt_mesh, fvars, scatt_file, vtkRes);

  /**
   **/
  std::set<int64_t> boundConnects;
  wgma::cmeshtools::FindDirichletConnects(scatt_mesh, boundConnects);

  int64_t indep_con_id_right{-1}, indep_con_id_left{-1};
  for(auto nm : nmodes){
    if(nm){
      indep_con_id_right = RestrictDofs(scatt_mesh, match_mesh, nm, boundConnects);
      if constexpr (!extendLeftDomain){
        indep_con_id_left = RestrictDofs(scatt_mesh, src_mesh, nm, boundConnects);
      }
    }

    auto scatt_an = wgma::scattering::Analysis(scatt_mesh, nThreads,
                                               optimizeBandwidth, filterBoundaryEqs);

  
    std::cout<<"nmodes on outgoing boundary: "<<nm<<std::endl;
    if(extendLeftDomain){
      auto src_mesh = src_an.GetMesh();
      //get id of source materials
      wgma::scattering::SourceWgma src;
      for(auto [id,mat] : src_mesh->MaterialVec()){
        src.id.insert(id);
      }
      src.modal_cmesh = src_mesh;

      TPZFMatrix<CSTATE> sol(scatt_mesh->NEquations(),1);
      const int nsol = source_coeffs.size();
      for(int isol = 0; isol < nsol; isol++){
        src_an.LoadSolution(isol);
        auto beta = std::sqrt(src_an.GetEigenvalues()[isol]);
        wgma::scattering::LoadSource1D(scatt_mesh, src);
        wgma::scattering::SetPropagationConstant(scatt_mesh, beta*source_coeffs[isol]);
        if(isol == 0){
          scatt_an.Assemble();
          if(nm){
            AddWaveguidePortContribution(scatt_an, indep_con_id_right, nm, wgbc_k_right, wgbc_f_right);
          }
        }else{
          scatt_an.AssembleRhs(src.id);
          TPZFMatrix<CSTATE> &rhs = scatt_an.Rhs();
          std::cout<<"rhs norm "<<Norm(rhs)<<std::endl;
        }
        scatt_an.Solve();
        scatt_an.LoadSolution();
        TPZFMatrix<CSTATE> &curr_sol = scatt_an.Solution();
        sol+=curr_sol;
      }
      scatt_an.LoadSolution(sol);
    }else{
      scatt_an.Assemble();
      //now we must add the waveguide port terms
      if(nm){
        AddWaveguidePortContribution(scatt_an, indep_con_id_right, nm, wgbc_k_right, wgbc_f_right);
      }
      AddWaveguidePortContribution(scatt_an, indep_con_id_left, nm, wgbc_k_left, wgbc_f_left);
      scatt_an.Solve();
    }

    vtk.Do();

    if(nm){
      wgma::cmeshtools::RemovePeriodicity(scatt_mesh);
    }    
    scatt_mesh->ComputeNodElCon();
    scatt_mesh->CleanUpUnconnectedNodes();
  }
}