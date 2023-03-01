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
#include <pzinterpolationspace.h>

#include <regex>//for string search


constexpr bool optimizeBandwidth{false};
constexpr bool filterBoundaryEqs{true};
constexpr int nThreads{8};
constexpr int nThreadsDebug{0};
constexpr int vtkRes{2};


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
                          const std::vector<int> &sources,
                          const std::set<int> &source_mats,
                          const std::vector<int> &nmodes,
                          const TPZVec<std::map<std::string, int>> &gmshmats,
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
  constexpr STATE lambda{1.5};

  
  constexpr STATE ncore{1.5};
  constexpr STATE nclad{1.000};

  constexpr STATE alphaPMLx {1.5};
  constexpr STATE alphaPMLy {1.0};
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


  constexpr int nEigenpairs_left{1};
  constexpr int nEigenpairs_right{50};

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
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs, alphaPMLx,
                                            alphaPMLy, modal_data, modal_dim);
    //we must now filter the 1D PMLs
    std::vector<wgma::pml::data>  pmlvec;
    for(const auto &pml : modal_data.pmlvec){
      const std::string pattern{"source_clad_left"};
      const auto rx = std::regex{pattern, std::regex_constants::icase };
    
      const bool found_pattern = std::regex_search(*(pml.names.begin()), rx);
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
  bool ortho{false};
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
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs, alphaPMLx,
                                            alphaPMLy, modal_data, modal_dim);

    //we must now filter the 1D PMLs
    std::vector<wgma::pml::data>  pmlvec;
    for(const auto &pml : modal_data.pmlvec){
      const std::string pattern{"source_clad_right"};
      const auto rx = std::regex{pattern, std::regex_constants::icase };
    
      const bool found_pattern = std::regex_search(*(pml.names.begin()), rx);
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
  //we will orthogonalise the modes for projecting the scattering solution in them
  ortho = true;
  
  ComputeModes(modal_r_an, modal_r_cmesh, ortho, nThreads);
  if(exportVtk){
    std::string modal_left_file{prefix+"_modal_right"};
    PostProcessModes(modal_r_an, modal_r_cmesh, modal_right_file, vtkRes);
  }
  //now we need to project the scattering solution over the computed modes

  //materials that will represent our source
  std::set<int> src_ids;
  const std::string srcMat[] = {"source_core_left", "source_clad_left"};
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
    scatt_bcs["scatt_bnd"] = wgma::bc::type::PEC;

    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, scatt_mats, scatt_bcs, alphaPMLx,
                                            alphaPMLy, scatt_data);
    
    //materials in which we would like to evaluate the solution
    const std::string probeMats[] = {"source_clad_right", "source_core_right"};

    for(auto &mat : probeMats){
      const auto matdim = 1;
      const auto id = gmshmats[matdim].at(mat);
      scatt_data.probevec.push_back({id,matdim});
    }
    
    return wgma::scattering::CMeshScattering2D(gmesh, mode, pOrder, scatt_data,src_ids,
                                               lambda,scale);
  }();
  /*********************
   * solve(scattering) *  
   *********************/  
  TPZSimpleTimer tscatt("Scattering");


  //index of the mode to be used as a source (left wg)
  std::vector<int> src_index = {0};
  //index of the number of modes to be used to restrict the dofs of the scatt mesh(right wg)
  std::vector<int> nmodes = {0,1,2,5,10,15,20,30};
  
  RestrictDofsAndSolve(scatt_cmesh, modal_l_an, modal_r_an,
                       src_index, src_ids, nmodes,gmshmats, prefix);
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
  const int krylovDim = std::max(5,2*neigenpairs);
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

void ComputeModes(wgma::wganalysis::WgmaPlanar &an,
                  TPZAutoPointer<TPZCompMesh> cmesh,
                  const bool orthogonalise,
                  const int nThreads)
{
  
  TPZSimpleTimer analysis("Modal analysis");
  an.Assemble(TPZEigenAnalysis::Mat::A);
  an.Assemble(TPZEigenAnalysis::Mat::B);
  static constexpr bool computeVectors{true};
  an.Solve(computeVectors);
  //let us discard all solutions with negative real part of beta

  auto &ev = an.GetEigenvalues();
  auto &evec = an.GetEigenvectors();

  int count = 0;
  for(auto w : ev){
    if(std::real(w) <= 0){break;}
    count++;
  }

  ev.Resize(count);
  const int neq = evec.Rows();
  evec.Resize(neq, count);
  //load all obtained modes into the mesh
  an.LoadAllSolutions();
  //leave empty for all valid matids
  std::set<int> matids {};
  auto ortho = wgma::post::OrthoSol<wgma::post::SingleSpaceIntegrator>(cmesh, matids, nThreads);
  //orthogonalise the modes
  auto normsol = ortho.Orthogonalise();
  //let us set the orthogonalised modes
  an.SetEigenvectors(normsol);

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

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          wgma::wganalysis::WgmaPlanar& src_an,
                          wgma::wganalysis::WgmaPlanar& match_an,
                          const std::vector<int> &sources,
                          const std::set<int> &source_mats,
                          const std::vector<int> &nmodes,
                          const TPZVec<std::map<std::string, int>> &gmshmats,
                          const std::string &prefix)
{
  auto ev = match_an.GetEigenvectors();
  const auto evr = ev.Rows();
  const auto evc = ev.Cols();

  auto gmesh = scatt_mesh->Reference();
  
  auto match_mesh = match_an.GetMesh();

  //set up post processing
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_imag",
    "Field_abs"};
  const std::string scatt_file = prefix+"_scatt";
  auto vtk = TPZVTKGenerator(scatt_mesh, fvars, scatt_file, vtkRes);
  
  for(auto nm : nmodes){

    if(nm){
      if(nm>evc){
        std::cout<<"nm "<<nm<<" bigger than n of computed modes "
                 <<evc<<" skipping.."<<std::endl;
        continue;
      }
      const TPZFMatrix<CSTATE> new_ev(evr,nm,ev.Elem(),evr*evc);
      match_an.SetEigenvectors(new_ev);
      //this "dummy" connect will impose the restriction over the 1D line
      auto indep_con_id = scatt_mesh->AllocateNewConnect(nm,1,1);
      auto &indep_con = scatt_mesh->ConnectVec()[indep_con_id];
      indep_con.IncrementElConnected();
      
      /*we will iterate over the elements of the modal analysis mesh
        and find their counterpart on the scattering mesh*/
      const auto nshape = indep_con.NShape();
      //we load all the modes in the comp mesh
      match_an.LoadAllSolutions();
      //the geometric mesh will point to the scattering mesh
      gmesh->ResetReference();
      scatt_mesh->LoadReferences();

      //so we can access the dofs
      TPZFMatrix<CSTATE> &modal_sol = match_mesh->Solution();
      auto &modal_block = match_mesh->Block();
      for(auto modal_el : match_mesh->ElementVec()){
        if(modal_el->Dimension()< match_mesh->Dimension()){continue;}
        auto scatt_el = modal_el->Reference()->Reference();
        if(!scatt_el){continue;}
        const auto ncon = scatt_el->NConnects();
        if(ncon != modal_el->NConnects()){
          DebugStop();
        }
        for(auto icon = 0; icon < ncon; icon++){
          auto &scatt_con = scatt_el->Connect(icon);
          auto &modal_con = modal_el->Connect(icon);
        
          const int64_t modal_dfseq = modal_con.SequenceNumber();
          const int64_t modal_pos = modal_block.Position(modal_dfseq);
          const int modal_dfvar = modal_block.Size(modal_dfseq);
          scatt_con.AddDependency(scatt_el->ConnectIndex(icon), indep_con_id, modal_sol, modal_pos, 0, modal_dfvar, nshape);
        }
      }
      scatt_mesh->ComputeNodElCon();
      scatt_mesh->CleanUpUnconnectedNodes();
    }

    auto scatt_an = wgma::scattering::Analysis(scatt_mesh, nThreads,
                                               optimizeBandwidth, filterBoundaryEqs);

    //in the first run we assemble the whole algebraic system. afterwards we can
    //just compute the rhs
    bool firstrun=true;

    /*********************
     * cmesh(scattering) *
     *********************/

    auto src_mesh = src_an.GetMesh();
    //get id of source materials
    wgma::scattering::SourceWgma src;
    src.id = source_mats;
    src.modal_cmesh = src_mesh;
    const int nsols = sources.size();
  
    for(int isol = 0; isol < nsols; isol++){
      std::cout<<"running source "<<isol+1<<" out of "<<nsols<<std::endl;
      src_an.LoadSolution(sources[isol]);
      auto beta = std::sqrt(src_an.GetEigenvalues()[isol]);
      wgma::scattering::LoadSource1D(scatt_mesh, src);
      wgma::scattering::SetPropagationConstant(scatt_mesh, beta);
      
      if(firstrun){
        firstrun=false;
        scatt_an.Assemble();
      }else{
        scatt_an.AssembleRhs(source_mats);
      }
      scatt_an.Solve();
      vtk.Do();
    }

    wgma::cmeshtools::RemovePeriodicity(scatt_mesh);
    scatt_mesh->ComputeNodElCon();
    scatt_mesh->CleanUpUnconnectedNodes();
  }

  match_an.SetEigenvectors(ev);
}