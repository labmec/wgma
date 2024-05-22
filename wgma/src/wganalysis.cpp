#include "wganalysis.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "slepcepshandler.hpp"
#include "cmeshtools_impl.hpp"

#include <TPZCutHillMcKee.h>
#include <TPZSpStructMatrix.h>
#include <TPZStructMatrixOMPorTBB.h>
#include <TPZKrylovEigenSolver.h>
#include <TPZQuadEigenSolver.h>
#include <Electromagnetics/TPZWgma.h>
#include <Electromagnetics/TPZAnisoWgma.h>
#include <Electromagnetics/TPZPeriodicWgma.h>
#include <Electromagnetics/TPZPlanarWgma.h>
#include <materials/solutionprojection.hpp>
#include <TPZNullMaterialCS.h>
#include <TPZSimpleTimer.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzelctemp.h>

#include <cassert>

namespace wgma::wganalysis{
  bool using_tbb_mat{false};
  //! Default suggestion for setting up the eigensolver for modal analysis.
  TPZAutoPointer<TPZEigenSolver<CSTATE>>
  SetupSolver(const CSTATE target, const int nEigen,
              TPZEigenSort sorting, bool &usingSLEPC, int krylovDim, bool verbose)
  {
    
#ifndef WGMA_USING_SLEPC
    if(usingSLEPC){
      std::cout<<"wgma was not configured with slepc. defaulting to: "
               <<"TPZKrylovSolver"<<std::endl;
      usingSLEPC = false;
    }
#endif

    TPZAutoPointer<TPZEigenSolver<CSTATE>> solver{nullptr};
    if (usingSLEPC){
      using namespace ::wgma::slepc;
      /*
        The following are suggested SLEPc settings.
        NOTE: -1 stands for PETSC_DECIDE
      */
    
      constexpr STATE eps_tol = 1e-13;//PETSC_DECIDE
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
      const bool eps_verbosity = verbose;
    
    
      auto eps_solver = new EPSHandler<CSTATE>;
      eps_solver->SetType(eps_solver_type);
      eps_solver->SetProblemType(eps_prob_type);
      eps_solver->SetEPSDimensions(nEigen, krylovDim, eps_mpd);
      eps_solver->SetTarget(target);
      eps_solver->SetTolerances(eps_tol,eps_max_its);
      eps_solver->SetConvergenceTest(eps_conv_test);
      eps_solver->SetKrylovOptions(eps_krylov_locking,eps_krylov_restart);
      eps_solver->SetVerbose(eps_verbosity);
      eps_solver->SetTrueResidual(eps_true_residual);
    
      eps_solver->SetLinearSolver(linsolver);
      eps_solver->SetLinearSolverTol(ksp_rtol,ksp_atol,ksp_dtol,ksp_max_its);
      eps_solver->SetPrecond(pc, 1e-14);

      solver = eps_solver;
    }else{
      krylovDim = krylovDim < nEigen ? 5*nEigen : krylovDim;
      auto krylov_solver = new TPZKrylovEigenSolver<CSTATE>;
      TPZSTShiftAndInvert<CSTATE> st;
      krylov_solver->SetSpectralTransform(st);
      krylov_solver->SetTarget(target);
      krylov_solver->SetKrylovDim(krylovDim);
      krylov_solver->SetNEigenpairs(nEigen);
      krylov_solver->SetAsGeneralised(true);
    
      solver = krylov_solver;
    }
  
    solver->SetEigenSorting(sorting);
    return solver;
  }
  
  Wgma::~Wgma(){}

  void Wgma::Run(bool compute_eigenvectors){

    std::cout<<"Assembling..."<<std::flush;
    {
      TPZSimpleTimer assemble("Assemble");
      Assemble();
    }
    std::cout<<"\rAssembled!"<<std::endl;
    
    TPZEigenSolver<CSTATE> *solv =
      dynamic_cast<TPZEigenSolver<CSTATE>*>(this->Solver());
    if(!solv){
      std::cerr<<__PRETTY_FUNCTION__
               <<"\nA solver has not been set.\n"
               <<"Check documentation of TPZKrylovEigenSolver"
               <<"\nor wgma::slepc::EPSHandler\nAborting..."<<std::endl;
      exit(-1);
    }

    AdjustSolver(solv);
    Solve(compute_eigenvectors);
    
  }
  
  void Wgma::Solve(bool compute_eigenvectors, bool verbose){
    TPZEigenAnalysisBase::SetComputeEigenvectors(compute_eigenvectors);
    TPZEigenSolver<CSTATE> *solv =
      dynamic_cast<TPZEigenSolver<CSTATE>*>(this->Solver());
    if(!solv){
      std::cerr<<__PRETTY_FUNCTION__
               <<"\nA solver has not been set.\n"
               <<"Check documentation of TPZKrylovEigenSolver"
               <<"\nor wgma::slepc::EPSHandler\nAborting..."<<std::endl;
      exit(-1);
    }


    TPZSimpleTimer tsolv("Solve");
    std::cout<<"Solving..."<<std::flush;
    this->Solve();
    std::cout<<"\rSolved!"<<std::endl;

    if(verbose){
      std::ios cout_state(nullptr);
      cout_state.copyfmt(std::cout);
    
      std::cout << std::setprecision(std::numeric_limits<STATE>::max_digits10);
      for(auto &w : this->GetEigenvalues()){
        std::cout<<w<<std::endl;
      }
      std::cout.copyfmt(cout_state);
    }
  }

  void Wgma::LoadSolution(const int isol)
  {
    const auto &ev = this->GetEigenvalues();
    const auto &eigenvectors = this->GetEigenvectors();
    
    if(eigenvectors.Rows() == 0){
      return;
    }
    if(isol > ev.size()){
      std::cout<<__PRETTY_FUNCTION__
               <<"\nThe solution "<<isol<<" was requested but\n"
               <<"only  "<<ev.size()<<" eigenvectors\n"
               <<"were calculated. Aborting...\n";
      DebugStop();
    }

    LoadSolutionInternal(isol, 1);
  }


  void Wgma::LoadAllSolutions()
  {
    const auto nsol = this->GetEigenvectors().Cols();
    LoadSolutionInternal(0, nsol);
  }
  

  Wgma2D::Wgma2D(const TPZVec<TPZAutoPointer<TPZCompMesh>> &meshvec,
                 const int n_threads, const bool reorder_eqs,
                 const bool filter_bound)
  {
    
    m_filter_bound = filter_bound;
    if(meshvec.size() != 3){
      std::cerr<<__PRETTY_FUNCTION__
               <<"\nThree computational meshes are required."
               <<"Aborting...\n";
      exit(-1);
    }
    //gets the multiphysics mesh (main mesh)
    m_cmesh_mf = meshvec[0];
    m_cmesh_h1 = meshvec[1 + TPZWgma::H1Index()];
    m_cmesh_hcurl = meshvec[1 + TPZWgma::HCurlIndex()];

    this->SetRenumber(new TPZCutHillMcKee());
    this->SetCompMeshInit(m_cmesh_mf.operator->(), reorder_eqs);

    TPZAutoPointer<TPZStructMatrix> strmtrx{nullptr};
    if(using_tbb_mat){
      auto mtrx = new TPZSpStructMatrix<CSTATE,TPZStructMatrixOMPorTBB<CSTATE>>(m_cmesh_mf);
      mtrx->SetShouldColor(false);
      mtrx->SetTBBorOMP(true);
      strmtrx = mtrx;
    }else{
      auto mtrx = new TPZSpStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>(m_cmesh_mf);
      strmtrx = mtrx;
    }

    strmtrx->SetNumThreads(n_threads);
    
  
    
    
    if(m_filter_bound){
      TPZVec<int64_t> activeEquations;
      int n_dofs_before = m_cmesh_mf->NEquations();
      wgma::cmeshtools::FilterBoundaryEquations(m_cmesh_mf, activeEquations,
                                                m_bound_cons);
      CountActiveEqs(m_n_dofs_mf,m_n_dofs_h1,m_n_dofs_hcurl);
      std::cout<<"neq(before): "<<n_dofs_before
               <<"\tneq(after): "<<m_n_dofs_mf<<std::endl;
      strmtrx->EquationFilter().SetActiveEquations(activeEquations);
    }else{
      CountActiveEqs(m_n_dofs_mf,m_n_dofs_h1,m_n_dofs_hcurl);
    }
    this->SetStructuralMatrix(strmtrx);
  }

  void Wgma2D::AdjustSolver(TPZEigenSolver<CSTATE> *solver){
    auto *krylov_solver =
      dynamic_cast<TPZKrylovEigenSolver<CSTATE>*> (solver);      
    if(krylov_solver && krylov_solver->KrylovInitialVector().Rows() == 0){
      /**this is to ensure that the eigenvector subspace is orthogonal to
         the spurious solutions associated with et = 0 ez != 0*/
      TPZFMatrix<CSTATE> initVec(m_n_dofs_mf, 1, 0.);
      const auto firstHCurl = m_n_dofs_h1 * TPZWgma::HCurlIndex();
      for (int i = 0; i < m_n_dofs_hcurl; i++) {
        initVec(firstHCurl + i, 0) = 1;
      }
      krylov_solver->SetKrylovInitialVector(initVec);
    }
  }
  
  void
  Wgma2D::CountActiveEqs(int &neq,int &nH1Equations, int &nHCurlEquations)
  {
    auto &cmesh = m_cmesh_mf;
    neq = nH1Equations = nHCurlEquations = 0;
    auto &cmeshHCurl = m_cmesh_hcurl;
    auto &cmeshH1 = m_cmesh_h1;
    auto &boundConnects = m_bound_cons;
    
    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
      bool isH1;
      if (boundConnects.find(iCon) == boundConnects.end()) {
        if (cmesh->ConnectVec()[iCon].HasDependency())
          continue;
        int seqnum = cmesh->ConnectVec()[iCon].SequenceNumber();
        int blocksize = cmesh->Block().Size(seqnum);
        if (TPZWgma::H1Index() == 0 && iCon < cmeshH1->NConnects()) {
          isH1 = true;
        } else if (TPZWgma::H1Index() == 1 && iCon >= cmeshHCurl->NConnects()) {
          isH1 = true;
        } else {
          isH1 = false;
        }
        for (int ieq = 0; ieq < blocksize; ieq++) {
          neq++;
          isH1 == true ? nH1Equations++ : nHCurlEquations++;
        }
      }
    }
    std::cout << "------\tactive eqs\t-------" << std::endl;
    std::cout << "# H1 equations: " << nH1Equations << std::endl;
    std::cout << "# HCurl equations: " << nHCurlEquations << std::endl;
    std::cout << "# equations: " << neq << std::endl;
    std::cout << "------\t----------\t-------" << std::endl;
    return;
  }
  
  void Wgma2D::LoadSolutionInternal(const int isol, const int ncols)
  { 
    const auto &ev = this->GetEigenvalues();
    const auto &eigenvectors = this->GetEigenvectors();

    
    
    TPZManVector<TPZAutoPointer<TPZCompMesh>,2> meshVecPost(2);
    meshVecPost[0] = m_cmesh_h1;
    meshVecPost[1] = m_cmesh_hcurl;

    //we dont have a proper way of setting multiple kz yet
    const CSTATE currentKz = [&ev,isol](){
      auto tmp = std::sqrt(-1.0*ev[isol]);
      constexpr auto epsilon = std::numeric_limits<STATE>::epsilon()/
        (10*std::numeric_limits<STATE>::digits10);
      //let us discard extremely small imag parts
      if (tmp.imag() < epsilon)
        {tmp = tmp.real();}
      return tmp;
    }();

    

    for(auto mat : m_cmesh_mf->MaterialVec()){
      auto id = mat.first;
      auto matPtr =
        dynamic_cast<TPZWgma *>(m_cmesh_mf->FindMaterial(id));
      if(!matPtr) continue;
      matPtr->SetKz(currentKz);
    }
    
    if(isol == 0 && ncols == eigenvectors.Cols()){
      this->LoadSolution(GetEigenvectors());
    }else{
      const auto neqOriginal = eigenvectors.Rows();
      TPZFMatrix<CSTATE> evector(neqOriginal, ncols, 0.);
      eigenvectors.GetSub(0, isol, neqOriginal, ncols, evector);
      this->LoadSolution(evector);
    }
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshVecPost, m_cmesh_mf);
  }
  
  void Wgma2D::WriteToCsv(std::string filename, STATE lambda){
    const auto &ev = this->GetEigenvalues();
    const int nev = ev.size();
    if(nev < 1){
      std::cout<<"There are no eigenvalues to write to .csv"<<std::endl;
      return;
    }
    std::cout<<"Exporting eigen info..."<<std::endl;
    std::ostringstream eigeninfo;
    constexpr auto pres = std::numeric_limits<STATE>::max_digits10;
    eigeninfo.precision(pres);
    REAL hSize = -1e12;
    TPZGeoMesh *gmesh = m_cmesh_mf->Reference();
    for(auto *el : gmesh->ElementVec()){
      if(el->HasSubElement()) continue;
      const auto elRadius = el->ElementRadius();
      hSize = elRadius > hSize ? elRadius : hSize;
    }
    const auto neq = m_cmesh_mf->NEquations();
    const auto nel = m_cmesh_mf->NElements();
    const auto porder = m_cmesh_mf->GetDefaultOrder();
    
    eigeninfo << std::fixed << neq << "," << nel << ",";
    eigeninfo << std::fixed << hSize << "," << porder<<",";
    eigeninfo << std::fixed << lambda<<",";
    eigeninfo << nev<<",";
    for(int i = 0; i < nev ; i++){
      eigeninfo<<std::fixed<<std::real(ev[i])<<",";
      eigeninfo<<std::fixed<<std::imag(ev[i]);
      if(i != nev - 1 ) {
        eigeninfo << ",";
      }
    }
    eigeninfo << std::endl;

    std::ofstream file(filename.c_str(),std::ios::app);
    file << eigeninfo.str();
    file.close();
  }


  WgmaAniso2D::WgmaAniso2D(const TPZVec<TPZAutoPointer<TPZCompMesh>> &meshvec,
                 const int n_threads, const bool reorder_eqs,
                 const bool filter_bound)
  {
    
    m_filter_bound = filter_bound;
    if(meshvec.size() != 3){
      std::cerr<<__PRETTY_FUNCTION__
               <<"\nThree computational meshes are required."
               <<"Aborting...\n";
      exit(-1);
    }
    //gets the multiphysics mesh (main mesh)
    m_cmesh_mf = meshvec[0];
    m_cmesh_h1 = meshvec[1 + TPZAnisoWgma::H1Index()];
    m_cmesh_hcurl = meshvec[1 + TPZAnisoWgma::HCurlIndex()];

    this->SetRenumber(new TPZCutHillMcKee());
    this->SetCompMeshInit(m_cmesh_mf.operator->(), reorder_eqs);

    TPZAutoPointer<TPZStructMatrix> strmtrx{nullptr};
    if(using_tbb_mat){
      auto mtrx = new TPZSpStructMatrix<CSTATE,TPZStructMatrixOMPorTBB<CSTATE>>(m_cmesh_mf);
      mtrx->SetShouldColor(false);
      mtrx->SetTBBorOMP(true);
      strmtrx = mtrx;
    }else{
      auto mtrx = new TPZSpStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>(m_cmesh_mf);
      strmtrx = mtrx;
    }
    

    strmtrx->SetNumThreads(n_threads);
    
  
    
    
    if(m_filter_bound){
      TPZVec<int64_t> activeEquations;
      int n_dofs_before = m_cmesh_mf->NEquations();
      wgma::cmeshtools::FilterBoundaryEquations(m_cmesh_mf, activeEquations,
                                                m_bound_cons);
      CountActiveEqs(m_n_dofs_mf,m_n_dofs_h1,m_n_dofs_hcurl);
      std::cout<<"neq(before): "<<n_dofs_before
               <<"\tneq(after): "<<m_n_dofs_mf<<std::endl;
      strmtrx->EquationFilter().SetActiveEquations(activeEquations);
    }else{
      CountActiveEqs(m_n_dofs_mf,m_n_dofs_h1,m_n_dofs_hcurl);
    }
    this->SetStructuralMatrix(strmtrx);
  }

  void WgmaAniso2D::AdjustSolver(TPZEigenSolver<CSTATE> *solver){
    auto *krylov_solver =
      dynamic_cast<TPZQuadEigenSolver<CSTATE>*> (solver);      
    if(krylov_solver && krylov_solver->KrylovInitialVector().Rows() == 0){
      /**this is to ensure that the eigenvector subspace is orthogonal to
         the spurious solutions associated with et = 0 ez != 0*/
      TPZFMatrix<CSTATE> initVec(2*m_n_dofs_mf, 1, 0.);
      const auto firstHCurl = m_n_dofs_h1 * TPZAnisoWgma::HCurlIndex();
      for (int i = 0; i < m_n_dofs_hcurl; i++) {
        initVec(firstHCurl + i, 0) = 1;
      }
      const auto target = krylov_solver->Target();
      for (int i = 0; i < m_n_dofs_hcurl; i++) {
        initVec(m_n_dofs_mf+firstHCurl + i, 0) = target;
      }
      
      krylov_solver->SetKrylovInitialVector(initVec);
    }
  }
  
  void
  WgmaAniso2D::CountActiveEqs(int &neq,int &nH1Equations, int &nHCurlEquations)
  {
    auto &cmesh = m_cmesh_mf;
    neq = nH1Equations = nHCurlEquations = 0;
    auto &cmeshHCurl = m_cmesh_hcurl;
    auto &cmeshH1 = m_cmesh_h1;
    auto &boundConnects = m_bound_cons;
    
    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
      bool isH1;
      if (boundConnects.find(iCon) == boundConnects.end()) {
        if (cmesh->ConnectVec()[iCon].HasDependency())
          continue;
        int seqnum = cmesh->ConnectVec()[iCon].SequenceNumber();
        int blocksize = cmesh->Block().Size(seqnum);
        if (TPZAnisoWgma::H1Index() == 0 && iCon < cmeshH1->NConnects()) {
          isH1 = true;
        } else if (TPZAnisoWgma::H1Index() == 1 && iCon >= cmeshHCurl->NConnects()) {
          isH1 = true;
        } else {
          isH1 = false;
        }
        for (int ieq = 0; ieq < blocksize; ieq++) {
          neq++;
          isH1 == true ? nH1Equations++ : nHCurlEquations++;
        }
      }
    }
    std::cout << "------\tactive eqs\t-------" << std::endl;
    std::cout << "# H1 equations: " << nH1Equations << std::endl;
    std::cout << "# HCurl equations: " << nHCurlEquations << std::endl;
    std::cout << "# equations: " << neq << std::endl;
    std::cout << "------\t----------\t-------" << std::endl;
    return;
  }
  
  void WgmaAniso2D::LoadSolutionInternal(const int isol, const int ncols)
  { 
    const auto &ev = this->GetEigenvalues();
    const auto &eigenvectors = this->GetEigenvectors();
    
    TPZManVector<TPZAutoPointer<TPZCompMesh>,2> meshVecPost(2);
    meshVecPost[0] = m_cmesh_h1;
    meshVecPost[1] = m_cmesh_hcurl;

    //we dont have a proper way of setting multiple kz yet
    const CSTATE currentKz = [&ev,isol](){
      auto tmp = CSTATE(-1.0*1i)*ev[isol];
      constexpr auto epsilon = std::numeric_limits<STATE>::epsilon()/
        (10*std::numeric_limits<STATE>::digits10);
      //let us discard extremely small imag parts
      if (tmp.imag() < epsilon)
        {tmp = tmp.real();}
      return tmp;
    }();

    

    
    for(auto mat : m_cmesh_mf->MaterialVec()){
      auto id = mat.first;
      auto matPtr =
        dynamic_cast<TPZWgma *>(m_cmesh_mf->FindMaterial(id));
      if(!matPtr) continue;
      matPtr->SetKz(currentKz);
    }
    
    if(isol == 0 && ncols == eigenvectors.Cols()){
      this->LoadSolution(GetEigenvectors());
    }else{
      const auto neqOriginal = eigenvectors.Rows();
      TPZFMatrix<CSTATE> evector(neqOriginal, ncols, 0.);
      eigenvectors.GetSub(0, isol, neqOriginal, ncols, evector);
      this->LoadSolution(evector);
    }
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshVecPost, m_cmesh_mf);
  }
  
  void WgmaAniso2D::WriteToCsv(std::string filename, STATE lambda){
    const auto &ev = this->GetEigenvalues();
    const int nev = ev.size();
    if(nev < 1){
      std::cout<<"There are no eigenvalues to write to .csv"<<std::endl;
      return;
    }
    std::cout<<"Exporting eigen info..."<<std::endl;
    std::ostringstream eigeninfo;
    constexpr auto pres = std::numeric_limits<STATE>::max_digits10;
    eigeninfo.precision(pres);
    REAL hSize = -1e12;
    TPZGeoMesh *gmesh = m_cmesh_mf->Reference();
    for(auto *el : gmesh->ElementVec()){
      if(el->HasSubElement()) continue;
      const auto elRadius = el->ElementRadius();
      hSize = elRadius > hSize ? elRadius : hSize;
    }
    const auto neq = m_cmesh_mf->NEquations();
    const auto nel = m_cmesh_mf->NElements();
    const auto porder = m_cmesh_mf->GetDefaultOrder();
    
    eigeninfo << std::fixed << neq << "," << nel << ",";
    eigeninfo << std::fixed << hSize << "," << porder<<",";
    eigeninfo << std::fixed << lambda<<",";
    eigeninfo << nev<<",";
    for(int i = 0; i < nev ; i++){
      eigeninfo<<std::fixed<<std::real(ev[i])<<",";
      eigeninfo<<std::fixed<<std::imag(ev[i]);
      if(i != nev - 1 ) {
        eigeninfo << ",";
      }
    }
    eigeninfo << std::endl;

    std::ofstream file(filename.c_str(),std::ios::app);
    file << eigeninfo.str();
    file.close();
  }

  WgmaPlanar::WgmaPlanar(TPZAutoPointer<TPZCompMesh> cmesh,
                                 const int n_threads, const bool reorder_eqs,
                                 const bool filter_bound) : m_cmesh(cmesh)
  {
    m_filter_bound = filter_bound;

    this->SetRenumber(new TPZCutHillMcKee());
    this->SetCompMeshInit(m_cmesh.operator->(),reorder_eqs);

    TPZAutoPointer<TPZStructMatrix> strmtrx{nullptr};
    if(using_tbb_mat){
      auto mtrx = new TPZSpStructMatrix<CSTATE,TPZStructMatrixOMPorTBB<CSTATE>>(m_cmesh);
      mtrx->SetShouldColor(false);
      mtrx->SetTBBorOMP(true);
      strmtrx=mtrx;
    }else{
      auto mtrx = new TPZSpStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>(m_cmesh);
      strmtrx=mtrx;
    }
    

    strmtrx->SetNumThreads(n_threads);
    auto ndofs = m_cmesh->NEquations();
    if(m_filter_bound){
      TPZVec<int64_t> activeEquations;
      const auto ndofs_before = m_cmesh->NEquations();
      std::set<int64_t> boundConnects;
    
      wgma::cmeshtools::FilterBoundaryEquations(m_cmesh, activeEquations,
                                                boundConnects);
      ndofs = activeEquations.size();
      std::cout<<"neq(before): "<<ndofs_before
               <<"\tneq(after): "<<ndofs<<std::endl;
      strmtrx->EquationFilter().SetActiveEquations(activeEquations);
    }else{
      std::cout<<"neq: "<<ndofs<<std::endl;
    }
    this->SetStructuralMatrix(strmtrx);
  }
  
  void WgmaPlanar::WriteToCsv(std::string filename, STATE lambda)
  {
    PZError<<__PRETTY_FUNCTION__
           <<"\nnot yet implemented...\n";
  }
  
  void WgmaPlanar::CountActiveEqs(int &neq)
  {
    auto eq_filt = this->StructMatrix()->EquationFilter();
    if(eq_filt.IsActive()){
      neq = eq_filt.NActiveEquations();
    }else{
      neq = m_cmesh->NEquations();
    }
  }

  
  void WgmaPlanar::LoadSolutionInternal(const int isol, const int nsol)
  {

    const auto &ev = this->GetEigenvalues();
    const auto &eigenvectors = this->GetEigenvectors();

    if(isol == 0 && nsol == eigenvectors.Cols()){
      this->LoadSolution(eigenvectors);
    }
    else{
      const auto neqOriginal = eigenvectors.Rows();
      TPZFMatrix<CSTATE> evector(neqOriginal, nsol, 0.);

      eigenvectors.GetSub(0, isol, neqOriginal, nsol, evector);
      this->LoadSolution(evector);
    }
  }


  void WgmaPeriodic2D::SetBeta(const CSTATE beta)
  {
    for(auto [id, mat] : m_cmesh->MaterialVec()){
        auto *mat_modal =
          dynamic_cast<TPZPeriodicWgma*>(mat);
        if(mat_modal){
          mat_modal->SetBeta(beta);
        }
      }
    m_beta = beta;
  }

  void WgmaPeriodic2D::AdjustSolver(TPZEigenSolver<CSTATE> *solv)
  {
    solv->SetTarget(m_beta*m_beta);
  }


  
  TPZAutoPointer<TPZCompMesh> CreateAtomicWgma2D(TPZAutoPointer<TPZGeoMesh> gmesh,
                                                 bool isH1,
                                                 int p,
                                                 const std::set<int> &volmats,
                                                 const std::set<int> &pmlmats,
                                                 
                                                 const std::vector<wgma::bc::data> &bcmats,
                                                 const std::vector<std::pair<int,int>> &probevec)
  {
    constexpr int dim = 2;
    constexpr bool isComplex{true};
  
    //for deRham compatibility
    const auto pOrder = isH1 ? p +1 : p;

    TPZAutoPointer<TPZCompMesh> cmesh =
      new TPZCompMesh(gmesh,isComplex);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    //number of state variables in the problem
    const int soldim = isH1 ? 1 : 3;

    std::set<int> allmats;

    TPZMaterialT<CSTATE> * dummyVolMat = nullptr;
    for(auto matid : volmats){
      auto *dummyMat =
        new wgma::materials::SolutionProjection<CSTATE>(matid,dim,soldim);
      cmesh->InsertMaterialObject(dummyMat);
      dummyVolMat = dummyMat;
      allmats.insert(matid);
    }
  
    for(auto id : pmlmats){
      auto *dummyMat =
        new wgma::materials::SolutionProjection<CSTATE>(id,dim,soldim);
        cmesh->InsertMaterialObject(dummyMat);
        allmats.insert(id);
    }

    for(auto [id,matdim] : probevec){
      auto *mat =
        new wgma::materials::SolutionProjection<CSTATE>(id,dim,soldim);
      cmesh->InsertMaterialObject(mat);
      allmats.insert(id);
    }

    TPZFNMatrix<1, CSTATE> val1(1, 1, 0);
    TPZManVector<CSTATE,1> val2(1, 0.);
    for(auto bc : bcmats){
      const int bctype = wgma::bc::to_int(bc.t);
      const int id = bc.id;
      const int volid = bc.volid;
      auto *dummyBC = dummyVolMat->CreateBC(dummyVolMat, id, bctype, val1, val2);
      cmesh->InsertMaterialObject(dummyBC);
      allmats.insert(id);
    }

    if(isH1){
      cmesh->SetAllCreateFunctionsContinuous();
    }else{
      cmesh->SetAllCreateFunctionsHCurl();
    }
    cmesh->AutoBuild(allmats);
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;
  }


  TPZAutoPointer<TPZCompMesh> CreateMfWgma2D(TPZAutoPointer<TPZGeoMesh> gmesh,
                                             cmeshtools::PhysicalData &data,
                                             const STATE lambda, const REAL &scale,
                                             const bool verbose)
  {
    constexpr int dim = 2;
    constexpr bool isComplex{true};
    
    auto &pmlDataVec = data.pmlvec;
    auto &bcDataVec = data.bcvec;
    
    const int nVolMats = data.matinfovec.size();
    const int nPmlMats = pmlDataVec.size();
    const int nBcMats = bcDataVec.size();
    
    TPZAutoPointer<TPZCompMesh> cmeshMF =
      new TPZCompMesh(gmesh,isComplex);
    cmeshMF->SetDimModel(dim);

    std::set<int> volmats;
    std::set<int> realvolmats;
    
    if(verbose){
      std::cout<<"inserting materials:\n";
    }
    for(auto [matid, er, ur] : data.matinfovec){
      auto *matWG = new TPZWgma(matid, er, ur, lambda, scale);
      cmeshMF->InsertMaterialObject(matWG);
      realvolmats.insert(matid);
      volmats.insert(matid);
      if(verbose){
        std::cout<<"\t id "<<matid<<" er "<<er<<" ur "<<ur<<'\n';
      }
    }

    if(verbose){
      std::cout<<"inserting pmls:\n";
    }
    //insert PML regions
    for(auto &pml : pmlDataVec){
      //skip PMLs of other dimensions
      if(pml->dim != cmeshMF->Dimension()){continue;}
      auto cart_pml = TPZAutoPointerDynamicCast<wgma::pml::cart::data>(pml);
      auto cyl_pml = TPZAutoPointerDynamicCast<wgma::pml::cyl::data>(pml);
      if(cart_pml){
        cart_pml->neigh =
          cmeshtools::AddRectangularPMLRegion<TPZWgma>(*cart_pml, realvolmats, gmesh, cmeshMF);
      }else if (cyl_pml){
        cyl_pml->neigh =
          cmeshtools::AddCylindricalPMLRegion<TPZWgma>(*cyl_pml, realvolmats, gmesh, cmeshMF);
      }else{
        DebugStop();
      }
      
      for(auto [id, _] : pml->neigh){
        volmats.insert(id);
      }
      if(verbose){
        std::cout<<"\tid (neighbour) ";
        for(auto [id, neigh] : pml->neigh){
          std::cout<<id<<" ("<<neigh<<") ";
        }
        if(cart_pml){
          std::cout<<"att dir "<<wgma::pml::cart::to_string(cart_pml->t)
                   <<" att coeff "<<cart_pml->alphax<<' '<<cart_pml->alphay<<std::endl;
        }else if(cyl_pml){
          std::cout<<"att dir "<<wgma::pml::cyl::to_string(cyl_pml->t)
                   <<" att coeff "<<cyl_pml->alphar<<std::endl;
        }
        else{
          DebugStop();
        }
      }
    }

    std::set<int> allmats = volmats;
    if(verbose){
      std::cout<<"inserting probes:\n";
    }
    for(auto [id,matdim] : data.probevec){
      static constexpr int nstate{1};
      auto *mat = new TPZNullMaterialCS<CSTATE>(id,matdim,nstate);
      cmeshMF->InsertMaterialObject(mat);
      allmats.insert(id);
      if(verbose){
        std::cout<<"\tid "<<id<<" dim "<<matdim<<std::endl;
      }
    }
    
    if(verbose){
      std::cout<<"inserting bcs:\n";
    }

    TPZFNMatrix<1, CSTATE> val1(1, 1, 0);
    TPZManVector<CSTATE,1> val2(1, 0.);
    
    for(auto bc : bcDataVec){
      const int bctype = wgma::bc::to_int(bc.t);
      const int id = bc.id;
      const int volid = bc.volid;
      auto *matWG =
        dynamic_cast<TPZMaterialT<CSTATE>*>(cmeshMF->FindMaterial(volid));
      if(!matWG){
        PZError<<__PRETTY_FUNCTION__
               <<"\n could not find material with id "<<id
               <<"\n.Is it a PML? Aborting..."<<std::endl;
        DebugStop();
      }
      auto *bcMat = matWG->CreateBC(matWG, id, bctype, val1, val2);
      cmeshMF->InsertMaterialObject(bcMat);
      allmats.insert(id);
      if(verbose){
        std::cout<<"\tid "<<id<<" vol mat "<<volid<<" type "<<wgma::bc::to_string(bc.t)<<std::endl;
      }
    }

    cmeshMF->SetDimModel(dim);
    cmeshMF->SetAllCreateFunctionsMultiphysicElem();

    cmeshMF->AutoBuild(allmats);
    cmeshMF->CleanUpUnconnectedNodes();

    return cmeshMF;
  }

  TPZAutoPointer<TPZCompMesh> CreateMfWgmaAniso2D(TPZAutoPointer<TPZGeoMesh> gmesh,
                                                  cmeshtools::PhysicalData &data,
                                                  const STATE lambda, const REAL &scale,
                                                  const bool verbose)
  {
    constexpr int dim = 2;
    constexpr bool isComplex{true};
    
    auto &pmlDataVec = data.pmlvec;
    auto &bcDataVec = data.bcvec;
    
    const int nVolMats = data.matinfovec.size();
    const int nPmlMats = pmlDataVec.size();
    const int nBcMats = bcDataVec.size();
    
    TPZAutoPointer<TPZCompMesh> cmeshMF =
      new TPZCompMesh(gmesh,isComplex);
    cmeshMF->SetDimModel(dim);

    std::set<int> volmats;
    std::set<int> realvolmats;
    
    if(verbose){
      std::cout<<"inserting materials:\n";
    }
    for(auto [matid, er, ur] : data.matinfovec){
      auto *matWG = new TPZAnisoWgma(matid, er, ur, lambda, scale);
      cmeshMF->InsertMaterialObject(matWG);
      realvolmats.insert(matid);
      volmats.insert(matid);
      if(verbose){
        std::cout<<"\t id "<<matid<<" er "<<er<<" ur "<<ur<<'\n';
      }
    }

    if(verbose){
      std::cout<<"inserting pmls:\n";
    }
    //insert PML regions
    for(auto &pml : pmlDataVec){
      //skip PMLs of other dimensions
      if(pml->dim != cmeshMF->Dimension()){continue;}

      auto cart_pml = TPZAutoPointerDynamicCast<wgma::pml::cart::data>(pml);
      auto cyl_pml = TPZAutoPointerDynamicCast<wgma::pml::cyl::data>(pml);
      if(cart_pml){
        cart_pml->neigh =
          cmeshtools::AddRectangularPMLRegion<TPZAnisoWgma>(*cart_pml, realvolmats, gmesh, cmeshMF);
      }else if (cyl_pml){
        cyl_pml->neigh =
          cmeshtools::AddCylindricalPMLRegion<TPZAnisoWgma>(*cyl_pml, realvolmats, gmesh, cmeshMF);
      }else{
        DebugStop();
      }
      for(auto id: pml->ids){
        volmats.insert(id);
      }
      if(verbose){
        std::cout<<"\tid (neighbour) ";
        for(auto [id, neigh] : pml->neigh){
          std::cout<<id<<" ("<<neigh<<") ";
        }
        if(cart_pml){
          std::cout<<"att dir "<<wgma::pml::cart::to_string(cart_pml->t)
                   <<" att coeff "<<cart_pml->alphax<<' '<<cart_pml->alphay<<std::endl;
        }else if(cyl_pml){
          std::cout<<"att dir "<<wgma::pml::cyl::to_string(cyl_pml->t)
                   <<" att coeff "<<cyl_pml->alphar<<std::endl;
        }
        else{
          DebugStop();
        }
      }
    }

    std::set<int> allmats = volmats;
    if(verbose){
      std::cout<<"inserting probes:\n";
    }
    for(auto [id,matdim] : data.probevec){
      static constexpr int nstate{1};
      auto *mat = new TPZNullMaterialCS<CSTATE>(id,matdim,nstate);
      cmeshMF->InsertMaterialObject(mat);
      allmats.insert(id);
      if(verbose){
        std::cout<<"\tid "<<id<<" dim "<<matdim<<std::endl;
      }
    }
    
    if(verbose){
      std::cout<<"inserting bcs:\n";
    }

    TPZFNMatrix<1, CSTATE> val1(1, 1, 0);
    TPZManVector<CSTATE,1> val2(1, 0.);
    
    for(auto bc : bcDataVec){
      const int bctype = wgma::bc::to_int(bc.t);
      const int id = bc.id;
      const int volid = bc.volid;

      allmats.insert(id);
      auto *matWG =
        dynamic_cast<TPZMaterialT<CSTATE>*>(cmeshMF->FindMaterial(volid));
      if(!matWG){
        PZError<<__PRETTY_FUNCTION__
               <<"\n could not find material with id "<<id
               <<"\n.Is it a PML? Aborting..."<<std::endl;
        DebugStop();
      }
      auto *bcMat = matWG->CreateBC(matWG, id, bctype, val1, val2);
      cmeshMF->InsertMaterialObject(bcMat);

      if(verbose){
        std::cout<<"\tid "<<id<<" vol mat "<<volid<<" type "<<wgma::bc::to_string(bc.t)<<std::endl;
      }
    }

    cmeshMF->SetDimModel(dim);
    cmeshMF->SetAllCreateFunctionsMultiphysicElem();

    cmeshMF->AutoBuild(allmats);
    cmeshMF->CleanUpUnconnectedNodes();

    return cmeshMF;
  }


  void SetupModalAnalysisMaterials(TPZAutoPointer<TPZGeoMesh> &gmesh, cmeshtools::PhysicalData& data,
                                   std::set<int>&volmats, std::set<int>&pmlmats)
  {
    constexpr int dim{2};

    // let us setup data for atomic meshes
    for(auto [matid, _, __] : data.matinfovec){
      volmats.insert(matid);
    }
    for(auto &pml : data.pmlvec){
      //skip PMLs of other dimensions
      if(pml->dim != dim){continue;}
      for(auto id : pml->ids){
        pmlmats.insert(id);
      }
    }
    std::set<int> allmats;

    std::set_union(volmats.begin(), volmats.end(),
                   pmlmats.begin(), pmlmats.end(),
                   std::inserter(allmats, allmats.begin()));
    /**let us associate each boundary with a given material.
       this is important for any non-homogeneous BCs*/
    for(auto &bc : data.bcvec){
      auto res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, bc.id, allmats);
      if(!res.has_value()){
        std::cout<<__PRETTY_FUNCTION__
                 <<"\nwarning: could not find neighbour of bc "<<bc.id<<std::endl;
      }
      bc.volid = res.value();
    }
  }
  
  TPZVec<TPZAutoPointer<TPZCompMesh>>
  CMeshWgma2D(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
              cmeshtools::PhysicalData &data,
              const STATE lambda, const REAL &scale,
              bool verbose)
  {
    TPZSimpleTimer timer ("Create cmesh");

    std::set<int> volmats, pmlmats;

    SetupModalAnalysisMaterials(gmesh,data,volmats,pmlmats);
    /*
      First we create the computational mesh associated with the H1 space
      (ez component)
    */
    bool ish1 = true;
    TPZAutoPointer<TPZCompMesh> cmeshH1 =
      CreateAtomicWgma2D(gmesh, ish1,pOrder,volmats,pmlmats, data.bcvec,data.probevec);

    /*
      Then we create the computational mesh associated with the HCurl space
    */
    ish1 = false;
    TPZAutoPointer<TPZCompMesh> cmeshHCurl =
      CreateAtomicWgma2D(gmesh, ish1,pOrder,volmats,pmlmats, data.bcvec,data.probevec);

    /*
      Now we create the MF mesh
    */
    TPZAutoPointer<TPZCompMesh> cmeshMF =
      CreateMfWgma2D(gmesh, data, lambda, scale, verbose);

    TPZManVector<TPZCompMesh*,3> meshVecIn(2);
    meshVecIn[TPZWgma::H1Index()] = cmeshH1.operator->();
    meshVecIn[TPZWgma::HCurlIndex()] = cmeshHCurl.operator->();

  
    TPZBuildMultiphysicsMesh::AddElements(meshVecIn, cmeshMF.operator->());
    TPZBuildMultiphysicsMesh::AddConnects(meshVecIn, cmeshMF.operator->());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshVecIn, cmeshMF.operator->());

    cmeshMF->ExpandSolution();
    cmeshMF->ComputeNodElCon();
    cmeshMF->CleanUpUnconnectedNodes();

    TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec(3,nullptr);
    meshVec[0] = cmeshMF;
    meshVec[1 + TPZWgma::H1Index()] = cmeshH1;
    meshVec[1 + TPZWgma::HCurlIndex()] = cmeshHCurl;
    return meshVec;
  
  }

  TPZVec<TPZAutoPointer<TPZCompMesh>>
  CMeshWgmaAniso2D(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
                   cmeshtools::PhysicalData &data,
                   const STATE lambda, const REAL &scale,
                   bool verbose)
  {
    TPZSimpleTimer timer ("Create cmesh");

    constexpr int dim{2};

    // let us setup data for atomic meshes
    std::set<int> volmats;
    std::set<int> pmlmats;

    SetupModalAnalysisMaterials(gmesh, data, volmats, pmlmats);
    /*
      First we create the computational mesh associated with the H1 space
      (ez component)
    */
    bool ish1 = true;
    TPZAutoPointer<TPZCompMesh> cmeshH1 =
      CreateAtomicWgma2D(gmesh, ish1,pOrder,volmats,pmlmats, data.bcvec,data.probevec);

    /*
      Then we create the computational mesh associated with the HCurl space
    */
    ish1 = false;
    TPZAutoPointer<TPZCompMesh> cmeshHCurl =
      CreateAtomicWgma2D(gmesh, ish1,pOrder,volmats,pmlmats, data.bcvec,data.probevec);

    /*
      Now we create the MF mesh
    */
    TPZAutoPointer<TPZCompMesh> cmeshMF =
      CreateMfWgmaAniso2D(gmesh, data, lambda, scale, verbose);

    TPZManVector<TPZCompMesh*,3> meshVecIn(2);
    meshVecIn[TPZAnisoWgma::H1Index()] = cmeshH1.operator->();
    meshVecIn[TPZAnisoWgma::HCurlIndex()] = cmeshHCurl.operator->();

  
    TPZBuildMultiphysicsMesh::AddElements(meshVecIn, cmeshMF.operator->());
    TPZBuildMultiphysicsMesh::AddConnects(meshVecIn, cmeshMF.operator->());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshVecIn, cmeshMF.operator->());

    cmeshMF->ExpandSolution();
    cmeshMF->ComputeNodElCon();
    cmeshMF->CleanUpUnconnectedNodes();

    TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec(3,nullptr);
    meshVec[0] = cmeshMF;
    meshVec[1 + TPZAnisoWgma::H1Index()] = cmeshH1;
    meshVec[1 + TPZAnisoWgma::HCurlIndex()] = cmeshHCurl;
    return meshVec;
  
  }
  

  TPZVec<TPZAutoPointer<TPZCompMesh>>
  CMeshWgma2DPeriodic(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
                      cmeshtools::PhysicalData &data,
                      const std::map<int64_t,int64_t> &periodic_els,
                      const STATE lambda, const REAL &scale,
                      bool verbose)
  {
    TPZSimpleTimer timer ("Create cmesh");

    std::set<int> volmats,pmlmats;

    SetupModalAnalysisMaterials(gmesh, data, volmats, pmlmats);

    
    /*
      First we create the computational mesh associated with the H1 space
      (ez component)
    */
    bool ish1 = true;
    TPZAutoPointer<TPZCompMesh> cmeshH1 =
      CreateAtomicWgma2D(gmesh, ish1,pOrder,volmats,pmlmats, data.bcvec,data.probevec);

    /*
      Now we add the periodicity
     */

    wgma::cmeshtools::SetPeriodic(cmeshH1,periodic_els);
    /*
      Then we create the computational mesh associated with the HCurl space
    */
    ish1 = false;
    TPZAutoPointer<TPZCompMesh> cmeshHCurl =
      CreateAtomicWgma2D(gmesh, ish1,pOrder,volmats,pmlmats, data.bcvec,data.probevec);
    /*
      Now we add the periodicity
     */
    wgma::cmeshtools::SetPeriodic(cmeshHCurl,periodic_els);
    /*
      Now we create the MF mesh
    */
    TPZAutoPointer<TPZCompMesh> cmeshMF =
      CreateMfWgma2D(gmesh, data, lambda, scale, verbose);

    TPZManVector<TPZCompMesh*,3> meshVecIn(2);
    meshVecIn[TPZWgma::H1Index()] = cmeshH1.operator->();
    meshVecIn[TPZWgma::HCurlIndex()] = cmeshHCurl.operator->();

  
    TPZBuildMultiphysicsMesh::AddElements(meshVecIn, cmeshMF.operator->());
    TPZBuildMultiphysicsMesh::AddConnects(meshVecIn, cmeshMF.operator->());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshVecIn, cmeshMF.operator->());

    cmeshMF->ExpandSolution();
    cmeshMF->ComputeNodElCon();
    cmeshMF->CleanUpUnconnectedNodes();

    TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec(3,nullptr);
    meshVec[0] = cmeshMF;
    meshVec[1 + TPZWgma::H1Index()] = cmeshH1;
    meshVec[1 + TPZWgma::HCurlIndex()] = cmeshHCurl;
    return meshVec;
  
  }
  
  TPZAutoPointer<TPZCompMesh>
  CMeshWgma1D(TPZAutoPointer<TPZGeoMesh> gmesh,
              wgma::planarwg::mode mode, int pOrder,
              wgma::cmeshtools::PhysicalData &data,
              const STATE lambda, const REAL scale)
  {
    return CMeshWgma1DPeriodic(gmesh,mode,pOrder,data,
                               {},lambda,scale);
  }

  TPZAutoPointer<TPZCompMesh>
  CMeshWgma1DPeriodic(TPZAutoPointer<TPZGeoMesh> gmesh,
                      wgma::planarwg::mode mode, int pOrder,
                      wgma::cmeshtools::PhysicalData &data,
                      const std::map<int64_t,int64_t> periodic_els,
                      const STATE lambda, const REAL scale)
  {
  
    static constexpr bool isComplex{true};
    static constexpr int dim{1};
    TPZAutoPointer<TPZCompMesh> cmeshH1 = new TPZCompMesh(gmesh, isComplex);
    cmeshH1->SetDimModel(dim);

    const int nvolmats = data.matinfovec.size();

    // insert volumetric mats
    std::set<int> volmats;
    std::set<int> allmats;
    TPZPeriodicWgma::ModeType matmode;
    switch (mode) {
    case wgma::planarwg::mode::TE:
      matmode = TPZScalarField::ModeType::TE;
      break;
    case wgma::planarwg::mode::TM:
      matmode = TPZScalarField::ModeType::TM;
      break;
    case wgma::planarwg::mode::Invalid:
      DebugStop();
      break;
    }
    for (auto [id, er, ur] : data.matinfovec) {
      auto *mat =
        new TPZPlanarWgma(id, er, ur, lambda, matmode, scale);
      cmeshH1->InsertMaterialObject(mat);
      // for pml
      volmats.insert(id);
      //for assembling only desired materials
      allmats.insert(id);
    }

    for (auto &pml : data.pmlvec) {
      //skip PMLs of other dimensions
      if(pml->dim != cmeshH1->Dimension()){continue;}
      auto cart_pml = TPZAutoPointerDynamicCast<wgma::pml::cart::data>(pml);
      auto cyl_pml = TPZAutoPointerDynamicCast<wgma::pml::cyl::data>(pml);
      if(cart_pml){
        cart_pml->neigh =
          cmeshtools::AddRectangularPMLRegion<TPZPlanarWgma>(*cart_pml, volmats, gmesh, cmeshH1);
      }else if (cyl_pml){
        cyl_pml->neigh =
          cmeshtools::AddCylindricalPMLRegion<TPZPlanarWgma>(*cyl_pml, volmats, gmesh, cmeshH1);
      }else{
        DebugStop();
      }
      
      for( auto id : pml->ids){
        allmats.insert(id);
      }
    }
    

    for(auto [id,matdim] : data.probevec){
      static constexpr int soldim{1};
      auto *mat = new wgma::materials::SolutionProjection<CSTATE>(id,matdim,soldim);
      cmeshH1->InsertMaterialObject(mat);
      allmats.insert(id);
    }
    TPZFNMatrix<1, CSTATE> val1(1, 1, 0);
    TPZManVector<CSTATE, 1> val2(1, 0.);

    /**let us associate each boundary with a given material.
       this is important for the source boundary*/
    for(auto &bc : data.bcvec){
      auto res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, bc.id, allmats);
      if(!res.has_value()){
        std::cout<<__PRETTY_FUNCTION__
                 <<"\nError: could not find neighbour of bc "<<bc.id<<std::endl;
        DebugStop();
      }else{
        bc.volid = res.value();
      }
    }

    // for(auto bc : data.bcvec){
    //   std::cout<<"bc "<<bc.id<<" mat "<<bc.volid<<std::endl;
    // }

  
    for (auto bc : data.bcvec) {
      const int bctype = wgma::bc::to_int(bc.t);
      const int id = bc.id;
      const int volmatid = bc.volid;
      auto *volmat =
        dynamic_cast<TPZMaterialT<CSTATE> *>(cmeshH1->FindMaterial(volmatid));
      auto *bcmat = volmat->CreateBC(volmat, id, bctype, val1, val2);
      cmeshH1->InsertMaterialObject(bcmat);
      allmats.insert(id);
    }

  
  

    cmeshH1->SetAllCreateFunctionsContinuous();
    cmeshH1->SetDefaultOrder(pOrder);
    cmeshH1->AutoBuild(allmats);

    if(!periodic_els.size()) return cmeshH1;

    wgma::cmeshtools::SetPeriodic(cmeshH1,periodic_els);
    return cmeshH1;
  }

  
  TPZAutoPointer<TPZCompMesh>
  CMeshWgmaPeriodic(TPZAutoPointer<TPZGeoMesh> gmesh,
                    const wgma::planarwg::mode mode, int pOrder,
                    wgma::cmeshtools::PhysicalData &data,
                    const std::map<int64_t,int64_t> periodic_els,
                    const STATE lambda,
                    const REAL scale)
  {
  
    static constexpr bool isComplex{true};
    static constexpr int dim{2};
    TPZAutoPointer<TPZCompMesh> cmeshH1 = new TPZCompMesh(gmesh, isComplex);
    cmeshH1->SetDimModel(dim);

    const int nvolmats = data.matinfovec.size();

    // insert volumetric mats
    std::set<int> volmats;
    std::set<int> allmats;
    TPZPeriodicWgma::ModeType matmode;
    switch (mode) {
    case wgma::planarwg::mode::TE:
      matmode = TPZPeriodicWgma::ModeType::TE;
      break;
    case wgma::planarwg::mode::TM:
      matmode = TPZPeriodicWgma::ModeType::TM;
      break;
    case wgma::planarwg::mode::Invalid:
      DebugStop();
      break;
    }
    for (auto [id, er, ur] : data.matinfovec) {
      auto *mat =
        new TPZPeriodicWgma(id, er, ur, lambda, matmode, scale);
      cmeshH1->InsertMaterialObject(mat);
      // for pml
      volmats.insert(id);
      //for assembling only desired materials
      allmats.insert(id);
    }

    for (auto &pml : data.pmlvec) {
      //skip PMLs of other dimensions
      if(pml->dim != cmeshH1->Dimension()){continue;}
      auto cart_pml = TPZAutoPointerDynamicCast<wgma::pml::cart::data>(pml);
      auto cyl_pml = TPZAutoPointerDynamicCast<wgma::pml::cyl::data>(pml);
      if(cart_pml){
        cart_pml->neigh =
          cmeshtools::AddRectangularPMLRegion<TPZWgma>(*cart_pml, volmats, gmesh, cmeshH1);
      }else if (cyl_pml){
        cyl_pml->neigh =
          cmeshtools::AddCylindricalPMLRegion<TPZWgma>(*cyl_pml, volmats, gmesh, cmeshH1);
      }else{
        DebugStop();
      }
      
      for( auto id : pml->ids){
        allmats.insert(id);
      }
    }

    for(auto [id,matdim] : data.probevec){
      static constexpr int soldim{1};
      auto *mat =
        new wgma::materials::SolutionProjection<CSTATE>(id,matdim,soldim);
      cmeshH1->InsertMaterialObject(mat);
      allmats.insert(id);
    }
    
    TPZFNMatrix<1, CSTATE> val1(1, 1, 0);
    TPZManVector<CSTATE, 1> val2(1, 0.);

    /**let us associate each boundary with a given material.
       this is important for the source boundary*/
    for(auto &bc : data.bcvec){
      auto res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, bc.id, allmats);
      if(!res.has_value()){
        std::cout<<__PRETTY_FUNCTION__
                 <<"\nError: could not find neighbour of bc "<<bc.id<<std::endl;
        DebugStop();
      }else{
        bc.volid = res.value();
      }
    }

    // for(auto bc : data.bcvec){
    //   std::cout<<"bc "<<bc.id<<" mat "<<bc.volid<<std::endl;
    // }

  
    for (auto bc : data.bcvec) {
      const int bctype = wgma::bc::to_int(bc.t);
      const int id = bc.id;
      const int volmatid = bc.volid;
      auto *volmat =
        dynamic_cast<TPZMaterialT<CSTATE> *>(cmeshH1->FindMaterial(volmatid));
      auto *bcmat = volmat->CreateBC(volmat, id, bctype, val1, val2);
      cmeshH1->InsertMaterialObject(bcmat);
      allmats.insert(id);
    }

  
  

    cmeshH1->SetAllCreateFunctionsContinuous();
    cmeshH1->SetDefaultOrder(pOrder);
    cmeshH1->AutoBuild(allmats);

    wgma::cmeshtools::SetPeriodic(cmeshH1, periodic_els);
    
    return cmeshH1;
  }
  
};
