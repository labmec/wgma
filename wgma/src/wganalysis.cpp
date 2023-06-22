#include "wganalysis.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "cmeshtools_impl.hpp"

#include <TPZSpStructMatrix.h>
#include <TPZKrylovEigenSolver.h>
#include <TPZQuadEigenSolver.h>
#include <Electromagnetics/TPZWgma.h>
#include <Electromagnetics/TPZAnisoWgma.h>
#include <Electromagnetics/TPZPeriodicWgma.h>
#include <Electromagnetics/TPZPlanarWgma.h>
#include <TPZNullMaterial.h>
#include <TPZNullMaterialCS.h>
#include <TPZSimpleTimer.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzelctemp.h>

#include <cassert>

namespace wgma::wganalysis{

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
  
  void Wgma::Solve(bool compute_eigenvectors){
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

    std::ios cout_state(nullptr);
    cout_state.copyfmt(std::cout);
    
    std::cout << std::setprecision(std::numeric_limits<STATE>::max_digits10);

    for(auto &w : this->GetEigenvalues()){
      std::cout<<w<<std::endl;
    }
    std::cout.copyfmt(cout_state);
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

    this->SetCompMeshInit(m_cmesh_mf.operator->(), reorder_eqs);

    TPZAutoPointer<TPZStructMatrix> strmtrx =
      new TPZSpStructMatrix<CSTATE>(m_cmesh_mf);

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

    const auto neqOriginal = eigenvectors.Rows();
    TPZFMatrix<CSTATE> evector(neqOriginal, ncols, 0.);
    
    TPZManVector<TPZAutoPointer<TPZCompMesh>,2> meshVecPost(2);
    meshVecPost[0] = m_cmesh_h1;
    meshVecPost[1] = m_cmesh_hcurl;
  
    const CSTATE currentKz = [&ev,isol](){
      auto tmp = std::sqrt(-1.0*ev[isol]);
      constexpr auto epsilon = std::numeric_limits<STATE>::epsilon()/
        (10*std::numeric_limits<STATE>::digits10);
      //let us discard extremely small imag parts
      if (tmp.imag() < epsilon)
        {tmp = tmp.real();}
      return tmp;
    }();

    

    eigenvectors.GetSub(0, isol, neqOriginal, ncols, evector);
    for(auto mat : m_cmesh_mf->MaterialVec()){
      auto id = mat.first;
      auto matPtr =
        dynamic_cast<TPZWgma *>(m_cmesh_mf->FindMaterial(id));
      if(!matPtr) continue;
      matPtr->SetKz(currentKz);
    }
    this->LoadSolution(evector);
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

    this->SetCompMeshInit(m_cmesh_mf.operator->(), reorder_eqs);

    TPZAutoPointer<TPZStructMatrix> strmtrx =
      new TPZSpStructMatrix<CSTATE>(m_cmesh_mf);

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

    const auto neqOriginal = eigenvectors.Rows();
    TPZFMatrix<CSTATE> evector(neqOriginal, ncols, 0.);
    
    TPZManVector<TPZAutoPointer<TPZCompMesh>,2> meshVecPost(2);
    meshVecPost[0] = m_cmesh_h1;
    meshVecPost[1] = m_cmesh_hcurl;
  
    const CSTATE currentKz = [&ev,isol](){
      auto tmp = CSTATE(-1.0*1i)*ev[isol];
      constexpr auto epsilon = std::numeric_limits<STATE>::epsilon()/
        (10*std::numeric_limits<STATE>::digits10);
      //let us discard extremely small imag parts
      if (tmp.imag() < epsilon)
        {tmp = tmp.real();}
      return tmp;
    }();

    

    eigenvectors.GetSub(0, isol, neqOriginal, ncols, evector);
    for(auto mat : m_cmesh_mf->MaterialVec()){
      auto id = mat.first;
      auto matPtr =
        dynamic_cast<TPZWgma *>(m_cmesh_mf->FindMaterial(id));
      if(!matPtr) continue;
      matPtr->SetKz(currentKz);
    }
    this->LoadSolution(evector);
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
    
    this->SetCompMeshInit(m_cmesh.operator->(),reorder_eqs);

    TPZAutoPointer<TPZStructMatrix> strmtrx =
      new TPZSpStructMatrix<CSTATE>(m_cmesh);

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
    
    const auto neqOriginal = eigenvectors.Rows();
    TPZFMatrix<CSTATE> evector(neqOriginal, nsol, 0.);

    eigenvectors.GetSub(0, isol, neqOriginal, nsol, evector);
    this->LoadSolution(evector);
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
    constexpr int nState = 1;


    TPZMaterialT<CSTATE> * dummyVolMat = nullptr;
    for(auto matid : volmats){
      auto *dummyMat = new TPZNullMaterial<CSTATE>(matid,dim,nState);
      cmesh->InsertMaterialObject(dummyMat);
      dummyVolMat = dummyMat;
    }
  
    for(auto id : pmlmats){
      auto *dummyMat = new TPZNullMaterial<CSTATE>(id,dim,nState);
        cmesh->InsertMaterialObject(dummyMat);
    }

    for(auto [id,matdim] : probevec){
      static constexpr int nstate{1};
      auto *mat = new TPZNullMaterial<CSTATE>(id,dim,nstate);
      cmesh->InsertMaterialObject(mat);
    }

    TPZFNMatrix<1, CSTATE> val1(1, 1, 0);
    TPZManVector<CSTATE,1> val2(1, 0.);
    for(auto bc : bcmats){
      const int bctype = wgma::bc::to_int(bc.t);
      const int id = bc.id;
      const int volid = bc.volid;
      auto *dummyBC = dummyVolMat->CreateBC(dummyVolMat, id, bctype, val1, val2);
      cmesh->InsertMaterialObject(dummyBC);
    }

    if(isH1){
      cmesh->SetAllCreateFunctionsContinuous();
    }else{
      cmesh->SetAllCreateFunctionsHCurl();
    }
    cmesh->AutoBuild();
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

    if(verbose){
      std::cout<<"inserting probes:\n";
    }
    for(auto [id,matdim] : data.probevec){
      static constexpr int nstate{1};
      auto *mat = new TPZNullMaterialCS<CSTATE>(id,matdim,nstate);
      cmeshMF->InsertMaterialObject(mat);
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

      if(verbose){
        std::cout<<"\tid "<<id<<" vol mat "<<volid<<" type "<<wgma::bc::to_string(bc.t)<<std::endl;
      }
    }

    cmeshMF->SetDimModel(dim);
    cmeshMF->SetAllCreateFunctionsMultiphysicElem();

    cmeshMF->AutoBuild();
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

    if(verbose){
      std::cout<<"inserting probes:\n";
    }
    for(auto [id,matdim] : data.probevec){
      static constexpr int nstate{1};
      auto *mat = new TPZNullMaterialCS<CSTATE>(id,matdim,nstate);
      cmeshMF->InsertMaterialObject(mat);
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

      if(verbose){
        std::cout<<"\tid "<<id<<" vol mat "<<volid<<" type "<<wgma::bc::to_string(bc.t)<<std::endl;
      }
    }

    cmeshMF->SetDimModel(dim);
    cmeshMF->SetAllCreateFunctionsMultiphysicElem();

    cmeshMF->AutoBuild();
    cmeshMF->CleanUpUnconnectedNodes();

    return cmeshMF;
  }
  
  TPZVec<TPZAutoPointer<TPZCompMesh>>
  CMeshWgma2D(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
              cmeshtools::PhysicalData &data,
              const STATE lambda, const REAL &scale,
              bool verbose)
  {
    TPZSimpleTimer timer ("Create cmesh");

    constexpr int dim{2};

    // let us setup data for atomic meshes
    std::set<int> volmats;
    for(auto [matid, _, __] : data.matinfovec){
      volmats.insert(matid);
    }
    std::set<int> pmlmats;
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
    for(auto [matid, _, __] : data.matinfovec){
      volmats.insert(matid);
    }
    std::set<int> pmlmats;
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
                      std::map<int64_t,int64_t> &periodic_els,
                      const STATE lambda, const REAL &scale,
                      bool verbose)
  {
    TPZSimpleTimer timer ("Create cmesh");

    constexpr int dim{2};

    // let us setup data for atomic meshes
    std::set<int> volmats;
    for(auto [matid, _, __] : data.matinfovec){
      volmats.insert(matid);
    }
    std::set<int> pmlmats;
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


    auto SetPeriodic = [&periodic_els](TPZAutoPointer<TPZCompMesh> cmesh){
      auto gmesh = cmesh->Reference();
      gmesh->ResetReference();
      cmesh->LoadReferences();
      //let us copy the connects
      for(auto [dep, indep] : periodic_els){
        //geometric elements
        auto *dep_gel = gmesh->Element(dep);
        const auto *indep_gel = gmesh->Element(indep);
        //computational element
        auto *indep_cel = indep_gel->Reference();
        auto *dep_cel = dep_gel->Reference();
        //number of connects
        const auto n_dep_con = dep_cel->NConnects();
        const auto n_indep_con = indep_cel->NConnects();
        //just to be sure
        assert(n_dep_con == n_indep_con);

        //now we create dependencies between connects
        for(auto ic = 0; ic < n_indep_con; ic++){
          const auto indep_ci = indep_cel->ConnectIndex(ic);
          const auto dep_ci = dep_cel->ConnectIndex(ic);

          auto &dep_con = dep_cel->Connect(ic);
          const auto ndof = dep_con.NDof(cmesh);
          if(ndof==0) {continue;}
          constexpr int64_t ipos{0};
          constexpr int64_t jpos{0};
      
          TPZFMatrix<REAL> mat(ndof,ndof);
          mat.Identity();
          dep_con.AddDependency(dep_ci, indep_ci, mat, ipos,jpos,ndof,ndof);
        } 
      }
      cmesh->CleanUpUnconnectedNodes();
      cmesh->ExpandSolution();
    };
    
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

    SetPeriodic(cmeshH1);
    /*
      Then we create the computational mesh associated with the HCurl space
    */
    ish1 = false;
    TPZAutoPointer<TPZCompMesh> cmeshHCurl =
      CreateAtomicWgma2D(gmesh, ish1,pOrder,volmats,pmlmats, data.bcvec,data.probevec);
    /*
      Now we add the periodicity
     */
    SetPeriodic(cmeshHCurl);
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
      static constexpr int nstate{1};
      auto *mat = new TPZNullMaterial<CSTATE>(id,matdim,nstate);
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
    return cmeshH1;
  }
  
  TPZAutoPointer<TPZCompMesh>
  CMeshWgmaPeriodic(TPZAutoPointer<TPZGeoMesh> gmesh,
                    const wgma::planarwg::mode mode, int pOrder,
                    wgma::cmeshtools::PhysicalData &data,
                    std::map<int64_t,int64_t> periodic_els, const STATE lambda,
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
      static constexpr int nstate{1};
      auto *mat = new TPZNullMaterial<CSTATE>(id,matdim,nstate);
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

  
    //let us copy the connects
    for(auto [dep, indep] : periodic_els){
      //geometric elements
      auto *dep_gel = gmesh->Element(dep);
      const auto *indep_gel = gmesh->Element(indep);
      //computational element
      //we need interpolatedelement for SideConnectLocId
      auto *indep_cel =
        dynamic_cast<TPZInterpolatedElement*>(indep_gel->Reference());
      auto *dep_cel =
        dynamic_cast<TPZInterpolatedElement*>(dep_gel->Reference());
      //number of connects
      const auto n_dep_con = dep_cel->NConnects();
      const auto n_indep_con = indep_cel->NConnects();
      //just to be sure
      assert(n_dep_con == n_indep_con);

      //now we create dependencies between connects
      for(auto ic = 0; ic < n_indep_con; ic++){
        const auto indep_ci = indep_cel->ConnectIndex(ic);
        const auto dep_ci = dep_cel->ConnectIndex(ic);

        auto &dep_con = dep_cel->Connect(ic);
        const auto ndof = dep_con.NDof(cmeshH1);
        if(ndof==0) {continue;}
        constexpr int64_t ipos{0};
        constexpr int64_t jpos{0};
      
        TPZFMatrix<REAL> mat(ndof,ndof);
        mat.Identity();
        dep_con.AddDependency(dep_ci, indep_ci, mat, ipos,jpos,ndof,ndof);
      } 
    }
    cmeshH1->CleanUpUnconnectedNodes();
    cmeshH1->ExpandSolution();
    return cmeshH1;
  }
  
};