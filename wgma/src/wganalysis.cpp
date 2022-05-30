#include "wganalysis.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "cmeshtools_impl.hpp"

#include <TPZSpStructMatrix.h>
#include <TPZKrylovEigenSolver.h>
#include <Electromagnetics/TPZWgma.h>
#include <Electromagnetics/TPZPeriodicWgma.h>
#include <TPZNullMaterial.h>
#include <TPZSimpleTimer.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzelctemp.h>

#include <cassert>

namespace wgma::wganalysis{

  Wgma::~Wgma(){}
  
  void Wgma::SetSolver(const TPZEigenSolver<CSTATE> & solv){
    m_an->SetSolver(solv);
  }
    
  TPZEigenSolver<CSTATE> & Wgma::GetSolver() const{
    return m_an->EigenSolver<CSTATE>();
  }

  void Wgma::Run(bool compute_eigenvectors){
    m_an->SetComputeEigenvectors(compute_eigenvectors);
    TPZEigenSolver<CSTATE> *solv =
      dynamic_cast<TPZEigenSolver<CSTATE>*>(m_an->Solver());
    if(!solv){
      std::cerr<<__PRETTY_FUNCTION__
               <<"\nA solver has not been set.\n"
               <<"Check documentation of TPZKrylovEigenSolver"
               <<"\nor wgma::slepc::EPSHandler\nAborting..."<<std::endl;
      exit(-1);
    }

    AdjustSolver(solv);
    
    {//scope for timer
      std::cout<<"Assembling..."<<std::flush;
      TPZSimpleTimer assemble("Assemble");
      m_an->Assemble();
      std::cout<<"\rAssembled!"<<std::endl;
    }
    {//scope for timer
      TPZSimpleTimer solv("Solve");
      std::cout<<"Solving..."<<std::flush;
      m_an->Solve();
      std::cout<<"\rSolved!"<<std::endl;
    }

    m_evalues = m_an->GetEigenvalues();

    std::ios cout_state(nullptr);
    cout_state.copyfmt(std::cout);
    
    std::cout << std::setprecision(std::numeric_limits<STATE>::max_digits10);

    for(auto &w : m_evalues){
      std::cout<<w<<std::endl;
    }
    std::cout.copyfmt(cout_state);
    
    if(m_an->ComputeEigenvectors()){
      m_evectors = m_an->GetEigenvectors();
    }
  }


  TPZVec<CSTATE> Wgma::GetEigenvalues() const
  {
    return m_evalues;
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
    m_cmesh_h1 = meshvec[TPZWgma::H1Index()];
    m_cmesh_hcurl = meshvec[TPZWgma::HCurlIndex()];
  
    m_an = new TPZEigenAnalysis(m_cmesh_mf, reorder_eqs);

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
    }
    CountActiveEqs(m_n_dofs_mf,m_n_dofs_h1,m_n_dofs_hcurl);
    m_an->SetStructuralMatrix(strmtrx);
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
  void Wgma2D::PostProcess(std::string filename,
                               const int vtk_res,
                               const bool print_real_part){
    if(!m_an->ComputeEigenvectors()){
      std::cout<<__PRETTY_FUNCTION__
               <<"\nOnly eigenvalues were calculated.\n"
               <<"Nothing to do here...\n";
      return;
    }

    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Ez_real");
    scalnames.Push("Ez_abs");
    vecnames.Push("Et_real");
    vecnames.Push("Et_abs");
    const std::string plotfile = filename+".vtk";
    constexpr int dim{2};
    m_an->DefineGraphMesh(dim, scalnames, vecnames,plotfile);

    const auto neqOriginal = m_evectors.Rows();
    TPZFMatrix<CSTATE> evector(neqOriginal, 1, 0.);
    
    auto &ev = m_evalues;
    auto &eigenvectors = m_evectors;
    
    const auto nev = ev.size();

    TPZManVector<TPZAutoPointer<TPZCompMesh>,2> meshVecPost(2);
    meshVecPost[0] = m_cmesh_h1;
    meshVecPost[1] = m_cmesh_hcurl;
  
    std::cout<<"Post processing..."<<std::endl;
  
    for (int iSol = 0; iSol < ev.size(); iSol++) {
      const CSTATE currentKz = [&ev,iSol](){
        auto tmp = std::sqrt(-1.0*ev[iSol]);
        constexpr auto epsilon = std::numeric_limits<STATE>::epsilon()/
          (10*std::numeric_limits<STATE>::digits10);
        //let us discard extremely small imag parts
        if (tmp.imag() < epsilon)
          {tmp = tmp.real();}
        return tmp;
      }();
      eigenvectors.GetSub(0, iSol, neqOriginal, 1, evector);
      for(auto mat : m_cmesh_mf->MaterialVec()){
        auto id = mat.first;
        auto matPtr =
          dynamic_cast<TPZWgma *>(m_cmesh_mf->FindMaterial(id));
        if(!matPtr) continue;
        matPtr->SetKz(currentKz);
      }
      std::cout<<"\rPost processing step "<<iSol+1<<" out of "<<ev.size()
               <<"(kz = "<<currentKz<<")"<<std::flush;
      m_an->LoadSolution(evector);
      TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshVecPost, m_cmesh_mf);
      m_an->PostProcess(vtk_res);
    }
    std::cout<<"\nFinished post processing"<<std::endl;
    std::cout<<std::endl;
  }
  
  void Wgma2D::WriteToCsv(std::string filename, STATE lambda){
    const int nev = m_evalues.size();
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
      eigeninfo<<std::fixed<<std::real(m_evalues[i])<<",";
      eigeninfo<<std::fixed<<std::imag(m_evalues[i]);
      if(i != nev - 1 ) {
        eigeninfo << ",";
      }
    }
    eigeninfo << std::endl;

    std::ofstream file(filename.c_str(),std::ios::app);
    file << eigeninfo.str();
    file.close();
  }


  WgmaPeriodic2D::WgmaPeriodic2D(TPZAutoPointer<TPZCompMesh> cmesh,
                                 const int n_threads, const bool reorder_eqs,
                                 const bool filter_bound) : m_cmesh(cmesh)
  {
    m_filter_bound = filter_bound;
    
    m_an = new TPZEigenAnalysis(m_cmesh, reorder_eqs);

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
    m_an->SetStructuralMatrix(strmtrx);
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

  void WgmaPeriodic2D::PostProcess(std::string filename, const int vtk_res)
  {
    if(!m_an->ComputeEigenvectors()){
      std::cout<<__PRETTY_FUNCTION__
               <<"\nOnly eigenvalues were calculated.\n"
               <<"Nothing to do here...\n";
      return;
    }
    ///vtk export
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Field_real");
    scalnames.Push("Field_imag");
    scalnames.Push("Field_abs");
    vecnames.Push("Deriv_real");
    vecnames.Push("Deriv_imag");
    vecnames.Push("Deriv_abs");
    const std::string plotfile = filename+".vtk";
    m_an->DefineGraphMesh(2,scalnames,vecnames,plotfile);
    const auto neqOriginal = m_evectors.Rows();
    TPZFMatrix<CSTATE> evector(neqOriginal, 1, 0.);

    m_evectors.GetSub(0, 0, neqOriginal, 1, evector);
    m_an->LoadSolution(evector);
    m_an->PostProcess(vtk_res);
  }
  
  void WgmaPeriodic2D::WriteToCsv(std::string filename, STATE lambda)
  {
    PZError<<__PRETTY_FUNCTION__
           <<"\nnot yet implemented...\n";
  }
  
  void WgmaPeriodic2D::CountActiveEqs(int &neq)
  {
    auto eq_filt = m_an->StructMatrix()->EquationFilter();
    if(eq_filt.IsActive()){
      neq = eq_filt.NActiveEquations();
    }else{
      neq = m_cmesh->NEquations();
    }
  }
  
  void WgmaPeriodic2D::AdjustSolver(TPZEigenSolver<CSTATE> *solv)
  {
    solv->SetTarget(m_beta*m_beta);
  }
  

  TPZVec<TPZAutoPointer<TPZCompMesh>>
  CMeshWgma2D(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
              cmeshtools::PhysicalData &data,
              const STATE lambda, const REAL &scale)
  {
    TPZSimpleTimer timer ("Create cmesh");
    constexpr int dim = 2;
    constexpr bool isComplex{true};
    auto &pmlDataVec = data.pmlvec;
    auto &bcDataVec = data.bcvec;
    const int nVolMats = data.matinfovec.size();
    const int nPmlMats = pmlDataVec.size();
    const int nBcMats = bcDataVec.size();
  
  
    /*
      First we create the computational mesh associated with the H1 space
      (ez component)*/
    TPZAutoPointer<TPZCompMesh> cmeshH1 =
      new TPZCompMesh(gmesh,isComplex);
    cmeshH1->SetDefaultOrder(pOrder +1);//for deRham compatibility
    cmeshH1->SetDimModel(dim);
    //number of state variables in the problem
    constexpr int nState = 1;


    std::set<int> volmats;
    for(auto regioninfo : data.matinfovec){
      auto matid = std::get<0>(regioninfo);
      auto *dummyMat = new TPZNullMaterial<CSTATE>(matid,dim,nState);
      cmeshH1->InsertMaterialObject(dummyMat);
      volmats.insert(matid);
    }
  
    for(auto pml : pmlDataVec){
      const auto matid = pml.id;
      auto *dummyMat = new TPZNullMaterial<CSTATE>(matid,dim,nState);
      cmeshH1->InsertMaterialObject(dummyMat);
    }



    /**let us associate each boundary with a given material.
       this is important for any non-homogeneous BCs*/
    for(auto &bc : bcDataVec){
      auto res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, bc.id, volmats);
      if(!res.has_value()){
        std::cout<<__PRETTY_FUNCTION__
                 <<"\nwarning: could not find neighbour of bc "<<bc.id<<std::endl;
      }
      bc.volid = res.value();
    }
  
    TPZFNMatrix<1, CSTATE> val1(1, 1, 0);
    TPZManVector<CSTATE,1> val2(1, 0.);
    for(auto bc : bcDataVec){
      const int bctype = wgma::bc::to_int(bc.t);
      const int id = bc.id;
      const int volid = bc.volid;
      auto *dummyMat =
        dynamic_cast<TPZMaterialT<CSTATE>*>(cmeshH1->FindMaterial(volid));
      auto *dummyBC = dummyMat->CreateBC(dummyMat, id, bctype, val1, val2);
      cmeshH1->InsertMaterialObject(dummyBC);
    }

    cmeshH1->SetAllCreateFunctionsContinuous();
    cmeshH1->AutoBuild();
    cmeshH1->CleanUpUnconnectedNodes();

    /*
      Then we create the computational mesh associated with the HCurl space
    */
    TPZAutoPointer<TPZCompMesh> cmeshHCurl =
      new TPZCompMesh(gmesh,isComplex);
    cmeshHCurl->SetDefaultOrder(pOrder);
    cmeshHCurl->SetDimModel(dim);
  
    for(auto regioninfo : data.matinfovec){
      auto matid = std::get<0>(regioninfo);
      auto *dummyMat = new TPZNullMaterial<CSTATE>(matid,dim,nState);
      cmeshHCurl->InsertMaterialObject(dummyMat);
    }
    for(auto pml : pmlDataVec){
      const auto dummyMatid = pml.id;
      auto *dummyMat = new TPZNullMaterial<CSTATE>(dummyMatid,dim,nState);
      cmeshHCurl->InsertMaterialObject(dummyMat);
    }

  
    for(auto bc : bcDataVec){
      const int bctype = wgma::bc::to_int(bc.t);
      const int id = bc.id;
      const int volid = bc.volid;
      auto *dummyMat =
        dynamic_cast<TPZMaterialT<CSTATE>*>(cmeshHCurl->FindMaterial(volid));
      auto *dummyBC = dummyMat->CreateBC(dummyMat, id, bctype, val1, val2);
      cmeshHCurl->InsertMaterialObject(dummyBC);
    }

    cmeshHCurl->SetAllCreateFunctionsHCurl();
    cmeshHCurl->AutoBuild();
    cmeshHCurl->CleanUpUnconnectedNodes();

  
    TPZAutoPointer<TPZCompMesh> cmeshMF =
      new TPZCompMesh(gmesh,isComplex);
    for(auto [matid, er, ur] : data.matinfovec){
      auto *matWG = new TPZWgma(matid, er, ur, lambda, scale);
      cmeshMF->InsertMaterialObject(matWG);
    }
  
    //insert PML regions
    for(auto pml : pmlDataVec){
      const auto id = pml.id;
      const auto alphax = pml.alphax;
      const auto alphay = pml.alphay;
      const auto type = pml.t;
      cmeshtools::AddRectangularPMLRegion<
        TPZWgma
        >(id, alphax, alphay, type, volmats, gmesh, cmeshMF);
    }
  
    for(auto bc : bcDataVec){
      const int bctype = wgma::bc::to_int(bc.t);
      const int id = bc.id;
      const int volid = bc.volid;
      auto *matWG =
        dynamic_cast<TPZMaterialT<CSTATE>*>(cmeshMF->FindMaterial(volid));
      auto *bcMat = matWG->CreateBC(matWG, id, bctype, val1, val2);
      cmeshMF->InsertMaterialObject(bcMat);
    }

    cmeshMF->SetDimModel(dim);
    cmeshMF->SetAllCreateFunctionsMultiphysicElem();

    cmeshMF->AutoBuild();
    cmeshMF->CleanUpUnconnectedNodes();

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
        new TPZPeriodicWgma(id, dim, er, ur, lambda, matmode, scale);
      cmeshH1->InsertMaterialObject(mat);
      // for pml
      volmats.insert(id);
      //for assembling only desired materials
      allmats.insert(id);
    }

    for (auto pml : data.pmlvec) {
      const auto id = pml.id;
      const auto alphax = pml.alphax;
      const auto alphay = pml.alphay;
      const auto type = pml.t;
      wgma::cmeshtools::AddRectangularPMLRegion<TPZPeriodicWgma>(
        id, alphax, alphay, type, volmats, gmesh, cmeshH1);
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