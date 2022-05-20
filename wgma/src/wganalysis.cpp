#include "wganalysis.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "cmeshtools_impl.hpp"

#include <TPZSpStructMatrix.h>
#include <TPZKrylovEigenSolver.h>
#include <Electromagnetics/TPZWaveguideModalAnalysis.h>
#include <TPZNullMaterial.h>
#include <TPZSimpleTimer.h>
#include <pzbuildmultiphysicsmesh.h>

namespace wgma::wganalysis{
  
  
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
    m_cmesh_h1 = meshvec[TPZWaveguideModalAnalysis::H1Index()];
    m_cmesh_hcurl = meshvec[TPZWaveguideModalAnalysis::HCurlIndex()];
  
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
      const auto firstHCurl = m_n_dofs_h1 * TPZWaveguideModalAnalysis::HCurlIndex();
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
        if (TPZWaveguideModalAnalysis::H1Index() == 0 && iCon < cmeshH1->NConnects()) {
          isH1 = true;
        } else if (TPZWaveguideModalAnalysis::H1Index() == 1 && iCon >= cmeshHCurl->NConnects()) {
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
    scalnames.Push("Ez");
    vecnames.Push("Et");
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
          dynamic_cast<TPZWaveguideModalAnalysis *>(m_cmesh_mf->FindMaterial(id));
        if(!matPtr) continue;
        matPtr->SetKz(currentKz);
        matPtr->SetPrintFieldRealPart(print_real_part);
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
  
    TPZFNMatrix<1, CSTATE> val1(1, 1, 1);
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
      auto *matWG = new TPZWaveguideModalAnalysis(matid, er, ur, lambda, scale);
      cmeshMF->InsertMaterialObject(matWG);
    }
  
    //insert PML regions
    for(auto pml : pmlDataVec){
      const auto id = pml.id;
      const auto alphax = pml.alphax;
      const auto alphay = pml.alphay;
      const auto type = pml.t;
      cmeshtools::AddRectangularPMLRegion<
        TPZWaveguideModalAnalysis
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
    meshVecIn[TPZWaveguideModalAnalysis::H1Index()] = cmeshH1.operator->();
    meshVecIn[TPZWaveguideModalAnalysis::HCurlIndex()] = cmeshHCurl.operator->();

  
    TPZBuildMultiphysicsMesh::AddElements(meshVecIn, cmeshMF.operator->());
    TPZBuildMultiphysicsMesh::AddConnects(meshVecIn, cmeshMF.operator->());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshVecIn, cmeshMF.operator->());

    cmeshMF->ExpandSolution();
    cmeshMF->ComputeNodElCon();
    cmeshMF->CleanUpUnconnectedNodes();

    TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec(3,nullptr);
    meshVec[0] = cmeshMF;
    meshVec[1 + TPZWaveguideModalAnalysis::H1Index()] = cmeshH1;
    meshVec[1 + TPZWaveguideModalAnalysis::HCurlIndex()] = cmeshHCurl;
    return meshVec;
  
  }

};