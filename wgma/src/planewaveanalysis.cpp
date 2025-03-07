#include "planewaveanalysis.hpp"
#include "gmeshtools.hpp"
#include "materials/planewavesolutions.hpp"
#include "materials/solutionprojection.hpp"


#include <TPZSpStructMatrix.h>
#include <TPZSSpStructMatrix.h>
#include <TPZStructMatrixOMPorTBB.h>
#include <pzstrmatrixot.h>
#include <TPZSimpleTimer.h>
#include <TPZNullMaterialCS.h>
#include <pzbuildmultiphysicsmesh.h>

#include <numeric>
#include <cassert>

void RenumberMultiphysicsMesh(TPZAutoPointer<TPZCompMesh>m_cmesh_h1,
                              TPZAutoPointer<TPZCompMesh>m_cmesh_hcurl,
                              TPZAutoPointer<TPZCompMesh>m_cmesh_mf);

namespace wgma::planewaveanalysis{
  bool using_tbb_mat{false};
  
  Analysis::Analysis(const TPZVec<TPZAutoPointer<TPZCompMesh>> &meshvec,
                     const int n_threads, const bool reorder_eqs,
                     const bool filter_bound) :
    TPZLinearAnalysis(),
    m_filter_bound(filter_bound)
  {
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

    for (auto [id,mat] : m_cmesh_mf->MaterialVec()){
      auto matplanewave =
        dynamic_cast<wgma::materials::PlaneWaveSolutions*>(mat);
      if(matplanewave){
        matplanewave->GetBeta(m_eigenvalues);
        break;
      }
    }
    //we do not reorder eqs on multiphysics mesh
    this->SetCompMesh(m_cmesh_mf.operator->(), false);
    if(reorder_eqs){
      RenumberMultiphysicsMesh(m_cmesh_h1, m_cmesh_hcurl, m_cmesh_mf);
    }

    TPZAutoPointer<TPZStructMatrix> strmtrx{nullptr};
    if(using_tbb_mat){
      auto mtrx = new TPZSpStructMatrix<CSTATE,TPZStructMatrixOMPorTBB<CSTATE>>(m_cmesh_mf);
      mtrx->SetShouldColor(false);
      mtrx->SetTBBorOMP(true);
      strmtrx = mtrx;
    }else{
      auto mtrx = new TPZSpStructMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>(m_cmesh_mf);
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
    std::cout << "------\tactive eqs\t-------" << std::endl;
    std::cout << "# H1 equations: " << m_n_dofs_h1 << std::endl;
    std::cout << "# HCurl equations: " << m_n_dofs_hcurl << std::endl;
    std::cout << "# equations: " << m_n_dofs_mf << std::endl;
    std::cout << "------\t----------\t-------" << std::endl;
    this->SetStructuralMatrix(strmtrx);
  }

  void Analysis::LoadSolution(){
    TPZAnalysis::LoadSolution();
    TPZManVector<TPZAutoPointer<TPZCompMesh>,2> meshVecPost(2);
    meshVecPost[TPZWgma::H1Index()] = m_cmesh_h1;
    meshVecPost[TPZWgma::HCurlIndex()] = m_cmesh_hcurl;
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshVecPost, m_cmesh_mf);
  }
  void Analysis::Assemble(){
    TPZSimpleTimer assemble("Assemble");
    //assembles the system
    TPZLinearAnalysis::Assemble();
  }
  
  void Analysis::AssembleRhs(std::set<int> matids){
    auto strmat = this->StructMatrix();
    auto matids_cp = strmat->MaterialIds();
    strmat->SetMaterialIds(matids);
    TPZLinearAnalysis::AssembleResidual();
    strmat->SetMaterialIds(matids_cp);
  }
  
  void Analysis::Solve(){
    TPZSimpleTimer solve("Solve");
    TPZLinearAnalysis::Solve();
  }
  void Analysis::Run(){
    Assemble();
    Solve();
  }
  
  void
  Analysis::CountActiveEqs(int64_t &neq,int64_t &nH1Equations, int64_t &nHCurlEquations)
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
        if(seqnum < 0){continue;}
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
    return;
  }

  TPZAutoPointer<TPZCompMesh> CreateMfPlaneWaveMats(TPZAutoPointer<TPZGeoMesh> gmesh,
                                                    cmeshtools::PhysicalData &data,
                                                    const STATE lambda,
                                                    const STATE lx,
                                                    const STATE ly,
                                                    const int max_k,
                                                    const REAL &scale,
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
      auto *matWG = new wgma::materials::PlaneWaveSolutions(
        matid,lambda,er,lx,ly,max_k);
      cmeshMF->InsertMaterialObject(matWG);
      realvolmats.insert(matid);
      volmats.insert(matid);
      if(verbose){
        std::cout<<"\t id "<<matid<<" er "<<er<<" ur "<<ur<<'\n';
      }
    }
    //insert PML regions
    for(auto &pml : pmlDataVec){
      PZError<<__PRETTY_FUNCTION__
             <<"\nThis is not supported for now! Aborting..."
             <<std::endl;
      DebugStop();
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
  CMeshPlaneWave2D(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
                   cmeshtools::PhysicalData &data,
                   const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &el_map,
                   const STATE lambda,
                   const STATE lx,
                   const STATE ly,
                   const int max_k,
                   const REAL &scale,
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

    for(auto periodic_els : el_map){
      wgma::cmeshtools::SetPeriodic(cmeshH1,periodic_els);
    }
    /*
      Then we create the computational mesh associated with the HCurl space
    */
    ish1 = false;
    TPZAutoPointer<TPZCompMesh> cmeshHCurl =
      CreateAtomicWgma2D(gmesh, ish1,pOrder,volmats,pmlmats, data.bcvec,data.probevec);
    /*
      Now we add the periodicity
    */
    for(auto periodic_els : el_map){
      wgma::cmeshtools::SetPeriodic(cmeshHCurl,periodic_els);
    }
    /*
      Now we create the MF mesh
    */
    TPZAutoPointer<TPZCompMesh> cmeshMF =
      CreateMfPlaneWaveMats(gmesh, data, lambda, lx, ly, max_k, scale, verbose);

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
};

void RenumberMultiphysicsMesh(TPZAutoPointer<TPZCompMesh>m_cmesh_h1,
                              TPZAutoPointer<TPZCompMesh>m_cmesh_hcurl,
                              TPZAutoPointer<TPZCompMesh>m_cmesh_mf)
{
#ifdef PZ_USING_METIS
  const auto renumtype = RenumType::EMetis;
#else
  const auto renumtype = RenumType::ECutHillMcKee;
#endif
  //we create TPZLinearAnalysis objects just to reorder the eq of the atomic meshes
  {
    TPZLinearAnalysis an(m_cmesh_h1,renumtype);
    m_cmesh_h1->ExpandSolution();
  }
  {
    TPZLinearAnalysis an(m_cmesh_hcurl,renumtype);
    m_cmesh_hcurl->ExpandSolution();
  }

  //first we compute the number of independent connects in each mesh
  const auto n_h1_con = m_cmesh_h1->NConnects();
  int64_t n_indep_h1{0};
  int64_t n_h1_eqs{0};
  for(auto &c : m_cmesh_h1->ConnectVec()){
    if (!c.HasDependency() && c.NElConnected() && !c.IsCondensed()) {
      n_indep_h1++;
      n_h1_eqs += c.NDof();
    }
  }
  const int64_t n_dep_h1 = n_h1_con - n_indep_h1;

    
  const auto n_hcurl_con = m_cmesh_hcurl->NConnects();
  int64_t n_indep_hcurl{0};
  int64_t n_hcurl_eqs{0};
  for(auto &c : m_cmesh_hcurl->ConnectVec()){
    if (!c.HasDependency() && c.NElConnected() && !c.IsCondensed()) {
      n_indep_hcurl++;
      n_hcurl_eqs += c.NDof();
    }
  }
  const int64_t n_dep_hcurl = n_hcurl_con - n_indep_hcurl;

    
  const auto first_h1_con = n_hcurl_con*TPZWgma::H1Index();
  const auto first_hcurl_con = n_h1_con*TPZWgma::HCurlIndex();

  const auto first_h1_seqnum = n_indep_hcurl*TPZWgma::H1Index();
  const auto first_hcurl_seqnum = n_indep_h1*TPZWgma::HCurlIndex();

  //seqnumber counter of dependent connects
  int64_t dep_count = n_indep_hcurl + n_indep_h1;
    

  auto AdaptConnects = [m_cmesh_mf, &dep_count](TPZCompMesh *atomic_mesh,
                                                const auto ncon,
                                                const auto first_con,
                                                const auto first_seqnum){
    for(auto ic = 0; ic < ncon;ic++){
      const auto &refc = atomic_mesh->ConnectVec()[ic];
      auto &c = m_cmesh_mf->ConnectVec()[first_con+ic];
      const auto atomic_seqnum = refc.SequenceNumber();
      const auto is_indep_con =
        !c.HasDependency() && c.NElConnected() && !c.IsCondensed();
      /*
        we want to keep the sequence number in two distinct scenarios
        1. the connect is independent, etc, it should go first
        2. the connect has a negative seqnum

        otherwise we sent it to the end
      */
      if (is_indep_con){
        c.SetSequenceNumber(first_seqnum+atomic_seqnum);
      }else if (atomic_seqnum == -1){
        c.SetSequenceNumber(-1);
      }else {
        c.SetSequenceNumber(dep_count++);
      }
      if(c.SequenceNumber() >= 0){
        m_cmesh_mf->Block().Set(c.SequenceNumber(),c.NDof());
      }
    }
  };
    
  AdaptConnects(m_cmesh_h1.operator->(),
                n_h1_con,first_h1_con,first_h1_seqnum);
  AdaptConnects(m_cmesh_hcurl.operator->(),
                n_hcurl_con,first_hcurl_con,first_hcurl_seqnum);

  m_cmesh_mf->InitializeBlock();
  return ;
}