#include "cmeshtools.hpp"
#include "pmltypes.hpp"
#include "cmeshtools_impl.hpp"

#include <pzgmesh.h>
#include <pzcmesh.h>
#include <TPZNullMaterial.h>
#include <Electromagnetics/TPZWaveguideModalAnalysis.h>
#include <Electromagnetics/TPZWaveguideModalAnalysisPML.h>
#include <pzbuildmultiphysicsmesh.h>
#include <TPZSimpleTimer.h>


#include <regex>//for string search

using namespace wgma;


void
cmeshtools::SetupGmshMaterialData(
  const TPZVec<std::map<std::string,int>> &gmshmats,
  const std::map<std::string,std::pair<CSTATE,CSTATE>> &matmap,
  const std::map<std::string,wgma::bc::type> &bcmap,
  const STATE alphaPMLx,
  const STATE alphaPMLy,
  wgma::cmeshtools::PhysicalData &data)
{
  auto &matinfo = data.matinfovec;
  auto &pmlvec = data.pmlvec;
  auto &bcvec = data.bcvec;
  pmlvec.resize(0);
    
  for(auto mat : gmshmats[2]){
    const std::string name = mat.first;
    const auto id = mat.second;
    constexpr auto pmlname{"pml"};

    //std::regex_constants::icase - ignores case
    const auto rx = std::regex{ pmlname ,std::regex_constants::icase };
    const bool ispml = std::regex_search(name, rx);
      
    if(ispml){
      
      const auto pos = pmlvec.size();
      bool found{false};
      const TPZVec<std::string> pmlnames =
        {"xpyp", "xmyp", "xmym", "xpym", "xp", "yp", "xm", "ym"};
      const TPZVec<wgma::pml::type> pmltypes =
        {wgma::pml::type::xpyp,wgma::pml::type::xmyp,
         wgma::pml::type::xmym,wgma::pml::type::xpym,
         wgma::pml::type::xp,wgma::pml::type::yp,
         wgma::pml::type::xm,wgma::pml::type::ym};
      
      for(int ipml = 0; ipml < pmlnames.size(); ipml++){
        const std::string pmlname = pmlnames[ipml];
        const auto rx = std::regex{ pmlname, std::regex_constants::icase };
        const bool test = std::regex_search(name, rx);
        if(test){
          pmlvec.resize(pos+1);
          pmlvec[pos].id = id;
          pmlvec[pos].alphax = alphaPMLx;
          pmlvec[pos].alphay = alphaPMLy;
          pmlvec[pos].t = pmltypes[ipml];
          found = true;
          break;
        }
      }
      //pml was not identified
      if(!found){
        std::cout<<"error: mat "<<name<<" id "<<id<<" not found"
                 <<"\nAborting..."<<std::endl;
      }
    }else{
      if(matmap.find(name) == matmap.end()){//material not found
        std::cout<<"error: mat "<<name<<" id "<<id<<" not found"
                 <<"\nAborting..."<<std::endl;
        DebugStop();
      }else{
        const auto pos = matinfo.size();
        matinfo.push_back(std::make_tuple(
                            id,
                            matmap.at(name).first,
                            matmap.at(name).second)
                          );
      }
    }
  }
  //bc materials
  {
    for(auto bc : gmshmats[1]){
      const auto name = bc.first;
      const auto id = bc.second;
      /*sometimes we need 1d materials with a given id for other purposes,
        such as circumference arcs. soo, not finding it is not a problem*/
      if(bcmap.find(name) != bcmap.end()){
        const int ibc = bcvec.size();
        bcvec.resize(ibc+1);
        bcvec[ibc].id = id;
        bcvec[ibc].t = bcmap.at(name);
      }
    }
  }
    
}


TPZVec<TPZAutoPointer<TPZCompMesh>>
cmeshtools::CreateCMesh(
  TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
  PhysicalData &data, const STATE lambda, const REAL &scale)
{
  TPZSimpleTimer timer ("Create cmesh");
  constexpr int dim = 2;
  constexpr bool isComplex{true};
  auto &pmlDataVec = data.pmlvec;
  auto &bcDataVec = data.bcvec;
  const int nVolMats = data.matinfovec.size();
  const int nPmlMats = pmlDataVec.size();
  const int nBcMats = bcDataVec.size();
  
  /**let us associate each boundary with a given material.
     this is important for any non-homogeneous BCs*/
  for(auto &bc : bcDataVec){
    for(auto *gel : gmesh->ElementVec()){
      if(gel->MaterialId() == bc.id){
        const auto maxside = gel->NSides() - 1;
        bc.volid = gel->Neighbour(maxside).Element()->MaterialId();
        break;
      }
    }
  }
  /*
   First we create the computational mesh associated with the H1 space
   (ez component)*/
  TPZAutoPointer<TPZCompMesh> cmeshH1 =
    new TPZCompMesh(gmesh,isComplex);
  cmeshH1->SetDefaultOrder(pOrder +1);//for deRham compatibility
  cmeshH1->SetDimModel(dim);
  //number of state variables in the problem
  constexpr int nState = 1;
  for(auto regioninfo : data.matinfovec){
    auto matid = std::get<0>(regioninfo);
    auto *dummyMat = new TPZNullMaterial<CSTATE>(matid,dim,nState);
    cmeshH1->InsertMaterialObject(dummyMat);
  }
  for(auto pml : pmlDataVec){
    const auto matid = pml.id;
    auto *dummyMat = new TPZNullMaterial<CSTATE>(matid,dim,nState);
    cmeshH1->InsertMaterialObject(dummyMat);
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
  std::set<int> volmats;
  for(auto [matid, er, ur] : data.matinfovec){
    volmats.insert(matid);
  }
  for(auto pml : pmlDataVec){
    const auto id = pml.id;
    const auto alphax = pml.alphax;
    const auto alphay = pml.alphay;
    const auto type = pml.t;
    AddRectangularPMLRegion<
      TPZWaveguideModalAnalysisPML,TPZWaveguideModalAnalysis
      >(id, alphax, alphay, type, volmats, gmesh, cmeshMF);
  }
  
  for(auto bc : bcDataVec){
    const int bctype = wgma::bc::to_int(bc.t);
    const int id = bc.id;
    const int volid = bc.volid;
    auto *matWG =
      dynamic_cast<TPZMaterialT<CSTATE>*>(cmeshH1->FindMaterial(volid));
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


int
cmeshtools::FindPMLNeighbourMaterial(
  TPZAutoPointer<TPZGeoMesh> gmesh,const int pmlId,
  const std::set<int> &volmats,
  const REAL boundPosX, const REAL boundPosY)
{
  TPZGeoEl * closestEl = nullptr;
  REAL dist = 1e16;
  for(auto &currentEl : gmesh->ElementVec()){
    if ( !currentEl ||
         currentEl->NSubElements() > 0  ||
         currentEl->Dimension() != 2 ||
         volmats.count(currentEl->MaterialId()) == 0) continue;
    TPZVec<REAL> qsi(2,-1);
    const int largerSize = currentEl->NSides() - 1;
    currentEl->CenterPoint(largerSize, qsi);
    TPZVec<REAL> xCenter(3,-1);
    currentEl->X(qsi, xCenter);
    const REAL currentDist = (xCenter[0]-boundPosX)*(xCenter[0]-boundPosX) +
      (xCenter[1]-boundPosY)*(xCenter[1]-boundPosY);
    if(currentDist < dist){
      dist = currentDist;
      closestEl = currentEl;
    }
  }

  if(!closestEl){
    PZError<<"Could not find pml neighbour, aborting...."<<std::endl;
    DebugStop();
  }
  return closestEl->MaterialId();
}

void
cmeshtools::CountActiveEquations(TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec,
                                 const std::set<int64_t> &boundConnects,
                                 int &neq,
                                 int &nH1Equations, int &nHCurlEquations)
{
  auto cmesh = meshVec[0];
  neq = nH1Equations = nHCurlEquations = 0;
  auto cmeshHCurl = meshVec[1];
  auto cmeshH1 = meshVec[2];
  
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

void
cmeshtools::FilterBoundaryEquations(TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec,
                                    TPZVec<int64_t> &activeEquations, int &neq,
                                    int &neqOriginal,
                                    int &nh1, int &nhcurl)
{
  TPZSimpleTimer timer ("Filter dirichlet eqs");
  auto cmesh = meshVec[0];
  TPZManVector<int64_t, 1000> allConnects;
  std::set<int64_t> boundConnects;

  for (int iel = 0; iel < cmesh->NElements(); iel++) {
    TPZCompEl *cel = cmesh->ElementVec()[iel];
    if (cel == nullptr) {
      continue;
    }
    if (cel->Reference() == nullptr) {//there is no associated geometric el
      continue;
    }
    TPZBndCond *mat = dynamic_cast<TPZBndCond *>(
        cmesh->MaterialVec()[cel->Reference()->MaterialId()]);
    if (mat && mat->Type() == 0) {//check for dirichlet bcs
      std::set<int64_t> boundConnectsEl;
      std::set<int64_t> depBoundConnectsEl;
      std::set<int64_t> indepBoundConnectsEl;
      cel->BuildConnectList(indepBoundConnectsEl, depBoundConnectsEl);
      cel->BuildConnectList(boundConnectsEl);
      for(auto val : boundConnectsEl){
        if (boundConnects.find(val) == boundConnects.end()) {
          boundConnects.insert(val);
        }
      }
    }
  }
  for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
    if (boundConnects.find(iCon) == boundConnects.end()) {
      TPZConnect &con = cmesh->ConnectVec()[iCon];
      if (con.HasDependency())
        continue;
      const auto seqnum = con.SequenceNumber();
      const auto pos = cmesh->Block().Position(seqnum);
      const auto blocksize = cmesh->Block().Size(seqnum);
      if (blocksize == 0)
        continue;

      const auto vs = activeEquations.size();
      activeEquations.Resize(vs + blocksize);
      for (auto ieq = 0; ieq < blocksize; ieq++) {
        activeEquations[vs + ieq] = pos + ieq;
      }
    }
  }

  neqOriginal = cmesh->NEquations();
  CountActiveEquations(meshVec, boundConnects, neq, nh1, nhcurl);
}