#include "cmeshtools.hpp"
#include "pmltypes.hpp"

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
  const STATE alphaPML,
  TPZVec<int> &volmatids,
  TPZVec<CSTATE> &ervec,
  TPZVec<CSTATE> &urvec,
  TPZVec<wgma::pml::data> &pmlvec,
  TPZVec<wgma::bc::data> &bcvec)
{
  volmatids.Resize(0);
  urvec.Resize(0);
  ervec.Resize(0);
  pmlvec.Resize(0);
    
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
          pmlvec.Resize(pos+1);
          pmlvec[pos].id = id;
          pmlvec[pos].alpha = alphaPML;
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
        const auto pos = volmatids.size();
        volmatids.Resize(pos+1);
        urvec.Resize(pos+1);
        ervec.Resize(pos+1);
        volmatids[pos] = id;
        ervec[pos] = matmap.at(name).first;
        urvec[pos] = matmap.at(name).second;
      }
    }
  }
  //bc materials
  const int nbcs = gmshmats[1].size();
  bcvec.Resize(nbcs);

  {
    int ibc = 0;
    for(auto bc : gmshmats[1]){
      const auto name = bc.first;
      const auto id = bc.second;
      if(bcmap.find(name) == bcmap.end()){//material not found
        std::cout<<"error: bc "<<name<<" id "<<id<<" not found"
                 <<"\nAborting..."<<std::endl;
        DebugStop();
      }else{
        bcvec[ibc].id = id;
        bcvec[ibc].t = bcmap.at(name);
        ibc++;
      }
    }
  }
    
}


TPZVec<TPZAutoPointer<TPZCompMesh>>
cmeshtools::CreateCMesh(
  TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
  const TPZVec<int> &volMatIdVec, const TPZVec<CSTATE> &urVec,
  const TPZVec<CSTATE> &erVec, const TPZVec<pml::data> &pmlDataVec,
  const TPZVec<bc::data> &bcDataVec, const STATE lambda, const REAL &scale)
{
  TPZSimpleTimer timer ("Create cmesh");
  constexpr int dim = 2;
  constexpr bool isComplex{true};

  const int nVolMats = volMatIdVec.size();
  const int nPmlMats = pmlDataVec.size();
  const int nBcMats = bcDataVec.size();
  

  /*
   First we create the computational mesh associated with the H1 space
   (ez component)*/
  auto * cmeshH1 =new TPZCompMesh(gmesh,isComplex);
  cmeshH1->SetDefaultOrder(pOrder +1);//for deRham compatibility
  cmeshH1->SetDimModel(dim);
  //number of state variables in the problem
  constexpr int nState = 1;
  TPZMaterialT<CSTATE> *dummyMat{nullptr};
  for(auto matid : volMatIdVec){
    dummyMat = new TPZNullMaterial<CSTATE>(matid,dim,nState);
    cmeshH1->InsertMaterialObject(dummyMat);
  }
  for(auto pml : pmlDataVec){
    const auto matid = pml.id;
    dummyMat = new TPZNullMaterial<CSTATE>(matid,dim,nState);
    cmeshH1->InsertMaterialObject(dummyMat);
  }

  
  TPZFNMatrix<1, CSTATE> val1(1, 1, 1);
  TPZManVector<CSTATE,1> val2(1, 0.);
  TPZBndCond *dummyBC = nullptr;
  for(auto bc : bcDataVec){
    const int bctype = wgma::bc::to_int(bc.t);
    const int id = bc.id;
    dummyBC = dummyMat->CreateBC(dummyMat, id, bctype, val1, val2);
    cmeshH1->InsertMaterialObject(dummyBC);
  }

  cmeshH1->SetAllCreateFunctionsContinuous();
  cmeshH1->AutoBuild();
  cmeshH1->CleanUpUnconnectedNodes();

  /*
    Then we create the computational mesh associated with the HCurl space
   */
  auto *cmeshHCurl = new TPZCompMesh(gmesh,isComplex);
  cmeshHCurl->SetDefaultOrder(pOrder);
  cmeshHCurl->SetDimModel(dim);
  
  for(auto matid : volMatIdVec){
    dummyMat = new TPZNullMaterial<CSTATE>(matid,dim,nState);
    cmeshHCurl->InsertMaterialObject(dummyMat);
  }
  for(auto pml : pmlDataVec){
    const auto matid = pml.id;
    dummyMat = new TPZNullMaterial<CSTATE>(matid,dim,nState);
    cmeshHCurl->InsertMaterialObject(dummyMat);
  }

  
  dummyBC = nullptr;
  for(auto bc : bcDataVec){
    const int bctype = wgma::bc::to_int(bc.t);
    const int id = bc.id;
    dummyBC = dummyMat->CreateBC(dummyMat, id, bctype, val1, val2);
    cmeshHCurl->InsertMaterialObject(dummyBC);
  }

  cmeshHCurl->SetAllCreateFunctionsHCurl();
  cmeshHCurl->AutoBuild();
  cmeshHCurl->CleanUpUnconnectedNodes();

  
  auto *cmeshMF =new TPZCompMesh(gmesh,isComplex);
  TPZWaveguideModalAnalysis *matWG = nullptr;
  for(auto i = 0; i < nVolMats; i++){
    matWG = new TPZWaveguideModalAnalysis(
      volMatIdVec[i], urVec[i], erVec[i], lambda, 1. / scale);
    cmeshMF->InsertMaterialObject(matWG);
  }
  
  //insert PML regions
  for(auto pml : pmlDataVec){
    const auto id = pml.id;
    const auto alpha = pml.alpha;
    const auto type = pml.t;
    AddRectangularPMLRegion(id, alpha, type, cmeshMF);
  }
  
  TPZBndCond *bcMat = nullptr;
  for(auto bc : bcDataVec){
    const int bctype = wgma::bc::to_int(bc.t);
    const int id = bc.id;
    bcMat = matWG->CreateBC(matWG, id, bctype, val1, val2);
    cmeshMF->InsertMaterialObject(bcMat);
  }

  cmeshMF->SetDimModel(dim);
  cmeshMF->SetAllCreateFunctionsMultiphysicElem();

  cmeshMF->AutoBuild();
  cmeshMF->CleanUpUnconnectedNodes();

  TPZManVector<TPZCompMesh*,3> meshVecIn(2);
  meshVecIn[TPZWaveguideModalAnalysis::H1Index()] = cmeshH1;
  meshVecIn[TPZWaveguideModalAnalysis::HCurlIndex()] = cmeshHCurl;

  
  TPZBuildMultiphysicsMesh::AddElements(meshVecIn, cmeshMF);
  TPZBuildMultiphysicsMesh::AddConnects(meshVecIn, cmeshMF);
  TPZBuildMultiphysicsMesh::TransferFromMeshes(meshVecIn, cmeshMF);

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
  const REAL boundPosX, const REAL boundPosY)
{
  TPZGeoEl * closestEl = nullptr;
  REAL dist = 1e16;
  for(auto &currentEl : gmesh->ElementVec()){
    if ( !currentEl ||
         currentEl->NSubElements() > 0  ||
         currentEl->Dimension() != 2 ||
         currentEl->MaterialId() == pmlId) continue;
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
cmeshtools::AddRectangularPMLRegion(const int matId, const int alpha,
                                    const wgma::pml::type type,
                                    TPZAutoPointer<TPZCompMesh> cmesh)
{

  //let us find the (xmin,xmax) and (ymin,ymax) of the PML region
  REAL xMax = -1e20, xMin = 1e20, yMax = -1e20, yMin = 1e20;
  TPZGeoMesh *gmesh = cmesh->Reference();
  for (auto geo : gmesh->ElementVec()){
    if (geo->MaterialId() == matId) {
      for (int iNode = 0; iNode < geo->NCornerNodes(); ++iNode) {
        TPZManVector<REAL, 3> co(3);
        geo->Node(iNode).GetCoordinates(co);
        const REAL &xP = co[0];
        const REAL &yP = co[1];
        if (xP > xMax) {
          xMax = xP;
        }
        if (xP < xMin) {
          xMin = xP;
        }
        if (yP > yMax) {
          yMax = yP;
        }
        if (yP < yMin) {
          yMin = yP;
        }
      }
    }
  }


  //now we compute xBegin, yBegin, attx, atty and d for the material ctor
  const bool attx = wgma::pml::attx(type);
  const bool atty = wgma::pml::atty(type);
  
  REAL xBegin{-1}, yBegin{-1}, dX{-01101991.}, dY{-01101991.};
  REAL boundPosX{-01101991.}, boundPosY{-01101991.};
  if(attx){
    const int xdir = wgma::pml::xinfo(type);
    dX = xMax - xMin;
    xBegin = xdir > 0 ? xMin : xMax;
    boundPosX = xBegin;
    boundPosY = (yMax + yMin)/2;
  }

  if(atty){
    const int ydir = wgma::pml::yinfo(type);
    dY = yMax - yMin;
    yBegin = ydir > 0 ? yMin : yMax;
    boundPosX = (xMax + xMin)/2;
    boundPosY = yBegin;
  }

  if(attx && atty){
    boundPosX = xBegin;
    boundPosY = yBegin;
  }
  //find the neighbouring material
  const auto neighMatId =
    FindPMLNeighbourMaterial(gmesh, matId, boundPosX, boundPosY);

  auto neighMat = dynamic_cast<TPZWaveguideModalAnalysis*>(
    cmesh->FindMaterial(neighMatId));
  if(!neighMat){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\n neighbouring material not found in mesh, aborting...."<<std::endl;
    DebugStop();
  }
  
  auto pmlMat = new TPZWaveguideModalAnalysisPML(matId, *neighMat);
  if(attx) pmlMat->SetAttX(xBegin, alpha, dX);
  if(atty) pmlMat->SetAttY(yBegin, alpha, dY);
  cmesh->InsertMaterialObject(pmlMat);
  
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
