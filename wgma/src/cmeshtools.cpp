#include "cmeshtools.hpp"

#include <pzgmesh.h>
#include <pzcmesh.h>
#include <TPZNullMaterial.h>
#include <Electromagnetics/TPZWaveguideModalAnalysis.h>
#include <pzbuildmultiphysicsmesh.h>
#include <TPZSimpleTimer.h>

using namespace wgma;


TPZVec<TPZAutoPointer<TPZCompMesh>>
cmeshtools::CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
            const TPZVec<int> &matIdVec, const CSTATE ur, const CSTATE er,
            const STATE lambda, const REAL &scale, bool usingSymmetry, bctype sym)
{
  TPZSimpleTimer timer ("Create cmesh");
  constexpr int dim = 2;
  constexpr bool isComplex{true};


  /*
   First we create the computational mesh associated with the H1 space
   (ez component)*/
  auto * cmeshH1 =new TPZCompMesh(gmesh,isComplex);
  cmeshH1->SetDefaultOrder(pOrder +1);//for deRham compatibility
  cmeshH1->SetDimModel(dim);

  const int volMatId = matIdVec[0];
  //number of state variables in the problem
  constexpr int nState = 1;

  auto dummyMat = new TPZNullMaterial<CSTATE>(volMatId,dim,nState);
  cmeshH1->InsertMaterialObject(dummyMat);

  
  TPZFNMatrix<1, CSTATE> val1(1, 1, 1);
  TPZManVector<CSTATE,1> val2(1, 0.);
  const int nMats = matIdVec.size();
  TPZBndCond *dummyBC = nullptr;
  for (int i = 1; i <nMats; i++) {
    //0 for dirichlet (PEC) and 1 for neumann (PMC)
    const int bcType = i==1 && usingSymmetry && sym == bctype::PMC ? 1 : 0;
    dummyBC =
      dummyMat->CreateBC(dummyMat, matIdVec[i], bcType, val1, val2);
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
  
  dummyMat = new TPZNullMaterial<CSTATE>(volMatId,dim,nState);
  cmeshHCurl->InsertMaterialObject(dummyMat);

  
  dummyBC = nullptr;
  for (int i = 1; i <nMats; i++) {
    //0 for dirichlet (PEC) and 1 for neumann (PMC)
    const int bcType = i==1 && usingSymmetry && sym == bctype::PMC ? 1 : 0;
    dummyBC =
      dummyMat->CreateBC(dummyMat, matIdVec[i], bcType, val1, val2);
    cmeshHCurl->InsertMaterialObject(dummyBC);
  }

  cmeshHCurl->SetAllCreateFunctionsHCurl();
  cmeshHCurl->AutoBuild();
  cmeshHCurl->CleanUpUnconnectedNodes();

  
  auto *cmeshMF =new TPZCompMesh(gmesh,isComplex);
  
  TPZWaveguideModalAnalysis *matWG  = new TPZWaveguideModalAnalysis(
      volMatId, ur, er, lambda, 1. / scale);
  cmeshMF->InsertMaterialObject(matWG);

  TPZBndCond *bcMat = nullptr;
  for (int i = 1; i <nMats; i++) {
    //0 for dirichlet (PEC) and 1 for neumann (PMC)
    const int bcType = i==1 && usingSymmetry && sym == bctype::PMC ? 1 : 0;
    bcMat =
      matWG->CreateBC(matWG, matIdVec[i], bcType, val1, val2);
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