#include "gmeshtools.hpp"

#include <pzgmesh.h>
#include <TPZRefPattern.h>
#include <TPZRefPatternDataBase.h>
#include <tpzarc3d.h>
#include <tpzgeoblend.h>


TPZAutoPointer<TPZGeoMesh>
wgma::gmeshtools::CreateStructuredMesh(
  const TPZVec<TPZVec<REAL>> &pointsVec, const TPZFMatrix<int> &quadPointsVec,
  const TPZVec<int> &matIdsQuads,
  const TPZVec<int> &nDivQsi, const TPZVec<int> &nDivEta,
  const TPZVec<bool> &side1NonLinearVec, const TPZVec<bool> &side2NonLinearVec,
  const TPZVec<bool> &side3NonLinearVec, const TPZVec<bool> &side4NonLinearVec,
  const TPZVec<TPZVec<REAL>> &thetaVec, const TPZVec<TPZVec<REAL>> &xcRef,
  const TPZVec<REAL> &rVec, const TPZVec<int> &matIdBoundVec,
  const TPZVec<REAL> &boundDistVec, const bool nonLinearMapping)
{

  auto map_quad_side_arc = [](const TPZVec<REAL> &theta ,const TPZVec<REAL> &xc, const REAL r, const REAL s) {
    TPZVec<REAL> point(2,0.);
    point[0] = xc[0] + r*cos((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
    point[1] = xc[1] + r*sin((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
    return point;
  };
  const int nQuads = side1NonLinearVec.size();

  TPZVec<EdgeData> edgesVec(0);
  TPZVec<QuadrilateralData *> quadVec(nQuads, nullptr);
  int iNonLinear = 0;
  for (int iQuad = 0; iQuad < nQuads; iQuad++) {
    quadVec[iQuad] = new QuadrilateralData(quadPointsVec.GetVal(iQuad, 0), quadPointsVec.GetVal(iQuad, 1),
                                           quadPointsVec.GetVal(iQuad, 2), quadPointsVec.GetVal(iQuad, 3),
                                           pointsVec[quadPointsVec.GetVal(iQuad, 0)],
                                           pointsVec[quadPointsVec.GetVal(iQuad, 1)],
                                           pointsVec[quadPointsVec.GetVal(iQuad, 2)],
                                           pointsVec[quadPointsVec.GetVal(iQuad, 3)]);
    if (side1NonLinearVec[iQuad]) {
      TPZVec<REAL> xc = xcRef[iNonLinear];
      quadVec[iQuad]->mapSide1 =
        [map_quad_side_arc, iNonLinear, thetaVec, xc, rVec](const REAL s)
        {
          return map_quad_side_arc(thetaVec[iNonLinear],xc,rVec[iNonLinear], s);
        };
      iNonLinear++;
    } else if (side2NonLinearVec[iQuad]) {
      TPZVec<REAL> xc = xcRef[iNonLinear];
      quadVec[iQuad]->mapSide2 =
        [map_quad_side_arc, iNonLinear, thetaVec, xc, rVec](const REAL s)
        {
          return map_quad_side_arc(thetaVec[iNonLinear],xc,rVec[iNonLinear], s);
        };
      iNonLinear++;
    } else if (side3NonLinearVec[iQuad]) {
      TPZVec<REAL> xc = xcRef[iNonLinear];
      quadVec[iQuad]->mapSide3 =
        [map_quad_side_arc, iNonLinear, thetaVec, xc, rVec](const REAL s)
        {
          return map_quad_side_arc(thetaVec[iNonLinear],xc,rVec[iNonLinear], s);
        };
      iNonLinear++;
    } else if (side4NonLinearVec[iQuad]) {
      TPZVec<REAL> xc = xcRef[iNonLinear];
      quadVec[iQuad]->mapSide4 =
        [map_quad_side_arc, iNonLinear, thetaVec, xc, rVec](const REAL s)
        {
          return map_quad_side_arc(thetaVec[iNonLinear],xc,rVec[iNonLinear], s);
        };
      iNonLinear++;
    }
    quadVec[iQuad]->SetMapsLinearityAndAdjustEdges(side1NonLinearVec[iQuad], side2NonLinearVec[iQuad],
                                     side3NonLinearVec[iQuad], side4NonLinearVec[iQuad],
                                     iQuad, edgesVec);
    quadVec[iQuad]->CreateQuadrilateralRegion(nDivQsi[iQuad], nDivEta[iQuad]);
  }

  ////////////////////////////////////////CREATE NODES///////////////////////////////////////
  ////////////////////////////////////////FOR INTERIOR///////////////////////////////////////
  ///////////////////////////////////////////POINTS//////////////////////////////////////////
  long nodeId = 0;
  TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
  TPZVec<TPZVec<long>> elNodeFaceVecs(nQuads, TPZVec<long>(0, -1));

  for (int el = 0; el < nQuads; el++) {
    const long nNodesOld = gmesh->NodeVec().NElements();
    const long nNodesEl = (nDivQsi[el] - 2) * (nDivEta[el] - 2);
    gmesh->NodeVec().Resize(nNodesOld + nNodesEl);
    elNodeFaceVecs[el].Resize(nDivEta[el] * nDivQsi[el], -1);
    //interior nodes only
    for (int i = 1; i < nDivEta[el] - 1; i++) {
      for (int j = 1; j < nDivQsi[el] - 1; j++) {
        TPZManVector<REAL, 3> node(3, 0.);
        node[0] = quadVec[el]->facePts(i * nDivQsi[el] + j, 0);
        node[1] = quadVec[el]->facePts(i * nDivQsi[el] + j, 1);
        gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
        elNodeFaceVecs[el][i * nDivQsi[el] + j] = nodeId;
        nodeId++;
      }
    }
  }

  //    std::cout<<"NUMBER OF EDGES "<<edgesVec.size()<<std::endl;
  //    for(int iEdge = 0; iEdge < edgesVec.size(); iEdge++){
  //        std::cout<<"Edge "<<iEdge<<std::setfill(' ');;
  //        std::cout<<"\tquad1:"<<std::setw(2)<<edgesVec[iEdge].quad1<<"\tquad2:"<<std::setw(2)<<edgesVec[iEdge].quad2;
  //        std::cout<<"\tNon-linear mapped:"<<edgesVec[iEdge].isNonLinearMapped;
  //        std::cout<<"\tBoundary:"<<edgesVec[iEdge].amIboundaryEdge<<std::endl;
  //    }
  const int nEdges = edgesVec.size();
  ////////////////////////////////////////CREATE NODES///////////////////////////////////////
  //////////////////////////////////////////FOR EDGE/////////////////////////////////////////
  ///////////////////////////////////////////POINTS//////////////////////////////////////////
  auto sidePos = [](const int side, const int &i, const int &nQsi, const int &nEta, const bool &reverse) {
    const int iP = reverse ? (side % 2 ? nQsi - 1 - i : nEta - 1 - i) : i;
    switch (side) {
    case 1 :
      return iP;
    case 2 :
      return (iP * nQsi + (nQsi - 1));
    case 3 :
      return (nQsi - 1 - iP + nQsi * (nEta - 1));
    case 4 :
      return (nEta - 1 - iP) * nQsi;
    default:
      DebugStop();
      return -1;
    }
  };

  auto createEdgeNodes = [sidePos, nDivEta, nDivQsi, quadVec, gmesh]
    (TPZVec<TPZVec<long>> &elNodeFaceVecs, const long &el1, const long &el2,
     const int &side1, const int &side2, long &nodeId,
     TPZVec<long> &side1pts, TPZVec<long> &side2pts, const bool &revEl2) {
    const int nPts = side1 % 2 ? nDivQsi[el1] : nDivEta[el1];
    long nNodesOld = gmesh->NodeVec().NElements();
    long nNodesEl = nPts;
    for (int i = 0; i < nPts; i++) {
      const int posEl1 = sidePos(side1, i, nDivQsi[el1], nDivEta[el1], false);
      const int posEl2 = sidePos(side2, i, nDivQsi[el2], nDivEta[el2], revEl2);
      if (elNodeFaceVecs[el1][posEl1] != -1) {
        elNodeFaceVecs[el2][posEl2] = elNodeFaceVecs[el1][posEl1];
        continue;
      }
      if (elNodeFaceVecs[el2][posEl2] != -1) {
        elNodeFaceVecs[el1][posEl1] = elNodeFaceVecs[el2][posEl2];
        continue;
      }
      TPZManVector<REAL, 3> node(3, 0.);
      node[0] = quadVec[el1]->facePts(posEl1, 0);
      node[1] = quadVec[el1]->facePts(posEl1, 1);
      long nNodesOld = gmesh->NodeVec().NElements();
      gmesh->NodeVec().Resize(nNodesOld + 1);
      gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
      elNodeFaceVecs[el1][posEl1] = nodeId;
      elNodeFaceVecs[el2][posEl2] = nodeId;
      nodeId++;
    }
    TPZFMatrix<REAL> *nonLinearPts = nullptr;
    if (quadVec[el1]->isSide1nonLinear && side1 == 1) { nonLinearPts = &quadVec[el1]->side1IntPts; }
    else if (quadVec[el1]->isSide2nonLinear && side1 == 2) { nonLinearPts = &quadVec[el1]->side2IntPts; }
    else if (quadVec[el1]->isSide3nonLinear && side1 == 3) { nonLinearPts = &quadVec[el1]->side3IntPts; }
    else if (quadVec[el1]->isSide4nonLinear && side1 == 4) { nonLinearPts = &quadVec[el1]->side4IntPts; }
    long nNonLinPts = 0;
    if (nonLinearPts) {
      nNonLinPts = nonLinearPts->Rows();
      side1pts.Resize(nNonLinPts);
      side2pts.Resize(nNonLinPts);
      nNodesOld = gmesh->NodeVec().NElements();
      nNodesEl = nNonLinPts;
      gmesh->NodeVec().Resize(nNodesOld + nNodesEl);
    }
    for (int i = 0; i < nNonLinPts; i++) {
      const int posEl1 = i;
      const int posEl2 = revEl2 ? nNonLinPts - 1 - i : i;
      TPZManVector<REAL, 3> node(3, 0.);
      node[0] += (*nonLinearPts)(posEl1, 0);
      node[1] += (*nonLinearPts)(posEl1, 1);
      gmesh->NodeVec()[nodeId].Initialize(node, *gmesh);
      side1pts[posEl1] = nodeId;
      side2pts[posEl2] = nodeId;
      nodeId++;
    }
  };

  TPZVec<TPZVec<long>> elNodeSide1Vecs(nQuads, TPZVec<long>(0, -1));
  TPZVec<TPZVec<long>> elNodeSide2Vecs(nQuads, TPZVec<long>(0, -1));
  TPZVec<TPZVec<long>> elNodeSide3Vecs(nQuads, TPZVec<long>(0, -1));
  TPZVec<TPZVec<long>> elNodeSide4Vecs(nQuads, TPZVec<long>(0, -1));
  for (int i = 0; i < nEdges; i++) {
    TPZVec<long> *side1pts = nullptr;
    TPZVec<long> *side2pts = nullptr;
    int quad1 = edgesVec[i].quad1, quad2 = edgesVec[i].quad2,
      side1 = edgesVec[i].side1, side2 = edgesVec[i].side2;
    bool revEl = !(edgesVec[i].amIboundaryEdge);

    switch (side1) {
    case 1:
      side1pts = &elNodeSide1Vecs[quad1];
      break;
    case 2:
      side1pts = &elNodeSide2Vecs[quad1];
      break;
    case 3:
      side1pts = &elNodeSide3Vecs[quad1];
      break;
    case 4:
      side1pts = &elNodeSide4Vecs[quad1];
      break;
    }

    switch (side2) {
    case 1:
      side2pts = &elNodeSide1Vecs[quad2];
      break;
    case 2:
      side2pts = &elNodeSide2Vecs[quad2];
      break;
    case 3:
      side2pts = &elNodeSide3Vecs[quad2];
      break;
    case 4:
      side2pts = &elNodeSide4Vecs[quad2];
      break;
    }
    createEdgeNodes(elNodeFaceVecs, quad1, quad2, side1, side2, nodeId,
                    *side1pts, *side2pts, revEl);
  }

  gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
  const auto refpQuad = gRefDBase.GetUniformRefPattern(EQuadrilateral);
  for (int quad = 0; quad < nQuads; quad++) {
    for (int i = 0; i < nDivEta[quad] - 1; i++) {
      for (int j = 0; j < nDivQsi[quad] - 1; j++) {

        const int node0 = i * nDivQsi[quad] + j;//Lower left vertex
        const int node1 = i * nDivQsi[quad] + j + 1;//Lower right vertex
        const int node2 = (i + 1) * nDivQsi[quad] + j + 1;//Upper right vertex
        const int node3 = (i + 1) * nDivQsi[quad] + j;//Upper left vertex

        TPZManVector<long, 4> nodesIdVec(4, -1);

        const int matId = matIdsQuads[quad];
        nodesIdVec[0] = elNodeFaceVecs[quad][node0];
        nodesIdVec[1] = elNodeFaceVecs[quad][node1];
        nodesIdVec[2] = elNodeFaceVecs[quad][node2];
        nodesIdVec[3] = elNodeFaceVecs[quad][node3];
        if (nodesIdVec[0] == -1 || nodesIdVec[1] == -1 || nodesIdVec[2] == -1 || nodesIdVec[3] == -1) {
          DebugStop();
        }
        TPZGeoEl *quadEl{nullptr};
        if ((i == 0 && quadVec[quad]->isSide1nonLinear) ||
            (j == nDivQsi[quad] - 2 && quadVec[quad]->isSide2nonLinear) ||
            (i == nDivEta[quad] - 2 && quadVec[quad]->isSide3nonLinear) ||
            (j == 0 && quadVec[quad]->isSide4nonLinear)) {
          quadEl =
            new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(nodesIdVec, matId,
                                                                             *gmesh);
        } else {
          quadEl =
            new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec, matId, *gmesh);
        }
        quadEl->SetRefPattern(refpQuad);
      }
    }
  }
  gRefDBase.InitializeUniformRefPattern(EOned);
  auto refpArc = gRefDBase.GetUniformRefPattern(EOned);
  for (int edge = 0; edge < nEdges; edge++) {
    int quad = edgesVec[edge].quad1, side = edgesVec[edge].side1;
    if (edgesVec[edge].amIboundaryEdge) {//boundary edge
      const int nArcs = side % 2 ? nDivQsi[quad] - 1 : nDivEta[quad] - 1;
      TPZVec<REAL> &pt1 = pointsVec[edgesVec[edge].coord1];
      TPZVec<REAL> &pt2 = pointsVec[edgesVec[edge].coord2];
      REAL tol = 1e-08;
      int matId = -666;
      if (std::abs(pt1[1] - pt2[1]) < tol) {//horizontal edge
        if (std::abs(pt1[1] - boundDistVec[2]) < tol) {
          matId = matIdBoundVec[2];
        } else {
          matId = matIdBoundVec[0];
        }
      } else if (std::abs(pt1[0] - pt2[0]) < tol) {//vertical edge
        if (std::abs(pt1[0] - boundDistVec[1]) < tol) {
          matId = matIdBoundVec[1];
        } else {
          matId = matIdBoundVec[3];
        }
      } else {
        DebugStop();
      }
      for (int i = 0; i < nArcs; i++) {
        TPZManVector<long, 3> nodesIdVec(2, 0.);
        const int vertex1 = sidePos(side, i, nDivQsi[quad], nDivEta[quad], false);
        const int vertex2 = sidePos(side, i + 1, nDivQsi[quad], nDivEta[quad], false);


        nodesIdVec[0] = elNodeFaceVecs[quad][vertex1];
        nodesIdVec[1] = elNodeFaceVecs[quad][vertex2];
        if (nodesIdVec[0] == -1 || nodesIdVec[1] == -1) {
          DebugStop();
        }

        TPZGeoElRefPattern<pzgeom::TPZGeoLinear> *arc =
          new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matId, *gmesh);
        arc->SetRefPattern(refpArc);
      }
      continue;
    }

    TPZVec<long> *intPointsCoord = nullptr;
    switch (side) {
    case 1:
      if (quadVec[quad]->isSide1nonLinear == false) continue;
      intPointsCoord = &(elNodeSide1Vecs[quad]);
      break;
    case 2:
      if (quadVec[quad]->isSide2nonLinear == false) continue;
      intPointsCoord = &(elNodeSide2Vecs[quad]);
      break;
    case 3:
      if (quadVec[quad]->isSide3nonLinear == false) continue;
      intPointsCoord = &(elNodeSide3Vecs[quad]);
      break;
    case 4:
      if (quadVec[quad]->isSide4nonLinear == false) continue;
      intPointsCoord = &(elNodeSide4Vecs[quad]);
      break;
    }

    const int nArcs = side % 2 ? nDivQsi[quad] - 1 : nDivEta[quad] - 1;
    //auto sidePos = [](const int side, const int &i, const int &nQsi, const int &nEta, const bool &reverse){
    if(nonLinearMapping == false) continue;
    for (int i = 0; i < nArcs; i++) {
      TPZManVector<long, 3> nodesIdVec(3, 0.);
      const int vertex1 = sidePos(side, i, nDivQsi[quad], nDivEta[quad], false);
      const int vertex2 = sidePos(side, i + 1, nDivQsi[quad], nDivEta[quad], false);


      nodesIdVec[0] = elNodeFaceVecs[quad][vertex1];
      nodesIdVec[1] = elNodeFaceVecs[quad][vertex2];
      nodesIdVec[2] = (*intPointsCoord)[i];//mid-point
      if (nodesIdVec[0] == -1 || nodesIdVec[1] == -1 || nodesIdVec[2] == -1) {
        DebugStop();
      }
      TPZGeoElRefPattern<pzgeom::TPZArc3D> *arc =
        new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, -24, *gmesh);//random material id
      arc->SetRefPattern(refpArc);
    }
  }

  return gmesh;
}

void
wgma::gmeshtools::SplitQuadMeshIntoTriangles(TPZAutoPointer<TPZGeoMesh> gmesh)
{

  ///creates refpattern
  TPZAutoPointer<TPZRefPattern> refp;
  constexpr char buf[] =
    "4 3 "
    "-50 Quad0000111111111 "
    "-1 -1 0 "
    "1 -1 0 "
    "1 1 0 "
    "-1 1 0 "
    "3 4 0  1  2  3 "
    "2 3 0  1  2 "
    "2 3 0  2  3 ";
  std::istringstream str(buf);
  refp = new TPZRefPattern(str);
  refp->GenerateSideRefPatterns();
  gRefDBase.InsertRefPattern(refp);
  if (!refp) {
    DebugStop();
  }
  auto nel = gmesh->NElements();
  TPZVec<TPZGeoEl *>sons;
  for (long iel = 0; iel < nel; iel++) {
    TPZGeoEl *gelQuad = gmesh->ElementVec()[iel];
    if(gelQuad->Type() == EQuadrilateral && gelQuad->HasSubElement() == 0){
      gelQuad->SetRefPattern(refp);
      gelQuad->Divide(sons);
    }
  }
  gmesh->BuildConnectivity();
}

using namespace wgma::gmeshtools;

QuadrilateralData::QuadrilateralData(const TPZVec<REAL> &co1,
                                     const TPZVec<REAL> &co2,
                                     const TPZVec<REAL> &co3,
                                     const TPZVec<REAL> &co4) :
  coord1(co1),coord2(co2),coord3(co3),coord4(co4)
{

  auto map_quad_side_linear = [](const TPZVec<REAL> &coord1 ,const TPZVec<REAL> &coord2,const REAL s) {
    TPZVec<REAL> point(2,0.);
    point[0] = (coord1[0] - s*coord1[0] + coord2[0] + s*coord2[0])/2.;
    point[1] = (coord1[1] - s*coord1[1] + coord2[1] + s*coord2[1])/2.;
    return point;
  };
  
  isSide1nonLinear = false;
  mapSide1 = [map_quad_side_linear,co1,co2](const REAL s){
    return map_quad_side_linear(co1,co2,s);
  };

  isSide2nonLinear = false;
  mapSide2 = [map_quad_side_linear,co2,co3](const REAL s){
    return map_quad_side_linear(co2,co3,s);
  };
  isSide3nonLinear = false;
  mapSide3 = [map_quad_side_linear,co3,co4](const REAL s){
    return map_quad_side_linear(co3,co4,s);
  };
  isSide4nonLinear = false;
  mapSide4 = [map_quad_side_linear,co4,co1](const REAL s){
    return map_quad_side_linear(co4,co1,s);
  };
}

QuadrilateralData::QuadrilateralData(const int p1, const int p2,
                                     const int p3, const int p4,
                                     const TPZVec<REAL> &co1,
                                     const TPZVec<REAL> &co2,
                                     const TPZVec<REAL> &co3,
                                     const TPZVec<REAL> &co4) :
  QuadrilateralData(co1,co2,co3,co4)
{
  pts.Resize(4);
  pts[0]=p1;
  pts[1]=p2;
  pts[2]=p3;
  pts[3]=p4;
}


void QuadrilateralData::SetMapsLinearity(const bool nl1, const bool nl2,
                                         const bool nl3, const bool nl4){
  isSide1nonLinear = nl1;
  isSide2nonLinear = nl2;
  isSide3nonLinear = nl3;
  isSide4nonLinear = nl4;
}

void QuadrilateralData::SetMapsLinearityAndAdjustEdges(
  const bool nl1,const bool nl2,const bool nl3,const bool nl4,
  const int quad, TPZVec<EdgeData> &allEdges)
{

  SetMapsLinearity(nl1,nl2,nl3,nl4);
  sides.Resize(4);
  invSide.Resize(4);
  for(int iSide = 0; iSide < 4; iSide++){
    const int &pt1 = pts[iSide];
    const int &pt2 = pts[(iSide+1)%4];
    int foundEdge = -1;
    for(int iEdge = 0; iEdge < allEdges.size(); iEdge++){
      if(allEdges[iEdge].coord1 == pt1 && allEdges[iEdge].coord2 == pt2) {
        std::cout<<"found existing edge with same orientation"<<std::endl;
        DebugStop();
      }
      if(allEdges[iEdge].coord1 == pt2 && allEdges[iEdge].coord2 == pt1) {
        foundEdge = iEdge;
        invSide = true;
        allEdges[foundEdge].quad2 = quad;
        allEdges[foundEdge].side2 = iSide + 1;
        allEdges[foundEdge].amIboundaryEdge = false;
      }
    }
    if(foundEdge == -1){
      foundEdge = allEdges.size();
      invSide[iSide] = false;
      allEdges.Resize(foundEdge + 1);
      allEdges[foundEdge].coord1 = pt1;
      allEdges[foundEdge].coord2 = pt2;
      allEdges[foundEdge].quad1 = quad;
      allEdges[foundEdge].side1 = iSide + 1;
      allEdges[foundEdge].quad2 = quad;
      allEdges[foundEdge].side2 = iSide + 1;
      allEdges[foundEdge].amIboundaryEdge = true;
      switch(iSide){
      case 0:
        allEdges[foundEdge].isNonLinearMapped = isSide1nonLinear;
        break;
      case 1:
        allEdges[foundEdge].isNonLinearMapped = isSide2nonLinear;
        break;
      case 2:
        allEdges[foundEdge].isNonLinearMapped = isSide3nonLinear;
        break;
      case 3:
        allEdges[foundEdge].isNonLinearMapped = isSide4nonLinear;
        break;
      default:
        DebugStop();
      }
    }
    sides[iSide] = foundEdge;
  }
}

void QuadrilateralData::CreateQuadrilateralRegion(const int nQsi, const int nEta)
{

  auto getPoints = [this](TPZFMatrix<REAL> &sidePts,const int nPtsQsi,const int nPtsEta,
                          const REAL iniQsi, const REAL endQsi, const REAL iniEta, const REAL endEta){

    sidePts.Resize(nPtsQsi * nPtsEta,2);
    // std::cout<<"nPtsQsi: "<<nPtsQsi <<" nPtsEta : "<< nPtsEta<<std::endl;
    // std::cout<<"iniQsi : "<< iniQsi  <<" endQsi : "<< endQsi<<std::endl;
    // std::cout<<"iniEta : "<< iniEta  <<" endEta : "<< endEta<<std::endl;

    auto mapFace = [this]
      (const REAL &qsi, const REAL &eta, REAL &x, REAL &y) {
      x=0;
      y=0;
      x += ((1.-qsi)/2)*((1.-eta)/2)*(coord1[0]);
      y += ((1.-qsi)/2)*((1.-eta)/2)*(coord1[1]);
      x += ((1.+qsi)/2)*((1.-eta)/2)*(coord2[0]);
      y += ((1.+qsi)/2)*((1.-eta)/2)*(coord2[1]);
      x += ((1.+qsi)/2)*((1.+eta)/2)*(coord3[0]);
      y += ((1.+qsi)/2)*((1.+eta)/2)*(coord3[1]);
      x += ((1.-qsi)/2)*((1.+eta)/2)*(coord4[0]);
      y += ((1.-qsi)/2)*((1.+eta)/2)*(coord4[1]);
    };

    auto s1 = [](const REAL &qsi, const REAL &eta) { return  1. * qsi;};
    auto s2 = [](const REAL &qsi, const REAL &eta) { return  1. * eta;};
    auto s3 = [](const REAL &qsi, const REAL &eta) { return -1. * qsi;};
    auto s4 = [](const REAL &qsi, const REAL &eta) { return -1. * eta;};

    auto mapLinear = [](const TPZVec<REAL>& coord1side, const TPZVec<REAL>& coord2side,const REAL &s, TPZVec<REAL> &x) {
      x[0] = coord1side[0] * (1.-s)/2 + coord2side[0] * (1+s)/2;
      x[1] = coord1side[1] * (1.-s)/2 + coord2side[1] * (1+s)/2;
    };

    auto correctionFactor1 = [](const REAL &qsi, const REAL &eta) { return 0.25 * (-1 + eta) * (-1 + eta);};
    auto correctionFactor2 = [](const REAL &qsi, const REAL &eta) { return 0.25 * ( 1 + qsi) * ( 1 + qsi);};
    auto correctionFactor3 = [](const REAL &qsi, const REAL &eta) { return 0.25 * ( 1 + eta) * ( 1 + eta);};
    auto correctionFactor4 = [](const REAL &qsi, const REAL &eta) { return 0.25 * (-1 + qsi) * (-1 + qsi);};
    int iPt = 0;
    for(int i = 0; i < nPtsEta ; i++){
      const REAL etaPt = nPtsEta > 1 ? ((REAL)(i) / (nPtsEta - 1)) : 1;
      for(int j = 0; j < nPtsQsi ; j++){
        const REAL qsiPt = nPtsQsi > 1 ? ((REAL)(j) / (nPtsQsi - 1)) : 1;
        const REAL qsi = iniQsi * (1. - qsiPt) + endQsi * qsiPt;
        const REAL eta = iniEta * (1. - etaPt) + endEta * etaPt;
        REAL xFace;
        REAL yFace;
        mapFace(qsi,eta,xFace,yFace);

        TPZVec<REAL> xSide(2,0.);

        sidePts(iPt,0) = xFace;
        sidePts(iPt,1) = yFace;
        mapLinear(coord1,coord2,s1(qsi,eta),xSide);
        sidePts(iPt,0) -= isSide1nonLinear ? correctionFactor1(qsi,eta) * (xSide[0] - mapSide1(s1(qsi,eta))[0]) : 0;
        sidePts(iPt,1) -= isSide1nonLinear ? correctionFactor1(qsi,eta) * (xSide[1] - mapSide1(s1(qsi,eta))[1]) : 0;
        mapLinear(coord2,coord3,s2(qsi,eta),xSide);
        sidePts(iPt,0) -= isSide2nonLinear ? correctionFactor2(qsi,eta) * (xSide[0] - mapSide2(s2(qsi,eta))[0]) : 0;
        sidePts(iPt,1) -= isSide2nonLinear ? correctionFactor2(qsi,eta) * (xSide[1] - mapSide2(s2(qsi,eta))[1]) : 0;
        mapLinear(coord3,coord4,s3(qsi,eta),xSide);
        sidePts(iPt,0) -= isSide3nonLinear ? correctionFactor3(qsi,eta) * (xSide[0] - mapSide3(s3(qsi,eta))[0]) : 0;
        sidePts(iPt,1) -= isSide3nonLinear ? correctionFactor3(qsi,eta) * (xSide[1] - mapSide3(s3(qsi,eta))[1]) : 0;
        mapLinear(coord4,coord1,s4(qsi,eta),xSide);
        sidePts(iPt,0) -= isSide4nonLinear ? correctionFactor4(qsi,eta) * (xSide[0] - mapSide4(s4(qsi,eta))[0]) : 0;
        sidePts(iPt,1) -= isSide4nonLinear ? correctionFactor4(qsi,eta) * (xSide[1] - mapSide4(s4(qsi,eta))[1]) : 0;

        //std::cout<<"x:   "<<sidePts(iPt,0)<<", y:   "<<sidePts(iPt,1)<<std::endl;
        iPt++;
      }
    }
  };


  //std::cout<<"--------------- intP ---------------"<<std::endl;
  getPoints(facePts,nQsi,nEta,-1., 1.,-1., 1.);
  const REAL deltaQsiInt = 1./(nQsi - 1);
  const REAL deltaEtaInt = 1./(nEta - 1);
  //std::cout<<"---------------side 1---------------"<<std::endl;
  isSide1nonLinear ?
    getPoints(side1IntPts,(nQsi-1),1,-1.+deltaQsiInt, 1.-deltaQsiInt,-1.,-1.) : (void)side1IntPts.Resize(0,0);
  //std::cout<<"---------------side 2---------------"<<std::endl;
  isSide2nonLinear ?
    getPoints(side2IntPts,1,(nEta-1), 1., 1.,-1.+deltaEtaInt, 1.-deltaEtaInt): (void)side2IntPts.Resize(0,0);
  //std::cout<<"---------------side 3---------------"<<std::endl;
  isSide3nonLinear ?
    getPoints(side3IntPts,(nQsi-1),1, 1.-deltaQsiInt,-1.+deltaQsiInt, 1., 1.): (void)side3IntPts.Resize(0,0);
  //std::cout<<"---------------side 4---------------"<<std::endl;
  isSide4nonLinear ?
    getPoints(side4IntPts,1,(nEta-1),-1.,-1., 1.-deltaEtaInt,-1.+deltaEtaInt): (void)side4IntPts.Resize(0,0);

}