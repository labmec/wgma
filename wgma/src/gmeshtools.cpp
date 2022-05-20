#include "gmeshtools.hpp"

#include <pzgmesh.h>
#include <TPZRefPattern.h>
#include <TPZRefPatternDataBase.h>
#include <tpzarc3d.h>
#include <tpzgeoblend.h>
#include <TPZVTKGeoMesh.h>
#include <TPZGmshReader.h>

#include <fstream>

void wgma::gmeshtools::PrintGeoMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                                    std::string filename)
{
  const std::string gmeshFileNameTxt{filename+".txt"};
  std::ofstream gmeshFileTxt(gmeshFileNameTxt);
  gmesh->Print(gmeshFileTxt);
  const std::string gmeshFileNameVtk{filename+".vtk"};
  std::ofstream gmeshFileVtk(gmeshFileNameVtk);
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, gmeshFileVtk, true);
}

/****
     The following are auxiliary functions for creating
     the structured mesh and are defined below
****/
std::map<std::pair<int,int>,int>
CreateEdgeMap(const TPZVec<wgma::gmeshtools::EdgeData> &edgesVec,
              const TPZVec<TPZVec<REAL>> &pointsVec);

TPZVec<TPZManVector<std::pair<int,bool>,4>>
IdentifyQuadEdges(const TPZVec<wgma::gmeshtools::QuadData> &quadsVec,
                  const std::map<std::pair<int,int>,int> &edge_map);

void
CheckNElEdges(const TPZVec<TPZManVector<std::pair<int,bool>,4>> &quadEdges,
               const TPZVec<wgma::gmeshtools::EdgeData> &edgesVec);

//linear mapping of an edge from coord1 to coord2
TPZVec<REAL> LinMapEdge(const TPZVec<REAL> &coord1, const TPZVec<REAL> &coord2,
                        const REAL s);
//linear mapping of a face with nodes {coord1,coord2,coord3,coord4}
TPZVec<REAL> LinMapFace(const TPZVec<REAL> &coord1, const TPZVec<REAL> &coord2,
                        const TPZVec<REAL> &coord3, const TPZVec<REAL> &coord4,
                        const REAL xi, const REAL eta);

/**Insert in the geometric mesh all points along the edges*/
void CreateEdgePts(const TPZVec<TPZVec<REAL>> &pointsVec,
                   const TPZVec<wgma::gmeshtools::EdgeData> &edgesVec,
                   TPZAutoPointer<TPZGeoMesh> gmesh,
                   TPZVec<TPZVec<int>> &edgePts,
                   int &nodeid);

/**Identify all nodes for all quads.
   Interior pts will be inserted in the geometric mesh, boundary points are located
   on the edges*/
void CreateQuadPts(const TPZVec<wgma::gmeshtools::QuadData> &quadsVec,
                   const TPZVec<TPZManVector<std::pair<int,bool>,4>> &quadEdges,
                   const TPZVec<wgma::gmeshtools::EdgeData> &edgesVec,
                   const TPZVec<TPZVec<REAL>>& pointsVec,
                   const TPZVec<TPZVec<int>> &edgePts,
                   TPZAutoPointer<TPZGeoMesh> gmesh,
                   TPZVec<TPZVec<int>> &quadPts,
                   int &nodeid);

TPZAutoPointer<TPZGeoMesh>
wgma::gmeshtools::CreateStructuredMesh(
  const TPZVec<TPZVec<REAL>> &pointsVec,
  const TPZVec<QuadData> &quadsVec,
  const TPZVec<EdgeData> &edgesVec,
  const bool nonLinearMapping)
{

  // given nodes <n1,n2> it returns the position of the edge in edgesVec
  std::map<std::pair<int,int>,int> edge_map
    = CreateEdgeMap(edgesVec,pointsVec);

  /**
     Checks which edges compose a given quuadrilateral and its orientation
     (true if it is the same, false if it is inverted)*/
  TPZVec<TPZManVector<std::pair<int,bool>,4>>
    quadEdges = IdentifyQuadEdges(quadsVec,edge_map);

  //Adjust nel for opposite edges of any quadrilateral
  CheckNElEdges(quadEdges, edgesVec);


  TPZAutoPointer<TPZGeoMesh> gmesh{new TPZGeoMesh()};
  
  const auto nEdges = edgesVec.size();
  ///Stores edge points for all edges
  TPZVec<TPZVec<int>> edgePts(nEdges,{});

  int nodeid = 0;
  CreateEdgePts(pointsVec,edgesVec,gmesh,edgePts,nodeid);

  const auto nQuads = quadsVec.size();
  /**Stores all points for all quads**/
  TPZVec<TPZVec<int>> quadPts(nQuads,{});
  CreateQuadPts(quadsVec,quadEdges,edgesVec,pointsVec,
                edgePts,gmesh,quadPts,nodeid);
  
  
  //Creating volumetric elements for each quadrilateral
  for(int iq = 0; iq < nQuads; iq++){
    const auto &quad = quadsVec[iq];
    const auto &qpts = quadPts[iq];
    const auto matid = quad.m_matid;
    const auto eltype = quad.m_type;

    const auto &qedges = quadEdges[iq];
    
    const TPZVec<bool> sign_edges=
      {qedges[0].second, qedges[1].second,qedges[2].second, qedges[3].second};
    
    const TPZVec<int> id_edges =
      {qedges[0].first, qedges[1].first,qedges[2].first, qedges[3].first};

    const auto nelx = edgesVec[id_edges[0]].m_nel;
    const auto nely = edgesVec[id_edges[1]].m_nel;

    const int nelnodes = eltype == ElType::Tri ? 3 : 4;

    TPZManVector<int64_t,4> nodeVec(nelnodes,0);
    for(int iely = 0; iely < nely; iely++){
      /**
         if the element is touching bottom or 
         top edges and the edge is non linear,
         a blend el should be created
      */
      const bool blendy =
        (iely == 0 && edgesVec[id_edges[0]].m_map)||
        (iely == nely-1 && edgesVec[id_edges[2]].m_map);
      
      
      for(int ielx = 0; ielx < nelx; ielx++){
        const bool blendx =
        (ielx == 0 && edgesVec[id_edges[3]].m_map)||
        (ielx == nelx-1 && edgesVec[id_edges[1]].m_map);
        
        
        const int node0 = iely * (nelx+1) + ielx;//Lower left vertex
        const int node1 = iely * (nelx+1) + ielx + 1;//Lower right vertex
        const int node2 = (iely + 1) * (nelx+1) + ielx + 1;//Upper right vertex
        const int node3 = (iely + 1) * (nelx+1) + ielx;//Upper left vertex
        
        if ((blendx || blendy) && nonLinearMapping){
          if(eltype == ElType::Tri){
            nodeVec = {qpts[node0],qpts[node1],qpts[node2]};
            new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(nodeVec, matid,*gmesh);
            nodeVec = {qpts[node2],qpts[node3],qpts[node0]};
            new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(nodeVec, matid,*gmesh);
          }else{
            nodeVec = {qpts[node0],qpts[node1],qpts[node2], qpts[node3]};
            new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>>(nodeVec, matid,*gmesh);
          }
        }else{
          if(eltype == ElType::Tri){
            nodeVec = {qpts[node0],qpts[node1],qpts[node2]};
            new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeVec, matid,*gmesh);
            nodeVec = {qpts[node2],qpts[node3],qpts[node0]};
            new TPZGeoElRefPattern<pzgeom::TPZGeoTriangle>(nodeVec, matid,*gmesh);
          }else{
            nodeVec = {qpts[node0],qpts[node1],qpts[node2], qpts[node3]};
            new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodeVec, matid,*gmesh);
          }
        }
      }
    } 
    
  }

  for(int ie = 0; ie < nEdges; ie++){
    const auto edge = edgesVec[ie];
    if (edge.m_create_el || (edge.m_map && nonLinearMapping)){//either non-linear or bc
      const auto matid = edge.m_matid;
      const auto nel = edge.m_nel;
      //non-linear elements have three nodes, linear elements just two
      const auto nelnodes = edge.m_map ? 3 : 2;

      TPZVec<int64_t> midnodes;
      if(edge.m_map){
        midnodes.Resize(nel);
        const auto nmeshnodes = gmesh->NodeVec().NElements();
        gmesh->NodeVec().Resize(nmeshnodes+nel);
        //size of each element in the parametric domain of the edge
        const REAL el_param_side = 2. / nel;
        for(int iel = 0; iel < nel; iel++){
          //center of each element on the parametric domain
          const auto s = -1 + el_param_side * iel + el_param_side/2;
          auto pt = edge.m_map(s);
          gmesh->NodeVec()[nodeid].Initialize(pt,*gmesh);
          midnodes[iel] = nodeid;
          nodeid++;
        }
      }
      
      TPZManVector<long, 3> nodesIdVec(nelnodes, 0.);
      for(int iel = 0; iel < nel; iel++){
        const int vertex1 = edgePts[ie][iel];
        const int vertex2 = edgePts[ie][iel+1];

        nodesIdVec[0] = vertex1;
        nodesIdVec[1] = vertex2;
        if(edge.m_map){
          nodesIdVec[2] = midnodes[iel];
          new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matid, *gmesh);
        }else{
          new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesIdVec, matid, *gmesh);
        }
      }
    }
  }
  return gmesh;
  
}

std::map<std::pair<int,int>,int>
CreateEdgeMap(const TPZVec<wgma::gmeshtools::EdgeData> &edgesVec,
              const TPZVec<TPZVec<REAL>> &pointsVec){

  //checks if two points are the same
  auto CheckPt = [](const TPZVec<REAL> &coord1,
                    const TPZVec<REAL> &coord2){
    constexpr REAL tol = std::numeric_limits<REAL>::epsilon();
    auto diff = 0;
    for(int i = 0; i < 3; i++){
      diff += (coord2[i] - coord1[i])*(coord2[i] - coord1[i]);
    }
    if(diff < tol) return true;
    else return false;
  };
  const auto nEdges = edgesVec.size();
  std::map<std::pair<int,int>,int> edge_map;
  /*
    here we will fill edge_map and check for mapping consistency
   */
  for (auto i = 0; i < nEdges; i++){
    auto &edge = edgesVec[i];
    std::pair<int,int> nodes(edge.m_nodes[0],edge.m_nodes[1]);
    std::pair<int,int> nodes_alt(edge.m_nodes[1],edge.m_nodes[0]);
    if(edge_map.find(nodes) == edge_map.end() &&
       edge_map.find(nodes_alt) == edge_map.end()){
      edge_map.insert({nodes,i});
    }else{
      PZError<<__PRETTY_FUNCTION__
             <<"\nFound duplicate edge with nodes"
             <<"\t"<<nodes.first<<" "<<nodes.second
             <<"\nAborting...\n";
      DebugStop();
    }
    if(edge.m_map){
      const auto &coord1 = pointsVec[edge.m_nodes[0]];
      const auto &coord2 = pointsVec[edge.m_nodes[1]];
      const auto map1 = edge.m_map(-1);
      const auto map2 = edge.m_map( 1);
      if(!CheckPt(map1,coord1) ||
         !CheckPt(map2,coord2)){
        PZError<<__PRETTY_FUNCTION__
               <<"\nEdge "<<i<<" with inconsistent mapping"
               <<"\nAborting...\n";
        DebugStop();
      }
    }
  }
  
  return std::move(edge_map);
}

TPZVec<TPZManVector<std::pair<int,bool>,4>>
IdentifyQuadEdges(const TPZVec<wgma::gmeshtools::QuadData> &quadsVec,
                  const std::map<std::pair<int,int>,int> &edge_map){
  const auto nQuads = quadsVec.size();

  /**
     Checks which edges compose a given quuadrilateral and its orientation
     (true if it is the same, false if it is inverted)*/
  TPZVec<TPZManVector<std::pair<int,bool>,4>>
    quadEdges(nQuads,
              {{-1,false},{-1,false},{-1,false},{-1,false}});

  for(int iq = 0; iq < nQuads; iq++){
    const auto &quad = quadsVec[iq];
    for(int ie = 0; ie < 4; ie++){
      const auto &n1 = quad.m_nodes[ie];
      const auto &n2 = quad.m_nodes[(ie+1)%4];
      if(edge_map.find({n1,n2}) != edge_map.end()){
        //same orientation
        quadEdges[iq][ie] = {edge_map.at({n1,n2}),true};
      }else if(edge_map.find({n2,n1}) != edge_map.end()){
        //other orientation
        quadEdges[iq][ie] = {edge_map.at({n2,n1}),false};
      }else{
        PZError<<__PRETTY_FUNCTION__
               <<"\n Edge "<<ie<<" not found for quad "<<iq
               <<"\nnodes: "<<n1<<" "<<n2<<"\nAborting...\n";
        DebugStop();
        //edge not found
      }
    }
    
  }
  return std::move(quadEdges);
}

void
CheckNElEdges(const TPZVec<TPZManVector<std::pair<int,bool>,4>> &quadEdges,
              const TPZVec<wgma::gmeshtools::EdgeData> &edgesVec)
{
  const int nQuads = quadEdges.size();
  bool check = false;
  while(!check){
    check = true;
    for(int iq = 0; iq < nQuads; iq++){
      const auto e0 = quadEdges[iq][0].first;
      const auto e1 = quadEdges[iq][1].first;
      const auto e2 = quadEdges[iq][2].first;
      const auto e3 = quadEdges[iq][3].first;
      if(edgesVec[e0].m_nel != edgesVec[e2].m_nel){
        check = false;
        std::cout<<__PRETTY_FUNCTION__
                 <<"\nEdges "<<e0<<" and "<<e3
                 <<"\nhave incompatible nels in\n"
                 <<"quad "<<iq<<". Adjusting...";
        edgesVec[e0].m_nel = edgesVec[e2].m_nel =
          std::max(edgesVec[e0].m_nel,edgesVec[e2].m_nel);
      }
      if(edgesVec[e1].m_nel != edgesVec[e3].m_nel){
        check = false;
        std::cout<<__PRETTY_FUNCTION__
                 <<"\nEdges "<<e1<<" and "<<e3
                 <<"\nhave incompatible ndivs in\n"
                 <<"quad "<<iq<<". Adjusting...";
        edgesVec[e1].m_nel = edgesVec[e3].m_nel =
          std::max(edgesVec[e1].m_nel,edgesVec[e3].m_nel);
      }
    }
  }
}

TPZVec<REAL> LinMapEdge(const TPZVec<REAL> &coord1, const TPZVec<REAL> &coord2,
                        const REAL s)
{
  TPZVec<REAL> point(3,0.);
  point[0] = (coord1[0] - s*coord1[0] + coord2[0] + s*coord2[0])/2.;
  point[1] = (coord1[1] - s*coord1[1] + coord2[1] + s*coord2[1])/2.;
  return point;
}


TPZVec<REAL> LinMapFace(const TPZVec<REAL> &coord1, const TPZVec<REAL> &coord2,
                        const TPZVec<REAL> &coord3, const TPZVec<REAL> &coord4,
                        const REAL xi, const REAL eta)
{

  TPZVec<REAL> x(3,0.);

  for(int i = 0; i < 3; i++){
    x[i] += ((1.-xi)/2)*((1.-eta)/2)*(coord1[i]);
  }
  for(int i = 0; i < 3; i++){
    x[i] += ((1.+xi)/2)*((1.-eta)/2)*(coord2[i]);
  }
  for(int i = 0; i < 3; i++){
    x[i] += ((1.+xi)/2)*((1.+eta)/2)*(coord3[i]);
  }
  for(int i = 0; i < 3; i++){
    x[i] += ((1.-xi)/2)*((1.+eta)/2)*(coord4[i]);
  }
  return x;
}


void CreateEdgePts(const TPZVec<TPZVec<REAL>> &pointsVec,
                   const TPZVec<wgma::gmeshtools::EdgeData> &edgesVec,
                   TPZAutoPointer<TPZGeoMesh> gmesh,
                   TPZVec<TPZVec<int>> &edgePts,
                   int &nodeid)
{
  const auto nEdges = edgePts.size();
  TPZVec<int> nodeIds(pointsVec.size(),-1);
  for(int ie = 0; ie < nEdges; ie++){
    const auto nel = edgesVec[ie].m_nel;
    const auto n1 = edgesVec[ie].m_nodes[0];
    const auto n2 = edgesVec[ie].m_nodes[1];
    auto edgemap = edgesVec[ie].m_map;
    if(!edgemap){
      const auto c1 = pointsVec[n1];
      const auto c2 = pointsVec[n2];
      edgemap = [c1,c2](const REAL s){
        return LinMapEdge(c1,c2,s);
      };
    }
    auto &epts = edgePts[ie];
    const auto npts = nel+1;
    epts.Resize(npts);

    //check if first point was inserted in mesh already
    int ipt = 0;
    if(nodeIds[n1] != -1){
      epts[ipt] = nodeIds[n1];
    }else{
      auto nel_mesh = gmesh->NodeVec().NElements();
      gmesh->NodeVec().Resize(nel_mesh+1);
      auto c1 = pointsVec[n1];
      gmesh->NodeVec()[nodeid].Initialize(c1,*gmesh);
      epts[ipt] = nodeid;
      nodeIds[n1] = nodeid;
      nodeid++;
    }
    //let us insert all mid points
    const auto nel_mesh = gmesh->NodeVec().NElements();
    gmesh->NodeVec().Resize(nel_mesh+nel-1);
    
    for(ipt=1; ipt < npts-1; ipt++){
      const REAL s = -1 + ((REAL)ipt)/(npts-1) * 2;
      auto node = edgemap(s);
      gmesh->NodeVec()[nodeid].Initialize(node,*gmesh);
      epts[ipt] = nodeid;
      nodeid++;
    }

    //check if last point was inserted in mesh already
    if(nodeIds[n2] != -1){
      epts[ipt] = nodeIds[n2];
    }else{
      auto nel_mesh = gmesh->NodeVec().NElements();
      gmesh->NodeVec().Resize(nel_mesh+1);
      auto c2 = pointsVec[n2];
      gmesh->NodeVec()[nodeid].Initialize(c2,*gmesh);
      epts[ipt] = nodeid;
      nodeIds[n2] = nodeid;
      nodeid++;
    }
  }
}

void CreateQuadPts(const TPZVec<wgma::gmeshtools::QuadData> &quadsVec,
                   const TPZVec<TPZManVector<std::pair<int,bool>,4>> &quadEdges,
                   const TPZVec<wgma::gmeshtools::EdgeData> &edgesVec,
                   const TPZVec<TPZVec<REAL>>& pointsVec,
                   const TPZVec<TPZVec<int>> &edgePts,
                   TPZAutoPointer<TPZGeoMesh> gmesh,
                   TPZVec<TPZVec<int>> &quadPts,
                   int &nodeid)
{

  typedef std::function<REAL (const REAL&, const REAL)> projfunc;

  //project interior point to a given edge
  projfunc side_proj[] = {
    [](const REAL &xi, const REAL &eta) { return  1. * xi;},
    [](const REAL &xi, const REAL &eta) { return  1. * eta;},
    [](const REAL &xi, const REAL &eta) { return -1. * xi;},
    [](const REAL &xi, const REAL &eta) { return -1. * eta;}
  };
  /**correction factor to aply for each non linear edge 
     to account for the deviation from the linear mapping*/
  projfunc corr_fact[] = {
    [](const REAL &xi, const REAL &eta) { return 0.25 * (-1 + eta) * (-1 + eta);},
    [](const REAL &xi, const REAL &eta) { return 0.25 * ( 1 + xi) * ( 1 + xi);},
    [](const REAL &xi, const REAL &eta) { return 0.25 * ( 1 + eta) * ( 1 + eta);},
    [](const REAL &xi, const REAL &eta) { return 0.25 * (-1 + xi) * (-1 + xi);}
  };
  
  const int nQuads = quadsVec.size();

  for(int iq = 0; iq < nQuads; iq++){
    const auto &quad = quadsVec[iq];
    const auto &qedges = quadEdges[iq];
    TPZVec<TPZVec<REAL>> coord_vec = {
      pointsVec[quad.m_nodes[0]],
      pointsVec[quad.m_nodes[1]],
      pointsVec[quad.m_nodes[2]],
      pointsVec[quad.m_nodes[3]]
      
    };
    const TPZVec<bool> sign_edges=
      {qedges[0].second, qedges[1].second,qedges[2].second, qedges[3].second};
    const TPZVec<int> id_edges =
      {qedges[0].first, qedges[1].first,qedges[2].first, qedges[3].first};

    //we can safely assume that ndiv for opposite edges are the same
    const auto nelx = edgesVec[id_edges[0]].m_nel;
    const auto nely = edgesVec[id_edges[1]].m_nel;
    const auto npts = (nelx+1)*(nely+1);
    /*since the edge points have already been created, we must create
      (nelx-1)*(nely-1) interior points * 
     */
    const auto n_newpts = (nelx-1)*(nely-1);

    const auto nmeshnodes = gmesh->NodeVec().NElements();
    gmesh->NodeVec().Resize(nmeshnodes+n_newpts);
    
    auto &qpts = quadPts[iq];
    qpts.Resize(npts,-1);

    auto ipt = 0;
    
    //first, let us assign the points along edge 0
    for(int ix = 0; ix < nelx+1; ix++){
      const auto &pts_edge = edgePts[id_edges[0]];
      const int pt_edge = sign_edges[0] ? ix : nelx - ix;
      qpts[ipt] = pts_edge[pt_edge];
      ipt++;
    }
    //now, some interior points and some points along edges 1 and 3
    for(int iy=1; iy < nely; iy++){
      const auto eta = -1 * (REAL)(nely-iy)/nely + 1 * (REAL)iy/nely;
      for (auto ix = 0; ix < nelx+1; ix++){
        if(ix == 0){//left edge, edge 3
          //edge 3 should be inverted
          const auto &pts_edge = edgePts[id_edges[3]];
          const int pt_edge = sign_edges[3] ? nely - iy : iy;
          qpts[ipt] = pts_edge[pt_edge];
          ipt++;
          continue;
        }else if (ix == nelx){//right edge, edge 1
          const int pt_edge = sign_edges[1] ? iy : nely - iy;
          const auto &pts_edge = edgePts[id_edges[1]];
          qpts[ipt] = pts_edge[pt_edge];
          ipt++;
          continue;
        }
        const auto xi = -1 * (REAL)(nelx-ix)/nelx + 1 * (REAL)ix/nelx;

        auto pt = LinMapFace(coord_vec[0], coord_vec[1],
                                coord_vec[2], coord_vec[3],
                                xi,eta);

        for(auto ie = 0; ie < 4; ie++){
          const auto edge_id = id_edges[ie];
          if(edgesVec[edge_id].m_map){//non linear map
            //check if edge orientation is inverted
            const int edge_sign = sign_edges[ie] ? 1 : -1;
            const auto n1 = sign_edges[ie] ? ie : (ie+1) % 4;
            const auto n2 = sign_edges[ie] ? (ie+1) % 4 : ie;
            const TPZVec<REAL> & v1 = coord_vec[n1];
            const TPZVec<REAL> & v2 = coord_vec[n2];
            //parametric coordinate of the edge
            const auto s_edge= edge_sign * side_proj[ie](xi,eta);
            //point if the edge was mapped linearly
            const auto lin_xedge = LinMapEdge(v1, v2, s_edge);
            const auto blend = corr_fact[ie](xi,eta);
            //actual point of the edge
            const auto nlin_xedge = edgesVec[edge_id].m_map(s_edge);
            for(auto x = 0 ; x < 3; x++){
              pt[x] -= blend*(lin_xedge[x]-nlin_xedge[x]);
            }
            
          }
        }
        gmesh->NodeVec()[nodeid].Initialize(pt,*gmesh);
        qpts[ipt] = nodeid;
        nodeid++;
        ipt++;
      }
    }

    //finally, the last points are along edge 2
    for(int ix = 0; ix < nelx+1; ix++){
      const auto &pts_edge = edgePts[id_edges[2]];
      //edge 2 has a different orientation
      const int pt_edge = sign_edges[2] ? nelx - ix : ix;
      qpts[ipt] = pts_edge[pt_edge];
      ipt++;
    }

    if(ipt != npts){
      PZError<<__PRETTY_FUNCTION__
             <<"\nERROR: wrong number of interior points.\n"
             <<"nelx: "<<nelx<<" nely: "<<nely<<"\n"
             <<"npts(read): "<<ipt<<" npts(calc): "<<npts
             <<" n_newpts: "<<n_newpts<<"\n"
             <<"Aborting...\n";
      DebugStop();
    }
  }
}


TPZAutoPointer<TPZGeoMesh>
wgma::gmeshtools::ReadGmshMesh(const std::string filename,
                               const REAL scale,
                               TPZVec<std::map<std::string,int>> & matids,
                               const bool verbose)
{
  std::map<int64_t,int64_t> dummy;
  return ReadPeriodicGmshMesh(filename,scale,matids,dummy, verbose);
}

TPZAutoPointer<TPZGeoMesh>
wgma::gmeshtools::ReadPeriodicGmshMesh(const std::string filename,
                                       const REAL scale,
                                       TPZVec<std::map<std::string,int>> & matids,
                                       std::map<int64_t,int64_t> &periodic_bcs,
                                       const bool verbose)
{
  
  TPZGmshReader meshReader;
  meshReader.SetCharacteristiclength(scale);
  TPZAutoPointer<TPZGeoMesh> gmesh =
    meshReader.GeometricGmshMesh(filename);

  matids = meshReader.GetDimNamePhysical();
  if(verbose){
    for(int i = 0; i < 3; i++){
      std::cout<<"materials with dim "<<i<<std::endl;
      for(auto &mat : matids[i]){
        std::cout<<"\t name: "<<mat.first <<" id: "<<mat.second<<std::endl;
      }
    }
  }
  
  periodic_bcs = meshReader.GetPeriodicEls();
  return gmesh;
}

std::optional<int>
wgma::gmeshtools::FindPMLNeighbourMaterial(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const int pmlId,
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
    return std::nullopt;
  }
  return closestEl->MaterialId();
}


std::optional<int>
wgma::gmeshtools::FindBCNeighbourMat(TPZAutoPointer<TPZGeoMesh> gmesh,
                                     const int mat,
                                     const std::set<int> &volmats)
{

  for(auto gel : gmesh->ElementVec()){
    if(!gel || gel->MaterialId() != mat){continue;}

    auto gelside = TPZGeoElSide(gel,gel->NSides()-1);
    auto neigh = gelside.Neighbour();
    while(neigh != gelside){
      
      if(auto el = neigh.Element(); el){
        const auto matid = el->MaterialId();
        if(volmats.count(matid)){//matid is in the set of valid mats
          return matid;
        }
      }
      neigh = neigh.Neighbour();
    }
  }
  
  return std::nullopt;
}

void wgma::gmeshtools::SetExactArcRepresentation(TPZAutoPointer<TPZGeoMesh> gmesh,
                                                 const TPZVec<ArcData> &circles)
{

  //creates mid node for arc elements
  auto CreateMidNode = [](const TPZVec<REAL> &x1, const TPZVec<REAL> &x2,
                          const REAL r,
                          const REAL xc, const REAL yc, const REAL zc){
    TPZVec<REAL> x3(3,0);

    //first we get its distance from xc
    x3[0] = (x1[0] + x2[0])/2 - xc;
    x3[1] = (x1[1] + x2[1])/2 - yc;
    x3[2] = (x1[2] + x2[2])/2 - zc;

    //norm of the vector
    const auto vecnorm = sqrt(x3[0]*x3[0] + x3[1]*x3[1] + x3[2]*x3[2]);

    //mid-arc coordinates
    x3[0] = xc + r * x3[0]/vecnorm;
    x3[1] = yc + r * x3[1]/vecnorm;
    x3[2] = zc + r * x3[2]/vecnorm;
    return x3;
  };

  
  //stores the first neighbour for each side before deleting element
  auto StoreNeighboursAndDeleteEl = [gmesh](TPZGeoEl *gel, bool skip_last_side){
    //skip last side: for 1D elements in 2D meshes, we want ALL the neighbours
    //for 2D elements in 2D meshes, we skip the last side
    const int diff = skip_last_side ? 1 : 0;
    auto nsides = gel->NSides() - diff;
    TPZVec<TPZGeoElSide> neighs(nsides);
    for(int i = 0; i < nsides; i++){
      neighs[i] = gel->Neighbour(i);
    }
    //now we delete the element and its connectivities
    gmesh->DeleteElement(gel);
    gel = nullptr;//just to be on the safeside
    return neighs;
  };
  
  std::map<int,int> arc_ids;
  std::map<int,bool> found_arcs;
  for(auto i = 0; i < circles.size(); i++){
    arc_ids[circles[i].m_matid] = i;
    found_arcs[circles[i].m_matid] = false;
  }

  
  for(auto el : gmesh->ElementVec()){
    //this way we avoid processing recently inserted elements
    if(!el || el->IsLinearMapping() == false) continue;
    const int matid = el->MaterialId();
    const bool is_arc = arc_ids.find(matid) != arc_ids.end();

    if(is_arc){//found arc
      found_arcs[matid] = true;
      const int arc_pos = arc_ids[matid];
      const REAL r = circles[arc_pos].m_radius;
      const REAL xc = circles[arc_pos].m_xc;
      const REAL yc = circles[arc_pos].m_yc;
      const REAL zc = circles[arc_pos].m_zc;
      
      
      TPZManVector<REAL,3> x1(3,0), x2(3,0), x3(3,0);
        
      el->Node(0).GetCoordinates(x1);
      el->Node(1).GetCoordinates(x2);

      x3 = CreateMidNode(x1, x2, r, xc, yc, zc);
        
      const auto nodeid = gmesh->NodeVec().AllocateNewElement();
      gmesh->NodeVec()[nodeid].Initialize(x3,*gmesh);
      TPZManVector<int64_t,3> nodesIdVec{el->Node(0).Id(), el->Node(1).Id(), nodeid};

      //we store the first neighbour of each side of the element to be deleted
      constexpr bool skip_side_arc{false};
      auto neighs = StoreNeighboursAndDeleteEl(el, skip_side_arc);
        
      //insert new element on the mesh
      auto arc = new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matid, *gmesh);
      //insert connectivities
      for(int i = 0; i< neighs.size();i++){
        TPZGeoElSide gelside(arc,i);
        neighs[i].SetConnectivity(gelside);
      }
        

      /**now we need to replace each neighbour of a blend element,
         so it can be deformed accordingly*/
      TPZGeoElSide gelside(arc,arc->NSides()-1);
      //now we iterate through all the neighbours of the linear side
      TPZGeoElSide neighbour = gelside.Neighbour();
      //let us store all the neighbours
      std::set<TPZGeoElSide> all_neighs;
      while(neighbour.Exists() && neighbour != gelside){
        all_neighs.insert(neighbour);
        neighbour = neighbour.Neighbour();
      }
      //let us replace all the neighbours
      for(auto neighbour : all_neighs){ 
        const auto neigh_side = neighbour.Side();
        auto neigh_el = neighbour.Element();

        /*let us take into account the possibility that
          one triangle might be neighbour of two arcs
         */
        TPZGeoEl *new_neigh{nullptr};
        if(!neigh_el->IsGeoBlendEl()){
          const auto neigh_matid = neigh_el->MaterialId();
          const auto neigh_nnodes = neigh_el->NNodes();
          TPZManVector<int64_t,4> neigh_nodes(neigh_nnodes,-1);
          //lets copy all the nodes
          for(int in = 0; in < neigh_nnodes; in++){
            neigh_nodes[in] = neigh_el->Node(in).Id();
          }
          
          const auto neigh_type = neigh_el->Type();

          constexpr bool skip_side_neigh{true};
          auto neighs = StoreNeighboursAndDeleteEl(neigh_el, skip_side_neigh);
          
          //create new element
          if(neigh_type == MElementType::ETriangle){
            new_neigh = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<
              pzgeom::TPZGeoTriangle>>(neigh_nodes, neigh_matid,*gmesh);
          }else if(neigh_type == MElementType::EQuadrilateral){
            new_neigh = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<
              pzgeom::TPZGeoQuad>>(neigh_nodes, neigh_matid,*gmesh);
          }
          //insert connectivities
          for(int i = 0; i<neighs.size();i++){
            TPZGeoElSide mygelside(new_neigh,i);
            neighs[i].SetConnectivity(mygelside);
          }
          //insert last remaining connectivity
          {
            TPZGeoElSide mygelside(new_neigh, new_neigh->NSides()-1);
            mygelside.SetConnectivity(mygelside);
          }
        }else{
          new_neigh = neigh_el;
        }
        new_neigh->SetNeighbourForBlending(neigh_side);
        new_neigh->BuildBlendConnectivity();
      }
    }
  }

  for(auto arc : found_arcs){
    if(!arc.second){
      PZError<<__PRETTY_FUNCTION__
             <<"\n arc "<<arc.first<<" not found in mesh"<<std::endl;
    }
  }
}

void wgma::gmeshtools::DirectionalRefinement(TPZAutoPointer<TPZGeoMesh> gmesh,
                           std::set<int> matids, const int nrefs)
{
  /*
    We initialise now the database of refinement patterns
  */

  const auto dim = gmesh->Dimension();
  for(auto ix = 1; ix <=dim; ix++){
    gRefDBase.InitializeRefPatterns(ix);
  }
    
  for(int iref = 0; iref < nrefs; iref++){
    const int nels = gmesh->NElements();
    for(int el = 0; el < nels; el++){
      auto *gel = gmesh->Element(el);
      if(gel && gel->NSubElements() == 0){
        TPZRefPatternTools::RefineDirectional(gel, matids);
      }
    }
  }
}