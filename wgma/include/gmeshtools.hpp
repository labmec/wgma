#ifndef _GMESHTOOLS_HPP_
#define _GMESHTOOLS_HPP_

#include <pzvec.h>
#include <pzfmatrix.h>

#include <functional>

namespace wgma::gmeshtools{
  //forward declarations
  struct EdgeData;
  struct QuadrilateralData;


  TPZAutoPointer<TPZGeoMesh>
  CreateStructuredMesh(
    const TPZVec<TPZVec<REAL>> &pointsVec, TPZVec<EdgeData> &edgesVec,
    const TPZFMatrix<int> &quadPointsVec, const TPZVec<int> &matIdsQuads,
    const TPZVec<int> &nDivQsi, const TPZVec<int> &nDivEta,
    const TPZVec<bool> &side1NonLinearVec, const TPZVec<bool> &side2NonLinearVec,
    const TPZVec<bool> &side3NonLinearVec, const TPZVec<bool> &side4NonLinearVec,
    const TPZVec<TPZVec<REAL>> &thetaVec, const TPZVec<TPZVec<REAL> *> &xcRef,
    const TPZVec<int> &matIdBoundVec, const TPZVec<REAL> &boundDistVec,
    const TPZVec<REAL> &rVec, const bool nonLinearMapping = true);

  /**
     @brief Splits a mesh formed of quadrilaterals into triangles*/
  void SplitQuadMeshIntoTriangles(TPZAutoPointer<TPZGeoMesh> gmesh);
  /**
     @brief data used to represent an edge when creating structured meshes.
  */
  struct EdgeData {
    //! node id for first vertex
    int coord1{-1};
    //! node id for second vertex
    int coord2{-1};
    //! whether the edge has a non-linear mapping
    bool isNonLinearMapped{-false};
    //! whether the edge belongs to a boundary
    bool amIboundaryEdge{false};
    //! id of first neighbouring quadrilateral element
    int quad1{-1};
    /**
       @brief one-dimensional side by which the first
       quadrilateral element is a neighbour of the edge
    */
    int side1{-1};
    //! id of second neighbouring quadrilateral element
    int quad2{-1};
    /**
       @brief one-dimensional side by which the second
       quadrilateral element is a neighbour of the edge
    */
    int side2{-1};
  };

  /** 
      @brief Represents a quadrilateral region of the domain that shall be
      meshed in a structured manner
  */
  struct QuadrilateralData {
    TPZVec<int> pts;
    TPZVec<int> sides;
    TPZVec<bool> invSide;

    const TPZVec<REAL> &coord1;
    const TPZVec<REAL> &coord2;
    const TPZVec<REAL> &coord3;
    const TPZVec<REAL> &coord4;

    std::function<TPZVec<REAL>(const REAL)> mapSide1;
    std::function<TPZVec<REAL>(const REAL)> mapSide2;
    std::function<TPZVec<REAL>(const REAL)> mapSide3;
    std::function<TPZVec<REAL>(const REAL)> mapSide4;

    bool isSide1nonLinear{false};
    bool isSide2nonLinear{false};
    bool isSide3nonLinear{false};
    bool isSide4nonLinear{false};

    TPZFMatrix<REAL> side1IntPts;
    TPZFMatrix<REAL> side2IntPts;
    TPZFMatrix<REAL> side3IntPts;
    TPZFMatrix<REAL> side4IntPts;
    TPZFMatrix<REAL> facePts;

    //! constructor taking nodes coordinates
    QuadrilateralData(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2,
                      const TPZVec<REAL> &p3, const TPZVec<REAL> &p4);
    
    //! constructor taking nodes coordinates and ids
    QuadrilateralData(const int p1, const int p2,
                      const int p3, const int p4,
                      const TPZVec<REAL> &co1,
                      const TPZVec<REAL> &co2,
                      const TPZVec<REAL> &co3,
                      const TPZVec<REAL> &co4);
    //! sets mapping linearity information for each edge of the quadrilateral
    void SetMapsLinearity(const bool nl1, const bool nl2,
                          const bool nl3, const bool nl4);
    //! sets mapping linearity info and adjust edges accordingly
    void SetMapsLinearityAndAdjustEdges(const bool nl1, const bool nl2,
                                        const bool nl3, const bool nl4,
                                        const int quad, TPZVec<EdgeData> &allEdges);
    //! meshes a quadrilateral region of the domain 
    void CreateQuadrilateralRegion(const int nQsi, const int nEta);
  };
  
};

#endif /* _GMESHTOOLS_HPP_ */
