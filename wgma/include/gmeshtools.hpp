#ifndef _GMESHTOOLS_HPP_
#define _GMESHTOOLS_HPP_

#include <pzvec.h>
#include <pzfmatrix.h>

#include <functional>

namespace wgma::gmeshtools{

  enum class ElType {Tri, Quad};
  //! Prints a geometric mesh in .txt and .vtk format
  void PrintGeoMesh(TPZAutoPointer<TPZGeoMesh> gmesh, std::string filename);

  /** @brief Parametric map of a circunference arc

      The arc goes from theta[0] to theta[1], over a circunference
      centered in xc with radius r.*/
  inline TPZVec<REAL> CircArcMap(const TPZVec<REAL> &theta ,
                                 const TPZVec<REAL> &xc,
                                 const REAL r, const REAL s)
  {
    TPZVec<REAL> point(3,0.);
    point[0] = xc[0] + r*cos((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
    point[1] = xc[1] + r*sin((theta[1] + s*theta[1] + theta[0] - s*theta[0])/2.);
    return point;
  };
  
  /**
   @brief Used to pass information about a quadrilateral region to
  wgma::gmeshtools::CreateStructuredMesh*/
  struct QuadData{
    /** 
        @brief Ids of corner nodes
        It should be defined such that m_nodes[i] and m_nodes[i+1] are
        vertices of an edge.
    */
    TPZManVector<int,4> m_nodes = {-1,-1,-1,-1};
    //! Material identifier of the quadrilateral region
    int m_matid{-999};

    ElType m_type{ElType::Tri};
  };
  /**
   @brief Used to pass information about an edge of a quadrilateral region to
   wgma::gmeshtools::CreateStructuredMesh */
  struct EdgeData{
    //! Ids of corner nodes
    TPZManVector<int,2> m_nodes = {-1,-1};
    /** @brief Mapping of the edge (should be left unchanged if edge is linear)
     It should be defined such that m_map(-1) = m_nodes[0] and m_map(1) = m_nodes[1]*/
    std::function<TPZVec<REAL>(const REAL)> m_map{nullptr};
    //! Number of elements of the edge
    int m_nel{-1};
    //! Material identifier of the edge (only relevant if edge belongs to boundary)
    int m_matid{-999};
    /** @brief Whether 1D elements should be created along the edge.
        It should be true for boundary elements, and elements with
        non linear mappings will be created regardless of this value.
    */
    bool m_create_el{false};
  };
  /**
     @brief Creates a structured geometrical mesh. 
     The parameters describe (possibly curved) quadrilaterals that will be then meshed.
     The parameters related to non-linear sides refer to the parametric representation
     of a circunference arc.
     @param[in] pointsVec array with coordinates of all quadrilateral vertices.
     @param[in] quadsVec array with quadrilateral data
     @param[in] edgesVec array with quadrilaterals' edges data
     @param[in] nonLinearMapping whether to actually generate curved elements
   */
  TPZAutoPointer<TPZGeoMesh>
  CreateStructuredMesh(
    const TPZVec<TPZVec<REAL>> &pointsVec,
    const TPZVec<QuadData> &quadsVec,
    const TPZVec<EdgeData> &edgesVec,
    const bool nonLinearMapping);

  /**
     @brief Reads a .msh mesh file, creates a TPZGeoMesh and
     stores the material identifiers (named physical regions) that were found.
     @param[in] filename name of the .msh file to be read
     @param[in] scale factor (characteristic length) to be applied in the mesh
     @param[out] matids on exit, it stores the pairs (name, matid) indexed by dimension
     @param[in] verbose whether to output information about read materials

     @note position [x] of matids will contain the materials of dimension x.
   */
  TPZAutoPointer<TPZGeoMesh>
  ReadGmshMesh(const std::string filename,
               const REAL scale,
               TPZVec<std::map<std::string,int>> & matids,
               bool verbose = true);
  /**
     @brief Reads a .msh mesh file, creates a TPZGeoMesh and
     stores the material identifiers (named physical regions) that were found.
     Also, returns periodicity info (corresponding elements' ids).
     @param[in] filename name of the .msh file to be read
     @param[in] scale factor (characteristic length) to be applied in the mesh
     @param[out] matids on exit, it stores the pairs (name, matid) indexed by dimension
     @param[out] periodic_els corresponding indexes of periodic elements.
     @param[in] verbose whether to output information about read materials
     @note position [x] of matids will contain the materials of dimension x.
   */
  TPZAutoPointer<TPZGeoMesh>
  ReadPeriodicGmshMesh(const std::string filename,
                       const REAL scale,
                       TPZVec<std::map<std::string,int>> & matids,
                       std::map<int64_t,int64_t> &periodic_els,
                       bool verbose = true);

  //! Stores data for allowing exact geometric representation of circumference arcs.
  struct ArcData{
    int m_matid{-10};//< material identifier
    REAL m_radius{-1};//< radius
    REAL m_xc{0};//< x-coordinate of circumference center
    REAL m_yc{0};//< y-coordinate of circumference center
    REAL m_zc{0};//< z-coordinate of circumference center
  };

  /**
     @brief Converts linear line elements to TPZArc3D elements,
     allowing for an exact geometric representation of curved geometries.
     The neighbouring elements will be converted to TPZGeoBlend<T> so that
     they take the curved sides into account.
     @param [in] gmesh geometric mesh to be transformed.
     @param [in] circles data of all circumference arcs.
     @note TPZGeoMesh::BuildConnectivity should have been called beforehand.
   */
  void SetExactArcRepresentation(TPZAutoPointer<TPZGeoMesh> gmesh,
                                 const TPZVec<ArcData> &circles);

  /**
     @brief Refines elements whose neighbours have materials in matids. 
     This function will call
     `TPZRefPatternTools::RefineDirectional(gel, matids);`
     in all elements of the mesh.
     @param gmesh geometrical mesh to be refined
     @param matids elements will be refined towards materials with these identifiers
     @param nrefs number of refinement stages
  **/
void DirectionalRefinement(TPZAutoPointer<TPZGeoMesh> gmesh,
                           std::set<int> matids,const int nrefs);
};

#endif /* _GMESHTOOLS_HPP_ */
