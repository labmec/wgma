#ifndef _GMESHTOOLS_HPP_
#define _GMESHTOOLS_HPP_

#include <pmltypes.hpp>

#include <pzvec.h>
#include <pzfmatrix.h>

#include <functional>
#include <optional>

struct SPZPeriodicData;

namespace wgma::gmeshtools{

  enum class ElType {Tri, Quad};
  //! Prints a geometric mesh in .txt and .vtk format
  void PrintGeoMesh(TPZAutoPointer<TPZGeoMesh>& gmesh, std::string filename, bool txt_too=true);

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
     @param[out] periodic_data (see SPZPeriodicData)
     @param[in] verbose whether to output information about read materials
     @note position [x] of matids will contain the materials of dimension x.
   */
  TPZAutoPointer<TPZGeoMesh>
  ReadPeriodicGmshMesh(const std::string filename,
                       const REAL scale,
                       TPZVec<std::map<std::string,int>> & matids,
                       TPZAutoPointer<SPZPeriodicData> &periodic_data,
                       bool verbose = true);

  /**
     @brief Finds the neighbouring material of a given cartesian pml region.
     @param[in] gmesh geometrical mesh.
     @param[in] pmlDim dimension of pml region
     @param[in] pmlId pml material identifier.
     @param[in] boundPosX x-coordinate of the pml with the inner domain
     @param[in] boundPosY y-coordinate of the pml with the inner domain
     @param[in] boundPosZ z-coordinate of the pml with the inner domain
  */
  [[nodiscard]] std::optional<int>
  FindCartPMLNeighbourMaterial(TPZAutoPointer<TPZGeoMesh>& gmesh,
                               const int pmlDim,
                               const int pmlId,
                               const std::set<int> &volmats,
                               const REAL boundPosX,
                               const REAL boundPosY,
                               const REAL boundPosZ);
  /**
     @brief Finds the neighbouring material of a given cylindrical pml region.
     @param[in] gmesh geometrical mesh.
     @param[in] pmlDim dimension of pml region
     @param[in] pmlId pml material identifier.
     @param[in] boundPosR radial coordinate of the pml with the inner domain
     @param[in] boundPosZ z-coordinate of the pml with the inner domain
  */
  [[nodiscard]] std::optional<int>
  FindCylPMLNeighbourMaterial(TPZAutoPointer<TPZGeoMesh>& gmesh,
                              const int pmlDim,
                              const int pmlId,
                              const std::set<int> &volmats,
                              const REAL boundPosR,
                              const REAL boundPosZ);
  
  /**
     @brief Finds the width of a given pml region
     @param[in] gmesh geometrical mesh
     @param[in] pmlId identifiers of pml region
     @param[in] type attenuation direction of the PML
     @param[out] boundPosX X coordinate of the PML interface with domain (for y-attenuating PMLS, X coordinate of the center of the PML)
     @param[out] pmlWidthX PML width in the X direction (only for relevant types)
     @param[out] boundPosY Y coordinate of the PML interface with domain (for x-attenuating PMLS, X coordinate of the center of the PML)
     @param[out] pmlWidthY PML width in the Y direction (only for relevant types)
     @param[out] boundPosZ Z coordinate of the PML interface with domain (for x-attenuating PMLS, X coordinate of the center of the PML)
     @param[out] pmlWidthZ PML width in the Z direction (only for relevant types)
  */
  void
  FindPMLWidth(TPZAutoPointer<TPZGeoMesh>& gmesh,
               const std::set<int> pmlId, const wgma::pml::cart::type type,
               REAL &boundPosX, REAL &pmlWidthX,
               REAL &boundPosY, REAL &pmlWidthY,
               REAL &boundPosZ, REAL &pmlWidthZ);

  /**
     @brief Finds the dimensions of a given cylindrical pml region
     @param[in] gmesh geometrical mesh
     @param[in] pmlId identifiers of pml region
     @param[in] type attenuation direction of the PML
     @param[out] rmin minimum radius
     @param[out] rmax maximum radius
     @param[out] boundPosZ Z coordinate of the PML interface with domain (for x-attenuating PMLS, X coordinate of the center of the PML)
     @param[out] pmlWidthZ PML width in the Z direction (only for relevant types)
  */
  void
  FindPMLWidth(TPZAutoPointer<TPZGeoMesh>& gmesh,
               const std::set<int> pmlId, const wgma::pml::cyl::type type,
               REAL &rmin, REAL &rmax,
               REAL &boundPosZ, REAL &pmlWidthZ);
  
  /**
     @brief Finds a valid neighbouring material for a given boundary material.
     Only neighbours with identifiers in volmats will be considered.
     @param[in] gmesh geometric mesh
     @param[in] mat material whose neighbour we are looking for
     @param[in] volmats possible ids for neighbour candidates
  */
  [[nodiscard]] std::optional<int> FindBCNeighbourMat(TPZAutoPointer<TPZGeoMesh>& gmesh,
                                                      const int mat,
                                                      const std::set<int> &volmats);

  //! Stores data for allowing exact geometric representation of circumference arcs.
  struct ArcData{
    int m_matid{-10};//< material identifier
    REAL m_radius{-1};//< radius
    REAL m_xc{0};//< x-coordinate of circumference center
    REAL m_yc{0};//< y-coordinate of circumference center
    REAL m_zc{0};//< z-coordinate of circumference center
  };
  //! Stores data for allowing exact geometric representation of cylinder walls
  struct CylinderData{
    int m_matid{-10};//< material identifier
    REAL m_radius{-1};//< radius
    REAL m_xc{0};//< x-coordinate of cylinder center
    REAL m_yc{0};//< y-coordinate of cylinder center
    REAL m_zc{0};//< z-coordinate of cylinder center
    REAL m_xaxis{0};//<x-coordinate of cylinder axis
    REAL m_yaxis{0};//<y-coordinate of cylinder axis
    REAL m_zaxis{0};//<z-coordinate of cylinder axis
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
  void SetExactArcRepresentation(TPZAutoPointer<TPZGeoMesh>& gmesh,
                                 const TPZVec<ArcData> &circles);

  /**
     @brief Converts linear 2D elements to TPZCylinderMap<T> elements,
     allowing for an exact geometric representation of curved walls of a cylinder.
     The neighbouring elements will be converted to TPZGeoBlend<T> so that
     they take the curved sides into account.
     @param [in] gmesh geometric mesh to be transformed.
     @param [in] cylinders data of all cylinders
     @note TPZGeoMesh::BuildConnectivity should have been called beforehand.
   */
  void SetExactCylinderRepresentation(TPZAutoPointer<TPZGeoMesh>& gmesh,
                                      const TPZVec<CylinderData> &cylinders);
  
  /**
     @brief Refines elements whose neighbours have materials in matids. 
     This function will call
     `TPZRefPatternTools::RefineDirectional(gel, matids);`
     in all elements of the mesh.
     @param gmesh geometrical mesh to be refined
     @param matids elements will be refined towards materials with these identifiers
     @param nrefs number of refinement stages
  **/
void DirectionalRefinement(TPZAutoPointer<TPZGeoMesh>& gmesh,
                           std::set<int> matids,const int nrefs);
  /**
     @brief Rotates all node-coordinates in the mesh around a given axis
     @param gmesh geometrical mesh
     @param axis rotation axis
     @param theta rotation angle
   */
  void RotateMesh(TPZAutoPointer<TPZGeoMesh>& gmesh,
                  const TPZVec<REAL> &axis, const REAL theta);

  /**
     @brief Compute map of periodic elements for materials given in desired_mats
     @param [in] gmesh geometrical mesh
     @param [in] desired_mats periodic materials of interest
     @param [in] data periodic data obtained from gmsh
     @param [out] periodic_els mapping of periodic elements
   */
  void GetPeriodicElements(TPZGeoMesh *gmesh,
                           const TPZVec<std::pair<int,int>> &desired_mats,
                           const TPZAutoPointer<SPZPeriodicData> &data,
                           TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &periodic_els);
};

#endif /* _GMESHTOOLS_HPP_ */
