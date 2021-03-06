#ifndef _CMESHTOOLS_HPP_
#define _CMESHTOOLS_HPP_

#include "bctypes.hpp"
#include "pmltypes.hpp"
#include "modetypes.hpp"
#include <pzreal.h>
#include <set>
#include <vector>
#include <map>

template<class T>
class TPZVec;
template<class T>
class TPZAutoPointer;

class TPZCompMesh;
class TPZGeoMesh;

//! Contains a set of routines for creating and manipulating computational meshes.
namespace wgma::cmeshtools{

  /** @brief Data needed for creating the computational mesh.
      @note for gmsh generated meshes, the routine SetupGmshMaterialData
      fills all the data in the correct format.
   */
  struct PhysicalData{
    //! each position contains the tuple (id,er,ur) for a given region of the mesh
    std::vector<std::tuple<int,CSTATE,CSTATE>> matinfovec;
    //! each position contains data for a given pml region
    std::vector<wgma::pml::data> pmlvec;
    //! each position contains data for a given boundary condition
    std::vector<wgma::bc::data> bcvec;
  };
  /**
   @brief This function associates materials (regions) read from .msh file
   with the materials to be created in the computational mesh.
   The output of this function is in the expected format for 
   wgma::cmeshtools::CMeshWgma2D.
   alphaPMLx/alphaPMLy are ignored if no PML regions are detected. 
   PML regions' name are mandatory 
   to contain pml and the appropriate pml type (both are case insensitive).
   See pmltypes.hpp for naming conventions.
   e.g. 
   pmlxp    -- ok
   pml_ym   -- ok
   PMLXPYP  -- ok
   PLMX     -- not ok
   
   @param [in] gmshmats vector returned from ReadGmshMesh
   @param [in] matmap for a given region name, contains the pairs (er,ur) for non pml regions
   @param [in] bcmap for a given boundary name, contains the bc type
   @param [in] alphaPMLx PML attenuation constant in the x direction
   @param [in] alphaPMLy PML attenuation constant in the y direction
   @param [out] data organised data for CMeshWgma2D routine
   @note The .msh file can be read with wgma::gmeshtools::ReadGmshMesh.
*/
void SetupGmshMaterialData(const TPZVec<std::map<std::string,int>> &gmshmats,
                           const std::map<std::string,std::pair<CSTATE,CSTATE>> &matmap,
                           const std::map<std::string,wgma::bc::type> &bcmap,
                           const STATE alphaPMLx,
                           const STATE alphaPMLy,
                           PhysicalData &data);

  /**
     @brief Adds a rectangular PML region to a computational mesh.
     All the other mesh regions (excluding BCs) should have been previously inserted,
     so the function can identify to which region the PML is associated.
     @tparam MATVOL material from which the PML inherits
     @param[in] matId pml material identifier.
     @param[in] alphax attenuation constant in the x-direction.
     @param[in] alphay attenuation constant in the y-direction.
     @param[in] type pml type.
     @param[in] volmats identifiers of valid mesh regions for pml neighbours
     @param[in] gmesh the geometric mesh.
     @param[in] cmesh the computational mesh.
     @return identifier of PML neighbour
     @note This method calls FindPMLNeighbourMaterial internally.
  */
  template<class MATVOL>
  int
  AddRectangularPMLRegion(const int matId,
                          const STATE alphax,const STATE alphay,
                          const wgma::pml::type type,
                          const std::set<int> &volmats,
                          TPZAutoPointer<TPZGeoMesh> gmesh,
                          TPZAutoPointer<TPZCompMesh> cmesh);
  
  /** @brief Gets the indices of the equations associated with dirichlet homogeneous BCs
      so they can be filtered out of the global system*/
  void FilterBoundaryEquations(TPZAutoPointer<TPZCompMesh> cmesh,
                               TPZVec<int64_t> &activeEquations,
                               std::set<int64_t>  &boundConnects);

  /**
     @brief Ensures correct destruction of a periodic computational mesh.
     This routine should be called before destroying a periodic mesh in
     order to avoid run-time errors.
   */
  void RemovePeriodicity(TPZAutoPointer<TPZCompMesh> cmesh);
};

#endif /* _CMESHTOOLS_HPP_ */
