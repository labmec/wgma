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
template<class T>
class TPZFMatrix;

class TPZCompMesh;
class TPZGeoMesh;

//! Contains a set of routines for creating and manipulating computational meshes.
namespace wgma::cmeshtools{

  /** @brief Data needed for creating the computational mesh.
      @note for gmsh generated meshes, the routine SetupGmshMaterialData
      fills all the data in the correct format.
      
      Materials related to the variational formulation are given in matinfovec.
      PML materials are in pmlvec, and boundary conditions are in bcvec.
      If no specific computation should take place for a given material id,
      but elements should be created anyway (i.e., probes for easier evaluation
      of the solution), insert the material ids on probevec.
   */
  struct PhysicalData{
    //! each position contains the tuple (id,er,ur) for a given region of the mesh
    std::vector<std::tuple<int,CSTATE,CSTATE>> matinfovec;
    //! each position contains data for a given pml region
    std::vector<TPZAutoPointer<wgma::pml::data>> pmlvec;
    //! each position contains data for a given boundary condition
    std::vector<wgma::bc::data> bcvec;
    //! each position contains ids and dimension of dummy materials (such as probes)
    std::vector<std::pair<int,int>> probevec;
  };

  /**
     @brief Prints the computational mesh in .txt format
     @param [in] cmesh Computational mesh
     @param [in] file name (no extension)
   */
  void PrintCompMesh(TPZAutoPointer<TPZCompMesh> cmesh,
                     std::string filename);
  
  /**
   @brief This function associates materials (regions) read from .msh file
   with the materials to be created in the computational mesh.
   The output of this function is in the expected format for 
   wgma::cmeshtools::CMeshWgma2D.
   alphaPMLx/alphaPMLy are ignored if no PML regions are detected. 
   PML regions' name are mandatory 
   to contain pml and the appropriate pml type (both are case insensitive).
   See cartesian_pml.hpp and cylindrical_pml.hpp for naming conventions.
   e.g. 
   pmlxp    -- ok
   pml_ym   -- ok
   PMLXPYP  -- ok
   PLMX     -- not ok
   
   @param [in] gmshmats vector returned from ReadGmshMesh
   @param [in] matmap for a given region name, contains the pairs (er,ur) for non pml regions
   @param [in] bcmap for a given boundary name, contains the bc type
   @param [in] alphaPML PML attenuation constant in the x,y and z direction (or r, r and z if cylindrical)
   @param [out] data organised data for CMeshWgma2D routine
   @param [in] maxdim maximum dimension to be considered. Leave default for automatic choice.
   @note The .msh file can be read with wgma::gmeshtools::ReadGmshMesh.
*/
void SetupGmshMaterialData(const TPZVec<std::map<std::string,int>> &gmshmats,
                           const std::map<std::string,std::pair<CSTATE,CSTATE>> &matmap,
                           const std::map<std::string,wgma::bc::type> &bcmap,
                           TPZVec<CSTATE> alphaPML,
                           PhysicalData &data,
                           int maxdim=-1);

  /**
     @brief Adds a rectangular PML region to a computational mesh.
     This PML region will consist of the same width and attenuation parameter,
     but can be composed of multiple regions (as in periodic meshes).
     If data.neigh has not been set, neighbours will be automatically identified.
     Therefore, all the other mesh regions (excluding BCs) should have been
     previously inserted, so the function can identify to which region 
     the each of the PML regions is associated.

     @tparam MATVOL material from which the PML inherits
     @param[in] data PML data.
     @param[in] volmats identifiers of valid mesh regions for pml neighbours
     @param[in] gmesh the geometric mesh.
     @param[in] cmesh the computational mesh.
     @return identifier of all PML neighbours
     @note This method calls FindPMLNeighbourMaterial internally.
  */
  template<class MATVOL>
  std::map<int,int>
  AddRectangularPMLRegion(const wgma::pml::cart::data data,
                          const std::set<int> &volmats,
                          TPZAutoPointer<TPZGeoMesh> gmesh,
                          TPZAutoPointer<TPZCompMesh> cmesh);

  /**
     @brief Adds a cylindrical PML region to a computational mesh.
     This PML region will consist of the same width and attenuation parameter,
     but can be composed of multiple regions (as in periodic meshes).
     If data.neigh has not been set, neighbours will be automatically identified.
     Therefore, all the other mesh regions (excluding BCs) should have been
     previously inserted, so the function can identify to which region 
     the each of the PML regions is associated.

     @tparam MATVOL material from which the PML inherits
     @param[in] data PML data.
     @param[in] volmats identifiers of valid mesh regions for pml neighbours
     @param[in] gmesh the geometric mesh.
     @param[in] cmesh the computational mesh.
     @return identifier of all PML neighbours
     @note This method calls FindPMLNeighbourMaterial internally.
  */
  template<class MATVOL>
  std::map<int,int>
  AddCylindricalPMLRegion(const wgma::pml::cyl::data data,
                          const std::set<int> &volmats,
                          TPZAutoPointer<TPZGeoMesh> gmesh,
                          TPZAutoPointer<TPZCompMesh> cmesh);

  /** @brief Finds all connects associated with a dirichlet boundary condition
      @param [in] cmesh Computational mesh
      @param [out] boundConnects boundary connects
      @praam [in] matIds if non-null, only gets connects matching these material identifiers
   */
  void FindDirichletConnects(TPZAutoPointer<TPZCompMesh> cmesh,
                             std::set<int64_t> &boundConnects,
                             const std::set<int> & matIds = {});
  
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

  /**
     @brief Extract solution from the intersection between two meshes
     The function will iterate over the elements from mesh_orig, find
     the correspondent elements in mesh_dest and copy the relevant dofs
     @note The common elements between meshes must have same polynomial order!


     @param [in] mesh_dest Mesh for which the solution is being extracted
     @param [in] mesh_orig Mesh from which the solution is being extracted
     @param [out] sol_dest Solution vector
   */
  void ExtractSolFromMesh(TPZAutoPointer<TPZCompMesh> mesh_dest,
                          TPZAutoPointer<TPZCompMesh> mesh_orig,
                          TPZFMatrix<CSTATE> &sol_dest);
};

#endif /* _CMESHTOOLS_HPP_ */
