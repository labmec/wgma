#ifndef _CMESHTOOLS_HPP_
#define _CMESHTOOLS_HPP_

#include "bctypes.hpp"
#include "pmltypes.hpp"
#include <pzreal.h>

#include <set>
#include <map>

template<class T>
class TPZVec;
template<class T>
class TPZAutoPointer;

class TPZCompMesh;
class TPZGeoMesh;

//! Contains a set of routines for creating and manipulating computational meshes.
namespace wgma::cmeshtools{

  /**
   @brief This function associates materials (regions) read from .msh file
   with the materials to be created in the computational mesh.
   The output of this function is in the expected format for wgma::cmeshtools::CreateCMesh.
   alphaPML is ignored if no PML regions are detected. PML regions' name are mandatory 
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
   @param [in] alphaPML PML attenuation constant
   @param [out] volmatids material identifiers of non-pml regions
   @param [out] ervec electric permittivity of non-pml regions
   @param [out] urvec magnetic permeability of non-pml regions
   @param [out] pmlvec data of pml regions
   @param [out] bcvec data of bc regins

   @note The .msh file can be read with wgma::gmeshtools::ReadGmshMesh.
*/
void SetupGmshMaterialData(const TPZVec<std::map<std::string,int>> &gmshmats,
                       const std::map<std::string,std::pair<CSTATE,CSTATE>> &matmap,
                       const std::map<std::string,wgma::bc::type> &bcmap,
                       const STATE alphaPML,
                       TPZVec<int> &volmatids,
                       TPZVec<CSTATE> &ervec,
                       TPZVec<CSTATE> &urvec,
                       TPZVec<wgma::pml::data> &pmlvec,
                       TPZVec<wgma::bc::data> &bcvec);
  
  /**
     @brief Creates the computational meshes used for approximating the waveguide EVP.
     Three meshes will be created: one for the H1 approximation space, one for the
     HCurl approximation space and one multiphysics mesh combining both spaces.
     If `pmlDataVec.size()>0`, the PML regions will be created and their 
     corresponding domain regions automatically identified.
     @param[in] gmesh geometrical mesh
     @param[in] pOrder polynomial order
     @param[in] volMatIdVec material identifiers of domain regions (no PML or BC)
     @param[in] urVec magnetic permeability of domain regions
     @param[in] erVec electric permittivity of domain regions
     @param[in] pmlDataVec array with PML data
     @param[in] bcDataVec array with BC data
     @param[in] lambda operational wavelength
     @param[in] scale geometric scaling for better floating point precision
  */
  TPZVec<TPZAutoPointer<TPZCompMesh>>
  CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
              const TPZVec<int> &volMatIdVec, const TPZVec<CSTATE> &urVec,
              const TPZVec<CSTATE> &erVec, const TPZVec<pml::data> &pmlDataVec,
              const TPZVec<bc::data> &bcDataVec, const STATE lambda, const REAL &scale);

  /**
     @brief Adds a rectangular PML region to a computational mesh.
     All the other mesh regions (excluding BCs) should have been previously inserted,
     so the function can identify to which region the PML is associated.
     @param[in] matId pml material identifier.
     @param[in] alpha attenuation constant.
     @param[in] type pml type.
     @param[in] cmesh the computational mesh.
     @return This method calls FindPMLNeighbourMaterial internally.
  */
  void
  AddRectangularPMLRegion(const int matId, const int alpha,
                          const wgma::pml::type type,
                          TPZAutoPointer<TPZCompMesh> cmesh);

  /**
     @brief Finds the neighbouring material of a given pml region.
     @param[in] gmesh geometrical mesh.
     @param[in] pmlId pml material identifier.
     @param[in] boundPosX x-coordinate of the pml with the inner domain
     @param[in] boundPosY y-coordinate of the pml with the inner domain
  */
  int
  FindPMLNeighbourMaterial(TPZAutoPointer<TPZGeoMesh> gmesh,const int pmlId,
                           const REAL boundPosX, const REAL boundPosY);
  
  /** @brief Counts active equations per approximation space*/
  void CountActiveEquations(TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec,
                            const std::set<int64_t> &boundConnects,
                            int &neq,
                            int &nH1Equations, int &nHCurlEquations);
  /** @brief Gets the indices of the equations associated with dirichlet homogeneous BCs
      so they can be filtered out of the global system*/
  void FilterBoundaryEquations(TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec,
                               TPZVec<int64_t> &activeEquations, int &neq,
                               int &neqOriginal, int& nh1, int &nhcurl);
};

#endif /* _CMESHTOOLS_HPP_ */
