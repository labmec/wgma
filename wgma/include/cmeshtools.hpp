#ifndef _CMESHTOOLS_HPP_
#define _CMESHTOOLS_HPP_

#include "bctypes.hpp"

#include <pzreal.h>

#include <set>

template<class T>
class TPZVec;
template<class T>
class TPZAutoPointer;

class TPZCompMesh;
class TPZGeoMesh;

//! Contains a set of routines for creating and manipulating computational meshes.
namespace wgma::cmeshtools{

  /**
     @brief Creates the computational meshes used for approximating the waveguide EVP.
     Three meshes will be created: one for the H1 approximation space, one for the
     HCurl approximation space and one multiphysics mesh combining both spaces.
     @note The material ids should be the same given to CreateGMeshRectangularWaveguide.
  */
  TPZVec<TPZAutoPointer<TPZCompMesh>>
  CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
              const TPZVec<int> &matIdVec, const CSTATE ur, const CSTATE er,
              const STATE lambda, const REAL &scale, bool usingSymmetry, bctype sym);
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
