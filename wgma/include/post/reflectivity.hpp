#ifndef _SOLREFLECT_HPP_
#define _SOLREFLECT_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  /** @brief  Computes the reflectivity  of a given FEM solution w.r.t
      a given source.
      The FEM solution is assumed to be the first column of the mesh solution,
      while the next column is the source
      This class will compute the parameter

      \int_k (E_c-E_s)\cdot E_s^* dA
      \int_k E_s\cdot E_s^* dA
      where E_c is the computed field and E_s the source field

   */
  template<class TSPACE>
  class SolutionReflectivity: public TSPACE{
  public:
    explicit SolutionReflectivity(TPZAutoPointer<TPZCompMesh> mesh,
                          std::set<int> matids = {},
                          int nThreads = 4) :
      TSPACE(mesh,matids,nThreads) {}
    explicit SolutionReflectivity(TPZAutoPointer<TPZCompMesh> mesh,
                          TPZVec<TPZCompEl*> elvec,
                          int nThreads = 4) :
      TSPACE(mesh,elvec,nThreads) {}
    //! Compute reflectivity of a solution w.r.t. a source
    CSTATE ComputeReflectivity();
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;
    //! Denominator results (one position for thread)
    TPZVec<CSTATE> m_denominator;
    //! Numerator results (one position for thread)
    TPZVec<CSTATE> m_numerator;
    //! For which modes the parameters should be computed
    TPZVec<int> m_modes;
  };
};

#endif /* _SPARAM_HPP_ */
