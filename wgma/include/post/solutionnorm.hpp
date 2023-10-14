#ifndef _SOLNORM_HPP_
#define _SOLNORM_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  //! Computes L2 norm of a given FEM solution (already loaded into the mesh)
  template<class TSPACE>
  class SolutionNorm: public TSPACE{
  public:
    explicit SolutionNorm(TPZAutoPointer<TPZCompMesh> mesh,
                          std::set<int> matids = {},
                          bool conj = true,
                          int nThreads = 4) :
      TSPACE(mesh,matids,nThreads), m_conj(conj) {}
    explicit SolutionNorm(TPZAutoPointer<TPZCompMesh> mesh,
                          TPZVec<TPZCompEl*> elvec,
                          bool conj = true,
                          int nThreads = 4) :
      TSPACE(mesh,elvec,nThreads), m_conj(conj) {}
    //! Compute norm of a given solution
    CSTATE ComputeNorm(int s);
    //! Compute norm of all solutions
    TPZVec<CSTATE> ComputeNorm();
    //! Normalise all solutions and return their norm before normalising
    TPZVec<CSTATE> Normalise();
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;
    //! Results (one position for thread)
    TPZVec<TPZVec<CSTATE>> m_res;
    //! Whether to compute actual norm (phi*conj(phi))  or not
    bool m_conj{true};
  };
};

#endif /* _SOLNORM_HPP_ */
