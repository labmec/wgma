#ifndef _SOLNORM_HPP_
#define _SOLNORM_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  //! Computes L2 norm of a given FEM solution (already loaded into the mesh)
  template<class TSPACE>
  class SolutionNorm: public TSPACE{
  public:
    SolutionNorm(TPZAutoPointer<TPZCompMesh> mesh,
                 std::set<int> matids = {},
                 int nThreads = 4) : TSPACE(mesh,matids,nThreads) {}
    SolutionNorm(TPZAutoPointer<TPZCompMesh> mesh,
                 TPZVec<TPZCompEl*> elvec,
                 int nThreads = 4) : TSPACE(mesh,elvec,nThreads) {}
    //! Compute norm of a given solution
    STATE ComputeNorm(int s);
    //! Compute norm of all solutions
    TPZVec<STATE> ComputeNorm();
    //! Normalise all solutions and return their norm before normalising
    TPZVec<STATE> Normalise();
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;
    //! Results (one position for thread)
    TPZVec<TPZVec<STATE>> m_res;
  };
};

#endif /* _SOLNORM_HPP_ */
