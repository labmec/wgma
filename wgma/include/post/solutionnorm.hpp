#ifndef _SOLNORM_HPP_
#define _SOLNORM_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  //! Computes L2 norm of a given FEM solution (already loaded into the mesh)
  class SolutionNorm: public SingleSpaceIntegrator{
  public:
    SolutionNorm(TPZAutoPointer<TPZCompMesh> mesh,
                 std::set<int> matids = {},
                 int nThreads = 4) : SingleSpaceIntegrator(mesh,matids,nThreads) {}
    SolutionNorm(TPZAutoPointer<TPZCompMesh> mesh,
                 TPZVec<TPZCompEl*> elvec,
                 int nThreads = 4) : SingleSpaceIntegrator(mesh,elvec,nThreads) {}
    //! Compute norm of a given solution
    STATE ComputeNorm(int s = 0);
  protected:
    //! Index of the solution to be used
    void SetSol(int s){m_which = s;}
    //! Gets index of the solution to be used
    int WhichSol() const {return m_which;}
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;
    //! Which solution is being calculated
    int m_which{-1};
    //! Results (one position for thread)
    TPZVec<REAL> m_res;
  };
};

#endif /* _SOLNORM_HPP_ */
