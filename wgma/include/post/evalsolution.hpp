#ifndef _EVALSOLUTION_HPP_
#define _EVALSOLUTION_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  /** @brief Eval the FEM solution at each integration point for given elements
      and computes minimum and maximum value.
  */
  template<class TSPACE>
  class EvalSolution: public TSPACE{
  public:
    
    EvalSolution(TPZAutoPointer<TPZCompMesh> mesh,
                   std::set<int> matids = {},
                   int nThreads = 4) : TSPACE(mesh,matids,nThreads) {}
    EvalSolution(TPZAutoPointer<TPZCompMesh> mesh,
                   TPZVec<TPZCompEl*> elvec,
                   int nThreads = 4) : TSPACE(mesh,elvec,nThreads) {}

    
    void EvalSolutionAtPoints(STATE &min, STATE &max);
  protected:
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;

    //!maximum solution at points
    TPZVec<STATE> m_max_sol;
    //!minimum solution at points
    TPZVec<STATE> m_min_sol;
    //!dimension of solution
    int m_sol_dim{-1};
  };
};

#endif /* _EVALSOLUTION_HPP_ */
