#ifndef _ORTHO_HPP_
#define _ORTHO_HPP_

#include <post/integrator.hpp>

namespace wgma::post{
  //! Orthogonalises a given set of solutions (already loaded into the mesh)
  class OrthoSol: public SingleSpaceIntegrator{
  public:
    OrthoSol(TPZAutoPointer<TPZCompMesh> mesh,
             std::set<int> matids = {},
             int nThreads = 4) :
      SingleSpaceIntegrator(mesh,matids,nThreads) {}
    OrthoSol(TPZAutoPointer<TPZCompMesh> mesh,
             TPZVec<TPZCompEl*> elvec,
             int nThreads = 4) :
      SingleSpaceIntegrator(mesh,elvec,nThreads) {}
    //! Orthogonalise all solutions and return them
    TPZFMatrix<CSTATE> Orthogonalise();
  protected:
    //! Index of the solution being orthogonalised
    void SetSol(int s){m_which = s;}
    //! Gets index of the solution being orthogonalised
    int WhichSol() const {return m_which;}
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;
    //! Which solution is being calculated
    int m_which{-1};
    //! Results (one row per thread, one column per solution)
    TPZFMatrix<CSTATE> m_res;
  };
};

#endif /* _ORTHO_HPP_ */
