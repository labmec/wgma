#ifndef _ORTHO_HPP_
#define _ORTHO_HPP_

#include <post/integrator.hpp>

namespace wgma::post{
  //! Orthogonalises a given set of solutions (already loaded into the mesh)
  template<class TSPACE>
  class OrthoSol: public TSPACE{
  public:
    explicit OrthoSol(TPZAutoPointer<TPZCompMesh> mesh,
             std::set<int> matids = {},
             bool conj=true,
             int nThreads = 4) :
      TSPACE(mesh,matids,nThreads),m_conj(conj) {}
    explicit OrthoSol(TPZAutoPointer<TPZCompMesh> mesh,
             TPZVec<TPZCompEl*> elvec,
             bool conj=true,
             int nThreads = 4) :
      TSPACE(mesh,elvec,nThreads), m_conj(conj) {}
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
    //! Whether to orthogonalise relative to complex conjugate
    bool m_conj{true};
  };
};

#endif /* _ORTHO_HPP_ */
