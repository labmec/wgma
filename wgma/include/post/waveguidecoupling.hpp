#ifndef _WAVEGUIDECOUPLING_HPP_
#define _WAVEGUIDECOUPLING_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  /** @brief Computes the orthogonality between different waveguide solutions.

      This class provides a tool for analysing the coupling between
      eigenmodes of a waveguide.
      The SetAdjoint option is used to distinguish between two use cases
      1. the mesh stores the solutions of a given waveguide problem
      2. the mesh stores the solution of a given problem AND of its adjoint problem
      Assuming that n is the number of computed modes, in the first case
      (SetAdjoint(false)), the mesh solution is expected to have n columns.
      Otherwise, it is expected to have 2n columns.
  */
  template<class TSPACE>
  class WaveguideCoupling: public TSPACE{
  public:
    WaveguideCoupling(TPZAutoPointer<TPZCompMesh> mesh,
                      std::set<int> matids = {},
                      bool conj=false, int nThreads = 4) :
      TSPACE(mesh,matids,nThreads), m_conj(conj) {}
    WaveguideCoupling(TPZAutoPointer<TPZCompMesh> mesh,
                      TPZVec<TPZCompEl*> elvec,
                      bool conj=false, int nThreads = 4) :
      TSPACE(mesh,elvec,nThreads), m_conj(conj) {}
    //! Computes coupling between original and adjoint problem
    void SetAdjoint(bool m){m_adj = m;}
    //! Computes coupling of all modes in mesh
    void ComputeCoupling();
    //! Gets coupling associated with each mode of the waveguide
    void GetCoupling(TPZFMatrix<CSTATE> &kii) const{
      kii = m_kii;
    }
    //! Initialises element data and sets material
    void InitData(TPZCompEl *el, ElData &data) override;
    //! Sets beta values (needed for 2d cross sections)
    void SetBeta(const TPZVec<CSTATE> &beta){m_beta = beta;}
    //! Gets beta values (needed for 2d cross sections)
    void GetBeta(TPZVec<CSTATE> &beta) const {beta = m_beta;}
  protected:
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;
    //! Temp mat results (one position for thread)
    TPZVec<TPZFMatrix<CSTATE>> m_k_scratch;
    //! Coupling results (coupling between i and j modes)
    TPZFMatrix<CSTATE> m_kii;
    //! Propagation constant beta (needed for 2D modal analysis)
    TPZVec<CSTATE> m_beta;
    //! Whether we are comparing solutions of one problem or problem + adjoint
    bool m_adj{false};
    //! Whether to compute the conjugate cross product
    bool m_conj{true};
  };
};

#endif /* _WAVEGUIDECOUPLING_HPP_ */
