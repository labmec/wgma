#ifndef _WAVEGUIDECOUPLING_HPP_
#define _WAVEGUIDECOUPLING_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  /** @brief Computes the cross product between different waveguide solutions.

      One can choose between computing the cross product or the conjugate cross product.
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
    //! Computes contribution of all modes in mesh
    void ComputeCoupling();
    //! Gets contributions associated with each mode of the waveguide
    void GetCoupling(TPZFMatrix<CSTATE> &kii) const{
      kii = m_kii;
    }
    //! Initialises element data and sets material
    void InitData(TPZCompEl *el, ElData &data) override;
    //! Sets beta values
    void SetBeta(const TPZVec<CSTATE> &beta){m_beta = beta;}
    //! Gets beta values
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
    //! Coefficients of source term as combination of modes
    TPZVec<CSTATE> m_coeff;
    //! Whether to compute the conjugate cross product
    bool m_conj{true};
  };
};

#endif /* _WAVEGUIDECOUPLING_HPP_ */
