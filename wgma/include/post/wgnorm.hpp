#ifndef _WGNORM_HPP_
#define _WGNORM_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  //! Computes L2 norm of a given FEM solution (already loaded into the mesh)
  template<class TSPACE>
  class WgNorm: public TSPACE{
  public:
    explicit WgNorm(TPZAutoPointer<TPZCompMesh> mesh,
                          std::set<int> matids = {},
                          bool conj = true,
                          int nThreads = 4) :
      TSPACE(mesh,matids,nThreads), m_conj(conj) {}
    explicit WgNorm(TPZAutoPointer<TPZCompMesh> mesh,
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
    //! Sets TE/TM
    void SetTE(bool te){m_te = te;}
    //! Gets TE/TM
    [[nodiscard]] bool GetTE() const{return m_te;}
    //! Sets beta values
    void SetBeta(const TPZVec<CSTATE> &beta){m_beta = beta;}
    //! Gets beta values
    void GetBeta(TPZVec<CSTATE> &beta) const {beta = m_beta;}
    //! Sets wl values
    void SetWavelength(const STATE wl){m_wl = wl;}
    //! Gets wl values
    void GetWavelength(STATE wl) const {wl = m_wl;}
  protected:
    //! Results (one position for thread)
    TPZVec<TPZVec<CSTATE>> m_res;
    //! Whether to compute actual norm (phi*conj(phi))  or not
    bool m_conj{true};
    //! Whether solving TE/TM case in 1D
    bool m_te{true};
    //! Propagation constant beta (needed for 2D modal analysis)
    TPZVec<CSTATE> m_beta;
    //! Wavelength
    STATE m_wl{1};
  };
};

#endif /* _SOLNORM_HPP_ */
