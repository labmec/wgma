#ifndef _WAVEGUIDEPORTBC_HPP_
#define _WAVEGUIDEPORTBC_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  /** @brief Computes the boundary terms for outgoing waveguide port boundary condition

      Combined with the restriction of the approximation space, this class provides
      the contribution for all modes to be added to the FEM matrix

      Using the notation of Jian-ming Jin's Finite Element Method for Electromagnetics,
      chapter Finite Element-Eigenfunction Expansion Methods,
      section on Waveguide Port Boundary Condition,
      this class will compute
      \int_\Gamma P(E)\cdot F d\Gamma

      provided that E,F are the basis functions restricted on the computed modes
      of the waveguide.
   */
  template<class TSPACE>
  class WaveguidePortBC: public TSPACE{
  public:
    WaveguidePortBC(TPZAutoPointer<TPZCompMesh> mesh,
                 std::set<int> matids = {},
                 int nThreads = 4) : TSPACE(mesh,matids,nThreads) {}
    WaveguidePortBC(TPZAutoPointer<TPZCompMesh> mesh,
                 TPZVec<TPZCompEl*> elvec,
                 int nThreads = 4) : TSPACE(mesh,elvec,nThreads) {}
    //! Computes contribution of all modes in mesh
    void ComputeContribution();
    //! Gets contributions associated with each mode of the waveguide
    void GetContribution(TPZVec<CSTATE> &km_vec) {km_vec = m_res;}
    //! Initialises element data and sets material
    void InitData(TPZCompEl *el, ElData &data) override;
    //! Sets beta values
    void SetBeta(const TPZVec<CSTATE> &beta){m_beta = beta;}
    //! Gets beta values
    void GetBeta(TPZVec<CSTATE> &beta) const {beta = m_beta;}
    //! Sets direction
    void SetPositiveZ(const bool b){m_pos_z = b;}
    //! Sets direction
    [[nodiscard]] bool GetPositiveZ() const{return m_pos_z;}
  protected:
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;
    //! Temp results (one position for thread)
    TPZVec<TPZVec<CSTATE>> m_scratch;
    //! Results (one position for solution)
    TPZVec<CSTATE> m_res;
    //! Propagation constant beta (needed for 2D modal analysis)
    TPZVec<CSTATE> m_beta;
    //! Whether boundary is +z or -z
    bool m_pos_z{true};
  };
};

#endif /* _WAVEGUIDEPORTBC_HPP_ */
