#ifndef _SOLCOUPLING_HPP_
#define _SOLCOUPLING_HPP_

#include <post/integrator.hpp>


namespace wgma::post{

  /** @brief Holds additional data on element level
      Used for debugging
   */
  class SolCouplData : public ElData{
  public:
    TPZFNMatrix<1000,CSTATE> m_elmat;
    int64_t m_elindex{-1};
  };
  
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
  class SolCoupling: public TSPACE{
  public:
    SolCoupling(TPZAutoPointer<TPZCompMesh> mesh,
                      std::set<int> matids = {},
                      int nThreads = 4) :
      TSPACE(mesh,matids,nThreads) {}
    SolCoupling(TPZAutoPointer<TPZCompMesh> mesh,
                      TPZVec<TPZCompEl*> elvec,
                      int nThreads = 4) :
      TSPACE(mesh,elvec,nThreads) {}
    //! Computes coupling of all modes in mesh
    void ComputeCoupling();
    //! Gets coupling associated with each mode of the waveguide
    void GetCoupling(TPZFMatrix<CSTATE> &kii) const{
      kii = m_kii;
    }
    ElData* CreateElData() override {return new SolCouplData;}
    //! Initialises element data and sets material
    void InitData(TPZCompEl *el, ElData &data) override;
    //! Exports element matrix
    void PostProcessData(ElData& data) override;
    //! Use this for printing element matrices
    void SetPrintMats(bool print){m_print_mats=print;}
    //! Prefix for file name for element matrices(normally output directory)
    void SetFilePrefix(const std::string &name){m_prefix=name;}
  protected:
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;
    //! Temp mat results (one position for thread)
    TPZVec<TPZFMatrix<CSTATE>> m_k_scratch;
    //! Coupling results (coupling between i and j modes)
    TPZFMatrix<CSTATE> m_kii;
    //! Whether to print element matrices (for debugging)
    bool m_print_mats{false};
    //! Prefix for element matrices files
    std::string m_prefix{""};
  };
};

#endif /* _SOLCOUPLING_HPP_ */
