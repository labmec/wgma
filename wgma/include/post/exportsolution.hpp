#ifndef _EXPORTSOLUTION_HPP_
#define _EXPORTSOLUTION_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  /** @brief Stores the FEM solution at each integration point for given elements
      for easy post processing.
      The solution format is
      geo_el_idx, glob_pt_idx, weight, sol1, sol2, ...

      where glob_pt_idx is the accumulated loc_pt_idx of each element
  */
  template<class TSPACE>
  class ExportSolution: public TSPACE{
  public:
    
    ExportSolution(TPZAutoPointer<TPZCompMesh> mesh,
                   std::set<int> matids = {},
                   int nThreads = 4) : TSPACE(mesh,matids,nThreads) {}
    ExportSolution(TPZAutoPointer<TPZCompMesh> mesh,
                   TPZVec<TPZCompEl*> elvec,
                   int nThreads = 4) : TSPACE(mesh,elvec,nThreads) {}
    
    void StoreSolutionAtPoints(const int sol_dim);
    const TPZVec<CSTATE> &GetSolutionAtPoints() const {return m_sol;}
    const TPZVec<REAL> &GetIntWeightAtPoints() const {return m_weights;}
  protected:
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;

    //!index corresponding of first point for a given element
    TPZVec<int64_t> m_first_pt_el;
    //!weight of a given integration point
    TPZVec<REAL> m_weights;
    //!solution vector, size = soldim*npts
    TPZVec<CSTATE> m_sol;
    //!dimension of solution
    int m_sol_dim{-1};
  };
};

#endif /* _EXPORTSOLUTION_HPP_ */
