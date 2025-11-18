#ifndef _EXPORTSOLUTION_HPP_
#define _EXPORTSOLUTION_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  /** @brief Stores the FEM solution at each integration point for given elements
      for easy post processing.
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
    const TPZVec<REAL> &GetCoordinatesAtPoints() const {return m_x;}
  protected:
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;

    //!index corresponding of first point for a given element
    TPZVec<int64_t> m_first_pt_el;
    //!weight of a given integration point
    TPZVec<REAL> m_weights;
    //!solution vector, size = soldim*npts
    TPZVec<CSTATE> m_sol;
    //!position vector, size = 3*npts
    TPZVec<REAL> m_x;
    //!dimension of solution
    int m_sol_dim{-1};
  };
};

#endif /* _EXPORTSOLUTION_HPP_ */
