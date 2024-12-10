#ifndef _SOLPLANEWAVE_HPP_
#define _SOLPLANEWAVE_HPP_

#include <post/integrator.hpp>


namespace wgma::post{
  /** @brief  Computes the decomposition of two given FEM solutions
      into a constant vector field in the x direction and a constant
      vector field in the y direction.

      This class is specially useful if one is looking for a planewave
      basis in 2D, since we have no control on the direction of the
      vector fields returned by the modal solver
   */
  template<class TSPACE>
  class PlanewaveDecomposition: public TSPACE{
  public:
    explicit PlanewaveDecomposition(TPZAutoPointer<TPZCompMesh> mesh,
                                    std::set<int> matids = {},
                                    int nThreads = 4) :
      TSPACE(mesh,matids,nThreads) {}
    explicit PlanewaveDecomposition(TPZAutoPointer<TPZCompMesh> mesh,
                                    TPZVec<TPZCompEl*> elvec,
                                    int nThreads = 4) :
      TSPACE(mesh,elvec,nThreads) {}
    /** @brief Compute coefficients of decomposition into planewaves.
     */
    void ComputeCoefficients(CSTATE &a1x,CSTATE &a2x, CSTATE &a1y, CSTATE &a2y);
  protected:
    //! Computes contribution at an integration point
    void Compute(const ElData &data, REAL weight, int thread) override;
    //! Matrix results (one matrix for thread)
    TPZVec<TPZFMatrix<CSTATE>> m_mat;
    //! Rhs results (one position for thread)
    TPZVec<TPZFMatrix<CSTATE>> m_rhs;
  };
};

#endif /* _SPARAM_HPP_ */
