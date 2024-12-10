#include "post/planewave.hpp"
#include <TPZMatrixWindow.h>
#include <numeric>

namespace wgma::post{

  template<class TSPACE>
  void PlanewaveDecomposition<TSPACE>::ComputeCoefficients(CSTATE &a1x,CSTATE &a2x,
                                                           CSTATE &a1y, CSTATE &a2y)
  {
    auto mesh = this->Mesh();
    const int size_res = std::max(this->NThreads(),1);
    //we need at least two solutions
    const int nsol = mesh->Solution().Cols();
    if(nsol < 2){DebugStop();}
    
    m_mat.Resize(size_res);
    m_rhs.Resize(size_res);
    for(auto it = 0; it < size_res; it++){
      m_mat[it].Redim(2,2);
      m_rhs[it].Redim(2,2);
    }
    this->Integrate(this->m_elvec);
    
    TPZFMatrix<CSTATE> mat(2,2,0), rhs(2,2,0);
    for(auto it = 0; it < size_res; it++){
      mat += m_mat[it];
      rhs += m_rhs[it];
    }

    mat.SolveDirect(rhs, ELU);
    a1x = rhs.GetVal(0, 0);
    a2x = rhs.GetVal(1, 0);
    a1y = rhs.GetVal(0, 1);
    a2y = rhs.GetVal(1, 1);
  }


  template<class TSPACE>
  void PlanewaveDecomposition<TSPACE>::Compute(const ElData &eldata, REAL weight, int index)
  {
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;
      const CSTATE cte = weight*fabs(data.detjac);

      const auto soldim = data.sol[0].size();
      //we dont care for other solutions
      constexpr int nsol{2};
      TPZFNMatrix<6,CSTATE> solmat(soldim, 2, 0);
      for(int isol = 0; isol < 2; isol++){
        const auto &sol = data.sol[isol];
        for(int ix = 0; ix < soldim; ix++){
          solmat.PutVal(ix,isol,sol[ix]);
        }
      }
      
      TPZFNMatrix<4,CSTATE> rhs(soldim,2,0);
      //first column is a vector field in x direction, second column in y direction
      rhs.PutVal(0,0,1);
      rhs.PutVal(1,1,1);


      constexpr int conj_trans{2};
      constexpr int no_trans{0};
      this->m_mat[index].AddContribution(0,0,solmat,conj_trans,solmat,no_trans,cte);
      this->m_rhs[index].AddContribution(0,0,solmat,conj_trans,rhs,no_trans,cte);
    }else{
      DebugStop();
    }
  }

  template
  class PlanewaveDecomposition<SingleSpaceIntegrator>;
  template
  class PlanewaveDecomposition<MultiphysicsIntegrator>;
  
};