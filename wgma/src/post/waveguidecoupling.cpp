#include "post/waveguidecoupling.hpp"
#include <TPZBndCond.h>
#include <TPZMaterialDataT.h>
#include <Electromagnetics/TPZScalarField.h>
#include <Electromagnetics/TPZWgma.h>
#include <TPZElectromagneticConstants.h>
#include <pzvec_extras.h>
#include <pzaxestools.h>
using namespace std::complex_literals;

namespace wgma::post{

  template<class TSPACE>
  void WaveguideCoupling<TSPACE>::ComputeCoupling(){
    const int size_res = std::max(this->NThreads(),1);
    const int ncols = this->Mesh()->Solution().Cols();
    //if we deal with the adjoint problem as well, we have 2n solutions
    const int nsol = m_adj ? ncols/2 : ncols;
    if constexpr (std::is_same_v<TSPACE,MultiphysicsIntegrator>){
      if(m_beta.size() != 0 && m_beta.size() != nsol){
        DebugStop();
      }
    }
    m_kii.Redim(nsol,nsol);

    m_k_scratch.Resize(size_res);
    for(auto &res : m_k_scratch) {res.Redim(nsol,nsol);}
    
    this->Integrate(this->m_elvec);
    
    for (int ires = 0; ires< size_res; ires++){
      for(auto isol = 0; isol < nsol; isol++){
        for(auto jsol = 0; jsol < nsol; jsol++){
          m_kii(isol,jsol) += m_k_scratch[ires](isol,jsol);
        }
      }
    }
    
    m_k_scratch.Resize(0);
  }

  template<class TSPACE>
  void WaveguideCoupling<TSPACE>::InitData(TPZCompEl *el, ElData &data)
  {
    TSPACE::InitData(el,data);
    data.SetMaterial(el->Material());
  }
  
  template<class TSPACE>
  void WaveguideCoupling<TSPACE>::Compute(const ElData &eldata, REAL weight, int index)
  {
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;
      const auto &sol =  data.sol;
      const auto nsol = m_adj ? sol.size()/2 : sol.size();
      const int firstj = m_adj ? nsol : 0;
      auto mat = dynamic_cast<const TPZScalarField*>(eldata.GetMaterial());
      TPZFNMatrix<9,CSTATE> ur;
      mat->GetPermeability(data.x, ur);
      const CSTATE urval = 1./ur.Get(1,1);
      const CSTATE cte = urval * weight * fabs(data.detjac);
      for(auto is = 0; is < nsol; is++){
        const auto isol = sol[is][0];
        for(auto js = 0; js < nsol; js++){
          const auto jsol = sol[firstj+js][0];
          const CSTATE kval = m_conj ? isol * std::conj(jsol) : isol * jsol;
          this->m_k_scratch[index](is,js) += cte * kval;
        }
      }
    }else{
      //NOT YET TESTED
      //we expect a TPZWgma material
      const TPZVec<TPZMaterialDataT<CSTATE>> &datavec = eldata;

      auto mat = dynamic_cast<const TPZWgma*>(eldata.GetMaterial());
      TPZFNMatrix<9,CSTATE> ur;
      mat->GetPermeability(datavec[0].x, ur);
      ur.Decompose(ELU);
      
      TPZFNMatrix<3,CSTATE> et_i(3,1,0.), et_j(3,1,0);
      TPZFNMatrix<1,CSTATE> tmp(3,3,0);
      
      const int solsize = datavec[0].sol.size();
      const auto nsol = m_adj ? solsize/2 : solsize;
      const int firstj = m_adj ? nsol : 0;
      const auto &axes = datavec[0].axes;
      const auto detjac = datavec[0].detjac;
      const CSTATE cte = weight*fabs(detjac);
      for(auto isol = 0; isol < nsol; isol++){
        auto &et_ref = datavec[ TPZWgma::HCurlIndex() ].sol[isol];
        CSTATE val = 0;
        const auto jbeta = 1i*m_beta[isol];
        for(int ix = 0; ix < 3; ix++) {
          et_i(ix,0)= et_ref[ix]/ m_beta[isol];
        }
        //et_i = ur^-1 et_i
        ur.Substitution(&et_i);
        for(auto jsol = 0; jsol < nsol; jsol++){
          auto &et_ref = datavec[ TPZWgma::HCurlIndex() ].sol[firstj+jsol];
          CSTATE val = 0;
          const auto jbeta = 1i*m_beta[jsol];
          for(int ix = 0; ix < 3; ix++) {
            et_j(ix,0)= et_ref[ix]/ m_beta[jsol];
          }
          this->m_k_scratch[index].AddContribution(isol,jsol,et_i,true,et_j,false,cte);
        }
      }
    }
  }

  template
  class WaveguideCoupling<SingleSpaceIntegrator>;
  template
  class WaveguideCoupling<MultiphysicsIntegrator>;
  
};