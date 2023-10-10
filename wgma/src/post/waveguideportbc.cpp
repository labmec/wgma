#include "post/waveguideportbc.hpp"
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
  void WaveguidePortBC<TSPACE>::ComputeContribution(){
    const int size_res = std::max(this->NThreads(),1);
    
    const int nsol = this->Mesh()->Solution().Cols();
    if(m_coeff.size() != 0 && m_coeff.size() != nsol){
      DebugStop();
    }
    m_kii.Redim(nsol,nsol);

    m_k_scratch.Resize(size_res);
    for(auto &res : m_k_scratch) {res.Redim(nsol,nsol);}

    m_fi.Resize(nsol,0.);
    m_f_scratch.Resize(size_res);
    for(auto &res : m_f_scratch) {res.Resize(nsol,0);}
    
    this->Integrate(this->m_elvec);
    
    for (int ires = 0; ires< size_res; ires++){
      for(auto isol = 0; isol < nsol; isol++){
        m_fi[isol] += m_f_scratch[ires][isol];
        for(auto jsol = 0; jsol < nsol; jsol++){
          m_kii(isol,jsol) += m_k_scratch[ires](isol,jsol);
        }
      }
    }
    
    m_k_scratch.Resize(0);
    m_f_scratch.Resize(0);
  }

  template<class TSPACE>
  void WaveguidePortBC<TSPACE>::InitData(TPZCompEl *el, ElData &data)
  {
    TSPACE::InitData(el,data);
    data.SetMaterial(el->Material());
  }
  
  template<class TSPACE>
  void WaveguidePortBC<TSPACE>::Compute(const ElData &eldata, REAL weight, int index)
  {
    const STATE sign = m_pos_z ? 1 : -1;
    const bool is_src = m_coeff.size() != 0;
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;
      const auto &sol =  data.sol;
      const auto solsize = sol.size();
      if(solsize != m_beta.size()){
        DebugStop();
      }
      STATE val = 0;
      for(auto isol = 0; isol < solsize; isol++){
        for(auto jsol = 0; jsol < solsize; jsol++){
          //matrix term
          const auto beta = m_beta[jsol]*sign;
          const CSTATE kval =  1i*beta* sol[jsol][0] * std::conj(sol[isol][0]);
          this->m_k_scratch[index](isol,jsol) += weight * fabs(data.detjac) * kval;
        }
        //load vector term
        if(is_src){
          const auto beta = m_beta[isol]*sign;
          const CSTATE kval =  1i*beta* sol[isol][0] * std::conj(sol[isol][0]);
          const CSTATE fval = 2.*m_coeff[isol]*kval;
          this->m_f_scratch[index][isol] += weight * fabs(data.detjac) * fval;
        }
      }
    }else{
      //we expect a TPZWgma material
      const TPZVec<TPZMaterialDataT<CSTATE>> &datavec = eldata;

      TPZManVector<CSTATE,3> et(3,0.);
      TPZFNMatrix<3,CSTATE> grad_ez(3,1,0.);

      const int solsize = datavec[0].sol.size();
      const auto &axes = datavec[0].axes;
      const auto detjac = datavec[0].detjac;
      STATE val = 0;
      for(auto isol = 0; isol < solsize; isol++){
        et = datavec[ TPZWgma::HCurlIndex() ].sol[isol];
        const TPZFMatrix<CSTATE> &dsoldaxes = datavec[ TPZWgma::H1Index()].dsol[isol];
        TPZAxesTools<CSTATE>::Axes2XYZ(dsoldaxes, grad_ez, axes);
        CSTATE val = 0;
        const auto beta = m_beta[isol];
        for(int ix = 0; ix < 3; ix++) {
          grad_ez(ix,0) *= 1i;
          et[ix] /= beta;
          val += (1i*beta*et[ix]+sign*grad_ez.Get(ix,0))*et[ix];
        }
        this->m_k_scratch[index](isol,isol) += weight * fabs(detjac) * val;
      }
    }
  }

  template
  class WaveguidePortBC<SingleSpaceIntegrator>;
  template
  class WaveguidePortBC<MultiphysicsIntegrator>;
  
};