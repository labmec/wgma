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
    m_k_scratch.Resize(size_res);
    
    const int nsol = this->Mesh()->Solution().Cols();
    m_kii.Resize(nsol);
    m_kii.Fill(0);
    for(auto &res : m_k_scratch) {res.Resize(nsol);}
    this->Integrate(this->m_elvec);
    
    for (const auto &res : m_k_scratch){
      for(auto isol = 0; isol < nsol; isol++){
        m_kii[isol] += res[isol];
      }
    }
    m_k_scratch.Resize(0);
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
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;
      const auto &sol =  data.sol;
      const auto solsize = sol.size();
      if(solsize != m_beta.size()){
        DebugStop();
      }
      STATE val = 0;
      for(auto isol = 0; isol < solsize; isol++){
        const auto beta = m_beta[isol];
        const CSTATE val =  1i*beta* sol[isol][0] * sol[isol][0];
        this->m_k_scratch[index][isol] += weight * fabs(data.detjac) * val;
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
        this->m_k_scratch[index][isol] += weight * fabs(detjac) * val;
      }
    }
  }

  template
  class WaveguidePortBC<SingleSpaceIntegrator>;
  template
  class WaveguidePortBC<MultiphysicsIntegrator>;
  
};