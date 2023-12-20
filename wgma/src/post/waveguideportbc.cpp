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
    const int ncols = this->Mesh()->Solution().Cols();
    //if we deal with the adjoint problem as well, we have 2n solutions
    const int nsol = m_adj ? ncols/2 : ncols;
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
  }


  #include <Electromagnetics/TPZScalarField.h>
  
  template<class TSPACE>
  void WaveguidePortBC<TSPACE>::Compute(const ElData &eldata, REAL weight, int index)
  {
    const STATE sign = m_pos_z ? 1 : -1;
    const bool is_src = m_coeff.size() != 0;
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
      constexpr bool m_conj = true;
      for(auto is = 0; is < nsol; is++){
        const auto isol = m_conj ? std::conj(sol[is][0]) : sol[is][0];
        for(auto js = 0; js < nsol; js++){
          const auto jsol = sol[firstj+js][0];
          const auto beta = m_beta[js];
          const CSTATE kval = 1i*beta*jsol*isol;
          this->m_k_scratch[index](is,js) += cte * kval;
          if(is_src){
            const CSTATE fval = -2.*1i*beta*m_coeff[js]*jsol*isol;
            this->m_f_scratch[index][is] += cte * fval;
          }
        }
      }
    }else{
      //we expect a TPZWgma material
      const TPZVec<TPZMaterialDataT<CSTATE>> &datavec = eldata;

      auto mat = dynamic_cast<const TPZWgma*>(eldata.GetMaterial());
      TPZFNMatrix<9,CSTATE> ur;
      mat->GetPermeability(datavec[0].x, ur);
      ur.Decompose(ELU);

      const auto nsol = datavec[0].sol.size();
      const auto &axes = datavec[0].axes;
      const auto detjac = datavec[0].detjac;
      const CSTATE cte = weight*fabs(detjac);

      TPZFNMatrix<3000,CSTATE> rot_et(3,nsol,0.);

      TPZFNMatrix<3000,CSTATE> rot_grad_ez(3,nsol,0.);
      TPZFNMatrix<2000,CSTATE> grad_ez_axes(2,nsol,0.);

      for(auto isol = 0; isol < nsol; isol++){
        const auto &et_ref = datavec[ TPZWgma::HCurlIndex() ].sol[isol];
        rot_et.Put(0,isol,et_ref[1]);
        rot_et.Put(1,isol,-et_ref[0]);
        auto &gradez_ref = datavec[ TPZWgma::H1Index() ].dsol[isol];
        grad_ez_axes.Put(0,isol,gradez_ref[0]);
        grad_ez_axes.Put(1,isol,gradez_ref[1]);
      }

      TPZAxesTools<CSTATE>::Axes2XYZ(grad_ez_axes, rot_grad_ez, axes);
      //rotate grad ez
      for(auto isol = 0; isol < nsol; isol++){
        const CSTATE v0 = rot_grad_ez.Get(0,isol);
        const CSTATE v1 = rot_grad_ez.Get(1,isol);
        rot_grad_ez.Put(0,isol,v1);
        rot_grad_ez.Put(1,isol,-v0);
      }

      TPZFNMatrix<3000,CSTATE> tmp;

      TPZFNMatrix<3000,CSTATE> test_func(3,nsol,0.);
      
      
      {
        const CSTATE *ptr_rot_et = rot_et.Elem();
        CSTATE *ptr_test_func = test_func.Elem();
        for(int j = 0; j < nsol; j++){
          for(int i = 0; i < 3; i++){
            *ptr_test_func++ = std::conj(*ptr_rot_et++);
          }
        }
      }
      
      tmp = rot_grad_ez;
      tmp *= sign;
      tmp += rot_et;
      ur.Substitution(&tmp);
      this->m_k_scratch[index].AddContribution(0, 0, test_func, true, tmp, false,cte);
      //src term
      if(is_src){
        //compute solution
        TPZFNMatrix<3,CSTATE> sol_mat(3,1,0.);
        for(auto is = 0; is < nsol; is++){
          const auto coeff = m_coeff[is];
          //only ZERO
          if(coeff == 0.){continue;}
          for(auto x = 0; x < 3; x++){
            const auto val = sol_mat.Get(x,0);
            sol_mat.Put(x,0,val + rot_et.Get(x,is)*coeff);
          }
          ur.Substitution(&sol_mat);
          TPZFMatrix<CSTATE> fmat(nsol,1,this->m_f_scratch[index].begin(),nsol);
          fmat.AddContribution(0, 0, test_func, true, sol_mat, false, cte*2.0);
        }
      }
    }
  }

  template
  class WaveguidePortBC<SingleSpaceIntegrator>;
  template
  class WaveguidePortBC<MultiphysicsIntegrator>;
  
};