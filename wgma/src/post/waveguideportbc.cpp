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
      //we expect a TPZPlanarWgma material
      const TPZMaterialDataT<CSTATE> &data = eldata;
      auto mat = dynamic_cast<const TPZScalarField*>(eldata.GetMaterial());
      TPZFNMatrix<9,CSTATE> coeff_mat;
      const bool is_te = m_te;
      if(is_te){
        mat->GetPermeability(data.x, coeff_mat);
      }else{
        mat->GetPermittivity(data.x, coeff_mat);
      }
      coeff_mat.Decompose(ELU);

      const auto nsol = data.sol.size();
      const auto detjac = data.detjac;
      const CSTATE cte = weight*fabs(detjac);

      TPZFNMatrix<3000,CSTATE> rot_et(3,nsol,0.);

      const auto &sol = data.sol;
      for(auto isol = 0; isol < nsol; isol++){
        rot_et.Put(1,isol,sol[isol][0]);
      }

      TPZFNMatrix<3000,CSTATE> tmp;

      TPZFNMatrix<3000,CSTATE> test_func(3,nsol,0.);
      
      /*
        the test functions correspond to the tangential field component
        in the bc, however, we always have 1i*beta*et
        so we perform this multiplication now, just to avoid
        yet another loop
      */
      {
        CSTATE *ptr_rot_et = rot_et.Elem();
        CSTATE *ptr_test_func = test_func.Elem();
        const auto nelem = 3*nsol;
        for(int isol = 0; isol < nsol; isol++){
          const auto beta = m_beta[isol];
          //first component
          ptr_test_func++;
          ptr_rot_et++;
          //second component
          *ptr_test_func++ = std::conj(*ptr_rot_et);
          *ptr_rot_et++ *= 1i*beta;
          ///third component
          ptr_test_func++;
          ptr_rot_et++;
        }
      }
      tmp = rot_et;
      coeff_mat.Substitution(&tmp);
      this->m_k_scratch[index].AddContribution(0, 0, test_func, true, tmp, false,cte);
      //src term
      if(is_src){
        //compute solution
        TPZFNMatrix<3,CSTATE> sol_mat(3,1,0.);
        for(auto is = 0; is < nsol; is++){
          const auto coeff = m_coeff[is];
          const auto z = data.x[0];
          //only ZERO
          if(coeff == 0.){continue;}
          const auto beta = m_beta[is];
          for(auto x = 0; x < 3; x++){
            const auto val = sol_mat.Get(x,0);
            const auto src = rot_et.Get(x,is);
            sol_mat.Put(x, 0, val + src*coeff*std::exp(sign*1i*beta*z));
          }
        }
        coeff_mat.Substitution(&sol_mat);
        TPZFMatrix<CSTATE> fmat(nsol,1,this->m_f_scratch[index].begin(),nsol);
        fmat.AddContribution(0, 0, test_func, true, sol_mat, false, cte*2.0);
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
      /*
        in the modal analysis we always compute the +z propagating modes
        we need to transform grad_ez at the -z boundary so we have the
        appropriate -z propagating modes
       */
      for(auto isol = 0; isol < nsol; isol++){
        const CSTATE v0 = sign*rot_grad_ez.Get(0,isol);
        const CSTATE v1 = sign*rot_grad_ez.Get(1,isol);
        rot_grad_ez.Put(0,isol,v1);
        rot_grad_ez.Put(1,isol,-v0);
      }

      TPZFNMatrix<3000,CSTATE> tmp;

      TPZFNMatrix<3000,CSTATE> test_func(3,nsol,0.);
      
      /*
        the test functions correspond to the tangential field component
        in the bc, however, we always have 1i*beta*et
        so we perform this multiplication now, just to avoid
        yet another loop
       */
      {
        CSTATE *ptr_rot_et = rot_et.Elem();
        CSTATE *ptr_test_func = test_func.Elem();
        const auto nelem = 3*nsol;
        for(int isol = 0; isol < nsol; isol++){
          const auto beta = m_beta[isol];
          for(int irow = 0; irow < 3; irow++){
            *ptr_test_func++ = std::conj(*ptr_rot_et);
            *ptr_rot_et++ *= 1i*beta;
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
          const auto z = datavec[0].x[2];
          //only ZERO
          if(coeff == 0.){continue;}
          const auto beta = m_beta[is];
          for(auto x = 0; x < 3; x++){
            const auto val = sol_mat.Get(x,0);
            const auto src = rot_et.Get(x,is)+sign*rot_grad_ez.Get(x,is);
            sol_mat.Put(x, 0, val + src*coeff*std::exp(sign*1i*beta*z));
          }
        }
        ur.Substitution(&sol_mat);
        TPZFMatrix<CSTATE> fmat(nsol,1,this->m_f_scratch[index].begin(),nsol);
        fmat.AddContribution(0, 0, test_func, true, sol_mat, false, cte*2.0);
      }
    }
  }

  template
  class WaveguidePortBC<SingleSpaceIntegrator>;
  template
  class WaveguidePortBC<MultiphysicsIntegrator>;
  
};