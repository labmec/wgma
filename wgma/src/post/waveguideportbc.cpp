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

    m_fi.Resize(nsol);
    m_fi.Fill(0);
    m_f_scratch.Resize(size_res);
    for(auto &res : m_f_scratch) {res.Resize(nsol);res.Fill(0);}
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
    constexpr int no_trans{0};
    constexpr int conj_trans{2};
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

      /*
        rot et is the rotated E_t
        et_beta is the rotated E_t times j*beta
       */
      TPZFNMatrix<3000,CSTATE> rot_et(3,nsol,0.), et_beta(3,nsol,0.);

      const auto &sol = data.sol;
      for(auto isol = 0; isol < nsol; isol++){
        const auto &solval = sol[isol][0];
        const auto beta = m_beta[isol];
        rot_et.Put(1,isol,solval);
        et_beta.Put(1,isol,1i*beta*solval);
      }
      
      TPZFNMatrix<3000,CSTATE> tmp(et_beta);
      coeff_mat.Substitution(&tmp);
      TPZFMatrix<CSTATE> &kmat = this->m_k_scratch[index];
      kmat.AddContribution(0, 0, rot_et, conj_trans, tmp, no_trans,cte);
      //src term
      if(is_src){
        //compute solution
        TPZFNMatrix<3,CSTATE> sol_mat(3,1,0.);
        for(auto is = 0; is < nsol; is++){
          const auto coeff = m_coeff[is];
          //only ZERO
          if(coeff == 0.){continue;}
          const auto beta = m_beta[is];
          for(auto x = 0; x < 3; x++){
            const auto val = sol_mat.Get(x,0);
            const auto src = et_beta.Get(x,is);
            sol_mat.Put(x, 0, val + src*coeff);
          }
        }
        coeff_mat.Substitution(&sol_mat);
        TPZFMatrix<CSTATE> fmat(nsol,1,this->m_f_scratch[index].begin(),nsol);
        fmat.AddContribution(0, 0, rot_et, conj_trans, sol_mat, no_trans, cte*2.0);
      }
    }else{
      const STATE sign = m_pos_z ? 1 : -1;
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

      TPZFNMatrix<3000,CSTATE> rot_et(3,nsol,0.), et_beta(3,nsol,0.);

      TPZFNMatrix<3000,CSTATE> rot_grad_ez(3,nsol,0.);
      TPZFNMatrix<2000,CSTATE> grad_ez_axes(2,nsol,0.);

      for(auto isol = 0; isol < nsol; isol++){
        const auto &et_ref = datavec[ TPZWgma::HCurlIndex() ].sol[isol];
        const auto &ex = et_ref[0];
        const auto &ey = et_ref[1];
        const auto beta = m_beta[isol];
        rot_et.Put(0,isol,ey);
        rot_et.Put(1,isol,-ex);
        et_beta.Put(0,isol,1i*beta*ey);
        et_beta.Put(1,isol,-1i*beta*ex);
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

      TPZFNMatrix<3000,CSTATE> tmp(et_beta);
      if(m_pos_z){
        tmp += rot_grad_ez;
      }else{
        tmp -= rot_grad_ez;
      }
      ur.Substitution(&tmp);
      this->m_k_scratch[index].AddContribution(0, 0, rot_et, conj_trans, tmp, no_trans,cte);
      //src term
      if(is_src){
        //compute solution
        TPZFNMatrix<3,CSTATE> sol_mat(3,1,0.);
        for(auto is = 0; is < nsol; is++){
          const auto coeff = m_coeff[is];
          //only ZERO
          if(coeff == 0.){continue;}
          const auto beta = m_beta[is];
          for(auto x = 0; x < 3; x++){
            const auto val = sol_mat.Get(x,0);
            const auto src = et_beta.Get(x,is)+sign*rot_grad_ez.Get(x,is);
            sol_mat.Put(x, 0, val + src*coeff);
          }
        }
        ur.Substitution(&sol_mat);
        TPZFMatrix<CSTATE> fmat(nsol,1,this->m_f_scratch[index].begin(),nsol);
        fmat.AddContribution(0, 0, rot_et, conj_trans, sol_mat, no_trans, cte*2.0);
      }
    }
  }

  template
  class WaveguidePortBC<SingleSpaceIntegrator>;
  template
  class WaveguidePortBC<MultiphysicsIntegrator>;
  
};