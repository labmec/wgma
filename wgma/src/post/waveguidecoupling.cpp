#include "post/waveguidecoupling.hpp"
#include <TPZBndCond.h>
#include <TPZMaterialDataT.h>
#include <Electromagnetics/TPZScalarField.h>
#include <Electromagnetics/TPZWgma.h>
#include <TPZElectromagneticConstants.h>
#include <pzvec_extras.h>
#include <pzaxestools.h>
using namespace std::complex_literals;

/*
  According to
  Operator Theory for Electromagnetics: An Introduction, G.W. Hanson,
  Theorem 4.30,

  given x_n as the eigenvectors of
  Ax = lBx
  and y_m as the eigenvectors of the adjoint problem
  we have
  (l_n-l_m)(Bx_n,y_m) = 0

  in case one wants to check this orthogonality,
  uncomment next line

  p.s. worth mentioning that, according to
  Quasimodal expansion of electromagnetic fields in open two-dimensional structures

  B. Vial, F. Zolla, A.Nicolet and M. Commandré,

  for homogeneous BCs, y_m = x_m^*
 */
// #define CHECK_ORTH

namespace wgma::post{

  template<class TSPACE>
  void WaveguideCoupling<TSPACE>::ComputeCoupling(){
    const int size_res = std::max(this->NThreads(),1);
    const int nsol = this->Mesh()->Solution().Cols();
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
    
    auto wgdata = dynamic_cast<WgCouplData*>(&data);
    if(!wgdata){DebugStop();}
    wgdata->m_elindex = el->Index();
    const auto nsol = el->Mesh()->Solution().Cols();
    if(m_print_mats){wgdata->m_elmat.Resize(nsol,nsol);}
  }


  template<class TSPACE>
  void WaveguideCoupling<TSPACE>::PostProcessData(ElData& data)
  {
    if(!m_print_mats) return;
    auto wgdata = dynamic_cast<WgCouplData*>(&data);
    
    std::ofstream matfile{m_prefix+"mat_"+std::to_string(wgdata->m_elindex)+".csv"};
    wgdata->m_elmat.Print("",matfile,ECSV);
  }
  
  template<class TSPACE>
  void WaveguideCoupling<TSPACE>::Compute(const ElData &eldataconst, REAL weight, int index)
  {
    ElData& eldata = (ElData&)eldataconst;
    auto wgdata = dynamic_cast<WgCouplData*>(&eldata);
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;

      auto mat = dynamic_cast<const TPZScalarField*>(eldata.GetMaterial());
      TPZFNMatrix<9,CSTATE> constitutive_mat;
      if(m_te){
        mat->GetPermeability(data.x, constitutive_mat);
      }else{
        mat->GetPermittivity(data.x, constitutive_mat);
      }
      constitutive_mat.Decompose(ELU);

      const int nsol = data.sol.size();
      const auto detjac = data.detjac;
      const CSTATE cte = weight*fabs(detjac);

      TPZFNMatrix<3000,CSTATE> et(3,nsol,0.);

      for(auto isol = 0; isol < nsol; isol++){
        const auto &et_ref = data.sol[isol];
        et.Put(1,isol,et_ref[0]);
      }
      TPZFNMatrix<3000,CSTATE> tmp;
      tmp = et;
      constitutive_mat.Substitution(&tmp);
      //whether to conjugate or only tanspose
      const int transp = m_conj ? 2 : 1;
      this->m_k_scratch[index].AddContribution(0, 0, et, transp, tmp, false,cte);
    }else{
      //we expect a TPZWgma material
      const TPZVec<TPZMaterialDataT<CSTATE>> &datavec = eldata;

      auto mat = dynamic_cast<const TPZWgma*>(eldata.GetMaterial());
      TPZFNMatrix<9,CSTATE> ur;
      mat->GetPermeability(datavec[0].x, ur);
      ur.Decompose(ELU);
#ifdef CHECK_ORTH
      TPZFNMatrix<9,CSTATE> er;
      mat->GetPermittivity(datavec[0].x, er);
#endif

      const int nsol = datavec[0].sol.size();
      const auto &axes = datavec[0].axes;
      const auto detjac = datavec[0].detjac;
      const CSTATE cte = weight*fabs(detjac);

      TPZFNMatrix<3000,CSTATE> rot_et(3,nsol,0.);

      
#ifdef CHECK_ORTH
      TPZFNMatrix<3000,CSTATE> ez(3,nsol,0.);
#endif
      TPZFNMatrix<3000,CSTATE> rot_grad_ez(3,nsol,0.);
      TPZFNMatrix<3000,CSTATE> grad_ez_axes(2,nsol,0.);

      for(auto isol = 0; isol < nsol; isol++){
        const auto &et_ref = datavec[ TPZWgma::HCurlIndex() ].sol[isol];
        rot_et.Put(0,isol,et_ref[1]);
        rot_et.Put(1,isol,-et_ref[0]);
#ifdef CHECK_ORTH
        const auto &ez_ref = datavec[ TPZWgma::H1Index() ].sol[isol];
        ez.Put(2,isol,ez_ref[0]);
#endif
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

      //Btt term
      tmp = rot_et;
      tmp += rot_grad_ez;
      if(m_conj){
        const auto sz = tmp.Rows()*tmp.Cols();
        CSTATE *tmp_ptr = tmp.Elem();
        for(int i = 0; i < sz; i++,tmp_ptr++){
          *tmp_ptr = std::conj(*tmp_ptr);
        }
      }
      ur.Substitution(&tmp);
      tmp *= cte;
      this->m_k_scratch[index].AddContribution(0, 0, rot_et, true, tmp, false);


#ifdef CHECK_ORTH
      //Azt term
      tmp = rot_grad_ez;
      ur.Substitution(&tmp);
      tmp *= cte;
      this->m_k_scratch[index].AddContribution(0,0,tmp, true, rot_et, false);
      //Azz term
      this->m_k_scratch[index].AddContribution(0,0,tmp,true,rot_grad_ez,false);
      //Czz term
      er.Multiply(ez,tmp);
      const auto wl = mat->GetWavelength();
      const auto sc = mat->GetScaleFactor();
      const auto k0 = sc*2*M_PI/wl;
      tmp *= -k0*k0*cte;
      this->m_k_scratch[index].AddContribution(0,0,tmp,true,ez,false);
#endif
    }
  }

  template
  class WaveguideCoupling<SingleSpaceIntegrator>;
  template
  class WaveguideCoupling<MultiphysicsIntegrator>;
  
};
