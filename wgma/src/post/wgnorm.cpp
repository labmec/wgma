#include "post/wgnorm.hpp"
#include <TPZBndCond.h>
#include <TPZMaterialDataT.h>
#include <Electromagnetics/TPZScalarField.h>
#include <Electromagnetics/TPZWgma.h>
#include <pzvec_extras.h>
#include <pzaxestools.h>
using namespace std::complex_literals;

namespace wgma::post{

  template<class TSPACE>
  TPZVec<CSTATE> WgNorm<TSPACE>::ComputeNorm()
  {
    auto mesh = this->Mesh();
    const int size_res = std::max(this->NThreads(),1);
    const int nsol = mesh->Solution().Cols();
    m_res.Resize(size_res,TPZVec<CSTATE>(nsol,0.));
    this->Integrate(this->m_elvec);
    
    TPZVec<CSTATE> res(nsol,0.);
    for (auto &it : m_res){
      for(int isol = 0; isol < nsol; isol++){
        res[isol] += it[isol];
      }
    }

    return res;
  }
  
  template<class TSPACE>
  TPZVec<CSTATE> WgNorm<TSPACE>::Normalise()
  {
    auto res = ComputeNorm();
    
    auto mesh = this->Mesh();
    TPZFMatrix<CSTATE> &evectors = mesh->Solution();
    const int nev = evectors.Cols();
    const int neq = mesh->NEquations();
    //we iterate through the eigenvectors
    for(int iev = 0; iev < nev; iev++){
      const auto norm = std::sqrt(res[iev]);
      const int offset = iev * neq;
      TPZFMatrix<CSTATE> ei(neq,1,evectors.Elem() + offset,neq);
      //let us avoid nasty divisions
      if(std::abs(norm) > 1e-12){ei *= M_SQRT1_2/norm;}
    }
    return res;
  }

  template<class TSPACE>
  CSTATE WgNorm<TSPACE>::ComputeNorm(int s){
    //we store the solution
    TPZFMatrix<CSTATE> solcp = this->Mesh()->Solution();
    const auto neq = solcp.Rows();
    const auto nsols = solcp.Cols();
    if(s >= nsols){
      DebugStop();
    }
    //we set only desired sol in mesh
    TPZFMatrix<CSTATE> desired_sol(neq,1,&solcp.s(0,s),neq);
    this->Mesh()->LoadSolution(desired_sol);
    //normalise it
    auto res = this->Normalise();
    //restore previous solution
    this->Mesh()->LoadSolution(std::move(solcp));
    return res[0];
  }

  template<class TSPACE>
  void WgNorm<TSPACE>::Compute(const ElData &eldata, REAL weight, int index)
  {
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
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
      const auto detjac = data.detjac;
      //speed of light times \mu_0
      const auto c_uo = 120*M_PI;
      const auto omega_uo = (2*M_PI/m_wl)*c_uo;
      const auto cte = weight*fabs(detjac);
      const int nsol = data.sol.size();
      TPZFNMatrix<3000,CSTATE> rot_e_field(3,nsol,0.), h_field(3,nsol,0.);
      for(int isol = 0; isol < nsol; isol++){
        const auto val = data.sol[isol][0];
        rot_e_field.Put(1,isol,-val);
        const auto beta = m_beta[isol];
        //we have -jbeta/(j\omega\mu_0)
        const auto ct = -1.0*beta/omega_uo;
        h_field.Put(1,isol,ct*val);
      }
      //apply constitutive param
      coeff_mat.Substitution(&h_field);
      //if conj
      if(m_conj){
        for(int isol = 0; isol < nsol; isol++){
          const auto val = h_field.Get(1,isol);
          h_field.Put(1,isol,std::conj(val));
        }
      }
      
      for(int isol = 0; isol < nsol; isol++){
        const auto e = rot_e_field.Get(1,isol);
        const auto h = h_field.Get(1,isol);
        
        this->m_res[index][isol] += e * h * cte;
      }
    }else{
      const TPZVec<TPZMaterialDataT<CSTATE>> &datavec = eldata;

      const auto nsol = datavec[0].sol.size();
      TPZFNMatrix<3000,CSTATE> rot_et(3,nsol,0.), rot_et_beta(3,nsol,0.),rot_grad_ez(3,nsol,0.);
      TPZFNMatrix<2000,CSTATE> grad_ez_axes(2,nsol,0.);

      for(auto isol = 0; isol < nsol; isol++){
        const auto &et_ref = datavec[ TPZWgma::HCurlIndex() ].sol[isol];
        const auto &ex = et_ref[0];
        const auto &ey = et_ref[1];
        const auto beta = m_beta[isol];
        rot_et.Put(0,isol,ey);
        rot_et.Put(1,isol,-ex);
        rot_et_beta.Put(0,isol,1i*beta*ey);
        rot_et_beta.Put(1,isol,-1i*beta*ex);
        auto &gradez_ref = datavec[ TPZWgma::H1Index() ].dsol[isol];
        grad_ez_axes.Put(0,isol,gradez_ref[0]);
        grad_ez_axes.Put(1,isol,gradez_ref[1]);
      }

      const auto &axes = datavec[0].axes;
      TPZAxesTools<CSTATE>::Axes2XYZ(grad_ez_axes, rot_grad_ez, axes);

      TPZFNMatrix<3000,CSTATE> h_field;
      h_field = rot_et_beta;
      h_field+= rot_grad_ez;

      auto mat = dynamic_cast<const TPZWgma*>(eldata.GetMaterial());
      TPZFNMatrix<9,CSTATE> ur;
      mat->GetPermeability(datavec[0].x, ur);
      ur.Decompose(ELU);
      ur.Substitution(&h_field);

      const auto c_uo = 120*M_PI;
      const auto omega_uo = (2*M_PI/m_wl)*c_uo;
      h_field *= 1i/omega_uo;

      TPZFMatrix<CSTATE> solmat(nsol,1,this->m_res[index].begin(),nsol);
      
      constexpr int no_trans{0};
      constexpr int conj_trans{2};
      const auto detjac = datavec[0].detjac;
      const CSTATE cte = weight*fabs(detjac);
      for(int is = 0; is < nsol; is++){
        const auto offset = 3*is;
        TPZFMatrix<CSTATE> my_e(3,1,rot_et.Elem()+offset,3);
        TPZFMatrix<CSTATE> my_h(3,1,h_field.Elem()+offset,3);
        solmat.AddContribution(is,0,my_h,conj_trans,my_e,no_trans,cte);
      }
    }
  }

  template
  class WgNorm<SingleSpaceIntegrator>;
  template
  class WgNorm<MultiphysicsIntegrator>;
  
};