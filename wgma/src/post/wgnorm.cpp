#include "post/wgnorm.hpp"
#include <TPZBndCond.h>
#include <TPZMaterialDataT.h>
#include <Electromagnetics/TPZScalarField.h>
#include <pzvec_extras.h>
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
      const auto cte = 0.5*weight*fabs(detjac);
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
      const auto nspace = datavec.size();
      for(int is = 0; is < nspace; is++){
        const auto &data = datavec[is];
        const int nsol = data.sol.size();
        for(int isol = 0; isol < nsol; isol++){
          const auto &cursol =  data.sol[isol];
          const auto solsize = cursol.size();
          CSTATE val = 0;
          for(auto ix = 0; ix < solsize; ix++){
            val += m_conj ?
              cursol[ix] * std::conj(cursol[ix]) : cursol[ix] * cursol[ix];
          }
          this->m_res[index][isol] += weight * fabs(data.detjac) * val;
        }
      }
    }
  }

  template
  class WgNorm<SingleSpaceIntegrator>;
  template
  class WgNorm<MultiphysicsIntegrator>;
  
};