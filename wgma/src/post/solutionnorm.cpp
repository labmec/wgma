#include "post/solutionnorm.hpp"
#include <TPZBndCond.h>
#include <TPZMaterialDataT.h>
#include <pzvec_extras.h>

namespace wgma::post{

  template<class TSPACE>
  TPZVec<CSTATE> SolutionNorm<TSPACE>::ComputeNorm()
  {
    auto mesh = this->Mesh();
    const int size_res = std::max(this->NThreads(),1);
    const int nsol = mesh->Solution().Cols();
    m_res.Resize(size_res);
    for (auto &it : m_res){
      it.Resize(nsol);
      it.Fill(0);
    }
    this->Integrate(this->m_elvec);
    
    TPZVec<CSTATE> res(nsol,0.);
    for (auto &it : m_res){
      for(int isol = 0; isol < nsol; isol++){
        res[isol] += it[isol];
      }
    }

    for(int isol = 0; isol < nsol; isol++){
      res[isol] = sqrt(res[isol]);
      // std::cout<<"computed norm of solution "<<isol<<": "<<res[isol]<<std::endl;
    }
    m_res.Resize(0);
    return res;
  }
  
  template<class TSPACE>
  TPZVec<CSTATE> SolutionNorm<TSPACE>::Normalise()
  {
    auto res = ComputeNorm();
    
    auto mesh = this->Mesh();
    TPZFMatrix<CSTATE> &evectors = mesh->Solution();
    const int nev = evectors.Cols();
    const int neq = mesh->NEquations();
    //we iterate through the eigenvectors
    for(int iev = 0; iev < nev; iev++){
      const auto norm = res[iev];
      const int offset = iev * neq;
      TPZFMatrix<CSTATE> ei(neq,1,evectors.Elem() + offset,neq);
      ei *= 1./norm;
    }
    return res;
  }

  template<class TSPACE>
  CSTATE SolutionNorm<TSPACE>::ComputeNorm(int s){
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
  void SolutionNorm<TSPACE>::Compute(const ElData &eldata, REAL weight, int index)
  {
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;
      
      const int nsol = data.sol.size();
      for(int isol = 0; isol < nsol; isol++){
        const auto &cursol =  data.sol[isol];
        const auto solsize = cursol.size();
        CSTATE val = 0;
        for(auto ix = 0; ix < solsize; ix++){
          val += m_conj ? cursol[ix] * std::conj(cursol[ix])
            : cursol[ix] * cursol[ix];
        }
        this->m_res[index][isol] += weight * fabs(data.detjac) * val;
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
  class SolutionNorm<SingleSpaceIntegrator>;
  template
  class SolutionNorm<MultiphysicsIntegrator>;
  
};