#include "post/solutionnorm.hpp"
#include <TPZBndCond.h>
#include <TPZMaterialDataT.h>
#include <pzvec_extras.h>

namespace wgma::post{

  template<class TSPACE>
  void SolutionNorm<TSPACE>::Normalise()
  {

    auto mesh = this->Mesh();

    TPZFMatrix<CSTATE> &evectors = mesh->Solution();
    const int nev = evectors.Cols();
    const int neq = mesh->NEquations();
    //we iterate through the eigenvectors
    for(int iev = 0; iev < nev; iev++){
      const auto norm = this->ComputeNorm(iev);
      const int offset = iev * neq;
      TPZFMatrix<CSTATE> ei(neq,1,evectors.Elem() + offset,neq);
      std::cout<<"computed norm of solution "<<iev<<": "<<norm<<std::endl;
      ei *= 1./norm;
    }
  }

  template<class TSPACE>
  STATE SolutionNorm<TSPACE>::ComputeNorm(int s){
    const int size_res = std::max(this->NThreads(),1);
    m_res.Resize(size_res);
    m_res.Fill(0.);
    this->SetSol(s);
    this->Integrate(this->m_elvec);
    STATE res = 0;
    for (auto &it : m_res){
      res += it;
    }
    return sqrt(res);
  }

  template<class TSPACE>
  void SolutionNorm<TSPACE>::Compute(const ElData &eldata, REAL weight, int index)
  {
    const auto which = WhichSol();
    
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;
      const auto &cursol =  data.sol[which];
      const auto solsize = cursol.size();
      STATE val = 0;
      for(auto ix = 0; ix < solsize; ix++){
         val += std::real(cursol[ix] * std::conj(cursol[ix]));
      }
      this->m_res[index] += weight * fabs(data.detjac) * val;
    }else{
      const TPZVec<TPZMaterialDataT<CSTATE>> &datavec = eldata;
      const auto nspace = datavec.size();
      for(int is = 0; is < nspace; is++){
        const auto &data = datavec[is];
        const auto &cursol =  data.sol[which];
        const auto solsize = cursol.size();
        STATE val = 0;
        for(auto ix = 0; ix < solsize; ix++){
          val += std::real(cursol[ix] * std::conj(cursol[ix]));
        }
        this->m_res[index] += weight * fabs(data.detjac) * val;
      }
    }
  }

  template
  class SolutionNorm<SingleSpaceIntegrator>;
  template
  class SolutionNorm<MultiphysicsIntegrator>;
  
};