#include "post/solutionnorm.hpp"
#include <TPZBndCond.h>
#include <TPZMaterialDataT.h>
#include <pzvec_extras.h>

namespace wgma::post{  
  STATE SolutionNorm::ComputeNorm(int s){
    const int size_res = std::max(NThreads(),1);
    m_res.Resize(size_res, 0);
    this->SetSol(s);
    this->Integrate(m_elvec);
    STATE res = 0;
    for (auto &it : m_res){
      res += it;
    }
    return res;
  }
  
  void SolutionNorm::Compute(const ElData &eldata, REAL weight, int index)
  {
    const TPZMaterialDataT<CSTATE> &data = eldata;
    const auto &which = WhichSol();
    const auto &cursol = data.sol[which];
    const auto solsize = cursol.size();
    for(auto ix = 0; ix < solsize; ix++){
      const auto val = std::real(cursol[ix] * std::conj(cursol[ix]));
      m_res[index] += weight * fabs(data.detjac) * val;
    }
  }
  
};