#include "post/evalsolution.hpp"
#include <pzcondensedcompel.h>
using namespace std::complex_literals;

namespace wgma::post{

  template<class TSPACE>
  void EvalSolution<TSPACE>::EvalSolutionAtPoints(STATE &min, STATE &max)
  {

    auto mesh = this->Mesh();
    const int size_res = std::max(this->NThreads(),1);
    //first column is the mode
    const int nsol = mesh->Solution().Cols();
    if(nsol != 1){DebugStop();}
    
    m_min_sol.Resize(size_res); m_min_sol.Fill(1e20);
    m_max_sol.Resize(size_res); m_max_sol.Fill(0);
    this->Integrate(this->m_elvec);

    min = std::sqrt(*std::min_element(m_min_sol.begin(),m_min_sol.end()));
    max = std::sqrt(*std::max_element(m_max_sol.begin(),m_max_sol.end()));

    
  }


  template<class TSPACE>
  void EvalSolution<TSPACE>::Compute(const ElData &eldata, REAL weight, int index)
  {
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;
      
      const auto &sol = data.sol[0];
      
      STATE sol_norm = 0;
      //numerator
      CSTATE val_num = 0;
      const int solsize = sol.size();
      for(auto ix = 0; ix < solsize; ix++){
        sol_norm +=(sol[ix]*std::conj(sol[ix])).real();
      }
      m_min_sol[index] = std::min(sol_norm,m_min_sol[index]);
      m_max_sol[index] = std::max(sol_norm,m_max_sol[index]);
    }else{
      DebugStop();
    }
  }

  template
  class EvalSolution<SingleSpaceIntegrator>;
  template
  class EvalSolution<MultiphysicsIntegrator>;
};