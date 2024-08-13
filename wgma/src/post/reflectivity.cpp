#include "post/reflectivity.hpp"
#include <TPZMaterialDataT.h>
#include <pzvec_extras.h>
#include <numeric>
namespace wgma::post{

  template<class TSPACE>
  CSTATE SolutionReflectivity<TSPACE>::ComputeReflectivity()
  {
    auto mesh = this->Mesh();
    const int size_res = std::max(this->NThreads(),1);
    //first column is the mode
    const int nsol = mesh->Solution().Cols();
    if(nsol != 2){DebugStop();}
    
    m_numerator.Resize(size_res); m_numerator.Fill(0);
    m_denominator.Resize(size_res); m_numerator.Fill(0);
    this->Integrate(this->m_elvec);
    
    CSTATE num{0}, den{0};
    for(auto it = 0; it < size_res; it++){
      num += m_numerator[it];
      den += m_denominator[it];
    }
    m_numerator.Resize(0);
    m_denominator.Resize(0);
    return num/den;
  }


  template<class TSPACE>
  void SolutionReflectivity<TSPACE>::Compute(const ElData &eldata, REAL weight, int index)
  {
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;
      
      const auto &sol = data.sol[0];
      const auto &src =  data.sol[1];
      //denominator
      CSTATE val_den = 0;
      //numerator
      CSTATE val_num = 0;
      const int solsize = sol.size();
      for(auto ix = 0; ix < solsize; ix++){
        val_num +=(sol[ix]-src[ix])*std::conj(src[ix]);
        val_den +=src[ix]*std::conj(src[ix]);
      }
      this->m_numerator[index] += weight * fabs(data.detjac) * val_num;
      this->m_denominator[index] += weight * fabs(data.detjac) * val_den;
    }else{
      DebugStop();
    }
  }

  template
  class SolutionReflectivity<SingleSpaceIntegrator>;
  template
  class SolutionReflectivity<MultiphysicsIntegrator>;
  
};