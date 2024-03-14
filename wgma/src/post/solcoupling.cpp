#include "post/solcoupling.hpp"

#include <TPZMaterialDataT.h>
#include <pzvec_extras.h>
#include <pzaxestools.h>


namespace wgma::post{

  template<class TSPACE>
  void SolCoupling<TSPACE>::ComputeCoupling(){
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
  void SolCoupling<TSPACE>::InitData(TPZCompEl *el, ElData &data)
  {
    TSPACE::InitData(el,data);
    const auto nsol = el->Mesh()->Solution().Cols();
    auto wgdata = dynamic_cast<SolCouplData*>(&data);
    if(!wgdata){DebugStop();}
    wgdata->m_elindex = el->Index();
    if(m_print_mats){wgdata->m_elmat.Resize(nsol,nsol);}
  }


  template<class TSPACE>
  void SolCoupling<TSPACE>::PostProcessData(ElData& data)
  {
    if(!m_print_mats) return;
    auto wgdata = dynamic_cast<SolCouplData*>(&data);
    std::ofstream matfile{m_prefix+"mat_"+std::to_string(wgdata->m_elindex)+".csv"};
    wgdata->m_elmat.Print("",matfile,ECSV);
  }
  
  template<class TSPACE>
  void SolCoupling<TSPACE>::Compute(const ElData &eldataconst, REAL weight, int index)
  {
    ElData& eldata = (ElData&)eldataconst;
    auto wgdata = dynamic_cast<SolCouplData*>(&eldata);
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;

      const int nsol = data.sol.size();
      const auto detjac = data.detjac;
      const CSTATE cte = weight*fabs(detjac);

      TPZFNMatrix<1000,CSTATE> et(1,nsol,0.),et_conj(1,nsol,0.);

      for(auto isol = 0; isol < nsol; isol++){
        const auto &et_ref = data.sol[isol];
        et.Put(0,isol,et_ref[0]);
        et_conj.Put(0,isol,std::conj(et_ref[0]));
      }
      this->m_k_scratch[index].AddContribution(0, 0, et_conj, true, et, false,cte);
    }else{
      DebugStop();
    }
  }

  template
  class SolCoupling<SingleSpaceIntegrator>;
  template
  class SolCoupling<MultiphysicsIntegrator>;
  
};
