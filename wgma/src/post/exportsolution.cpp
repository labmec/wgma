#include "post/exportsolution.hpp"
#include <pzcondensedcompel.h>
using namespace std::complex_literals;

namespace wgma::post{

  template<class TSPACE>
  void ExportSolution<TSPACE>::StoreSolutionAtPoints(const int sol_dim)
  {

    auto GetElPts = [](auto el){
      auto condensed = dynamic_cast<TPZCondensedCompEl*>(el);
      if(condensed){
        return condensed->ReferenceCompEl()->GetIntegrationRule().NPoints();
      }
      return el->GetIntegrationRule().NPoints();
    };
    m_sol_dim = sol_dim;
    const auto nel = this->m_elvec.size();
    m_first_pt_el.Resize(nel);
    m_first_pt_el[0] = 0;
    int64_t total_pts{0};
    for(auto iel = 0; iel < nel-1; iel++){
      auto el = this->m_elvec[iel];
      const int npts = GetElPts(el);
      total_pts += npts;
      m_first_pt_el[iel+1] = total_pts;
    }
    //now we get the remaining pts for the last element
    total_pts += GetElPts(this->m_elvec[nel-1]);
    m_weights.Resize(total_pts);
    m_sol.Resize(total_pts*m_sol_dim);
    this->Integrate(this->m_elvec);
  }


  template<class TSPACE>
  void ExportSolution<TSPACE>::Compute(const ElData &eldata, REAL weight, int index)
  {
    const auto el_idx = eldata.GetElIndex();
    const auto first_pt = m_first_pt_el[el_idx];
    int64_t loc_id{-1};
    REAL detjac{0};

    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;
      loc_id = data.intLocPtIndex;
      detjac = data.detjac;
    }else{
      const TPZVec<TPZMaterialDataT<CSTATE>> &datavec = eldata;
      loc_id = datavec[0].intLocPtIndex;
      detjac = datavec[0].detjac;
    }

    m_weights[first_pt+loc_id] = weight * fabs(detjac);
    
    int64_t pt = (first_pt+loc_id)*m_sol_dim;

    auto WriteSol = [this, &pt](const auto data){
      const auto nsol = data.sol.size();
      for(auto isol = 0; isol < nsol; isol++){
        const auto &sol = data.sol[isol];
        const auto solsz = sol.size();
        for(auto x = 0; x < solsz; x++){
          m_sol[pt++] = sol[x];
        }
      }
    };
    
    if constexpr(std::is_same_v<TSPACE,SingleSpaceIntegrator>){
      const TPZMaterialDataT<CSTATE> &data = eldata;
      WriteSol(data);
    }else{
      const TPZVec<TPZMaterialDataT<CSTATE>> &datavec = eldata;
      const auto ndata = datavec.size();
      for(auto idata = 0; idata < ndata; idata++){
        const TPZMaterialDataT<CSTATE> &data = datavec[idata];
        WriteSol(data);
      }
    }
  }

  template
  class ExportSolution<SingleSpaceIntegrator>;
  template
  class ExportSolution<MultiphysicsIntegrator>;
};