#ifndef _ELDATA_HPP_
#define _ELDATA_HPP_

#include <TPZMaterialDataT.h>

namespace wgma::post{

  //!Contains data for computing the solution of a given element
  class ElData {
  public:
    ElData() = default;
    ElData(TPZMaterialDataT<CSTATE>& ar) : m_data(ar){};
    ElData(TPZVec<TPZMaterialDataT<CSTATE>>& ar) : m_datavec(ar){};
    operator TPZMaterialDataT<CSTATE>&() {return m_data;}
    operator TPZVec<TPZMaterialDataT<CSTATE>>&() {return m_datavec;}
    operator const TPZMaterialDataT<CSTATE>&() const{return m_data;}
    operator const TPZVec<TPZMaterialDataT<CSTATE>>&() const {return m_datavec;}
  protected: 
    TPZMaterialDataT<CSTATE> m_data;
    TPZManVector<TPZMaterialDataT<CSTATE>,10> m_datavec;
  };
};

#endif /* _ELDATA_HPP_ */
