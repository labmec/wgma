#ifndef _ELDATA_HPP_
#define _ELDATA_HPP_

#include <TPZMaterialDataT.h>
#include <TPZMaterial.h>
namespace wgma::post{

  //!Contains data for computing the solution of a given element
  class ElData {
  public:
    ElData() = default;
    virtual ~ElData() = default;
    ElData(const ElData&) = default;
    ElData(ElData&&) = default;
    ElData& operator=(const ElData&) = default;
    ElData& operator=(ElData&&) = default;
    ElData(TPZMaterialDataT<CSTATE>& ar) : m_data(ar){};
    ElData(TPZVec<TPZMaterialDataT<CSTATE>>& ar) : m_datavec(ar){};
    operator TPZMaterialDataT<CSTATE>&() {return m_data;}
    operator TPZVec<TPZMaterialDataT<CSTATE>>&() {return m_datavec;}
    operator const TPZMaterialDataT<CSTATE>&() const{return m_data;}
    operator const TPZVec<TPZMaterialDataT<CSTATE>>&() const {return m_datavec;}
    const TPZMaterial *GetMaterial() const {return m_mat;}
    void SetMaterial(TPZMaterial *mat){m_mat = mat;}
  protected: 
    TPZMaterialDataT<CSTATE> m_data;
    TPZManVector<TPZMaterialDataT<CSTATE>,10> m_datavec;
    TPZMaterial * m_mat{nullptr};
  };
};

#endif /* _ELDATA_HPP_ */
