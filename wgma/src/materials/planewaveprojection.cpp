#include "materials/planewaveprojection.hpp"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"

using namespace wgma::materials;

template<class TVar>
PlanewaveProjection<TVar>::PlanewaveProjection(int id, REAL theta) :
  TPZRegisterClassId(&PlanewaveProjection::ClassId),
  TBase(id), m_theta(theta)
{
}

template<class TVar>
TPZMaterial * PlanewaveProjection<TVar>::NewMaterial() const{
	return new PlanewaveProjection(*this);
}

template<class TVar>
void PlanewaveProjection<TVar>::Contribute(const TPZMaterialDataT<TVar> &data,
                                          REAL weight,
                                          TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef){
  //nshape
  const int nr = data.phi.Rows();
  //shape dim
  const int nc = data.phi.Cols();
  
  const TPZFNMatrix<300,TVar> phi = [&data, &nc, &nr]() ->TPZFMatrix<TVar>{
    if constexpr (std::is_same_v<TVar,STATE>){
      return data.phi;
    }else{
      TPZFNMatrix<300,TVar> phi_cplx(nr,nc,0);
      for(int i = 0; i < nr; i++){
        for(int j = 0; j < nc; j++){
          phi_cplx.Put(i,j,data.phi.Get(i,j));
        }
      }
      return phi_cplx;
    }
  }();

  TPZFNMatrix<3,TVar> sol(3,1,0);
  sol.PutVal(0,0,cos(m_theta));
  sol.PutVal(1,0,sin(m_theta));
  sol.PutVal(2,0,0);
  
  ek.AddContribution(0, 0, phi, false, phi, true,weight*m_scale);
  ef.AddContribution(0, 0, phi, false, sol, false, weight*m_scale);
}

template<class TVar>
int PlanewaveProjection<TVar>::VariableIndex(const std::string &name) const{
	if(!strcmp("Solution",name.c_str())){
    return ESol;
  }
  if(!strcmp("Derivative",name.c_str())) return EDerivative;
	return TPZMaterial::VariableIndex(name);
}

template<class TVar>
int PlanewaveProjection<TVar>::NSolutionVariables(int var) const{
	if(var == ESol) {return 3;}
  if(var == EDerivative) {return 1;}
  return TPZMaterial::NSolutionVariables(var);
}

template<class TVar>
void PlanewaveProjection<TVar>::Solution(const TPZMaterialDataT<TVar> &data,
                                        int var, TPZVec<TVar> &solOut)
{
  const auto &sol = data.sol[0];
	if(var == ESol){
    solOut.Resize(sol.size());
    for (int i=0; i<sol.size(); i++) {
      solOut[i] = std::real(sol[i])/m_scale;
    }
		return;
	}
  if(var == EDerivative){
    const auto &dsol = data.dsol[0];
    solOut.Resize(fDim);
    for (int i=0; i<fDim; i++) {
      solOut[i] = std::real(dsol.GetVal(i,0))/m_scale;
    }
    return;
  }
  DebugStop();
}

template<class TVar>
int PlanewaveProjection<TVar>::ClassId() const{
  return Hash("PlanewaveProjection") ^ TBase::ClassId() << 1;
}


template class wgma::materials::PlanewaveProjection<STATE>;
template class wgma::materials::PlanewaveProjection<CSTATE>;
