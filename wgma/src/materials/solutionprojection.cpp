#include "materials/solutionprojection.hpp"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"

using namespace wgma::materials;

template<class TVar>
SolutionProjection<TVar>::SolutionProjection(int id, int dim) :
  TPZRegisterClassId(&SolutionProjection::ClassId),
  TBase(id), fDim(dim)
{
}

template<class TVar>
TPZMaterial * SolutionProjection<TVar>::NewMaterial() const{
	return new SolutionProjection(*this);
}

template<class TVar>
void SolutionProjection<TVar>::FillDataRequirements(TPZMaterialData &data) const
{
  data.fNeedsSol = true;
}

template<class TVar>
void SolutionProjection<TVar>::Contribute(const TPZMaterialDataT<TVar> &data,
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

  const int sz = data.sol[0].size();
  if(sz!=nc){
    DebugStop();
  }
  TPZFNMatrix<3,TVar> sol(sz,1,0);
  for(int i = 0; i < sz; i++){ sol.Put(i,0,data.sol[0][i]);}
  
  ek.AddContribution(0, 0, phi, false, phi, true,weight*fScale);
  ef.AddContribution(0, 0, phi, false, sol, false, weight*fScale);
}

template<class TVar>
void SolutionProjection<TVar>::ContributeBC(const TPZMaterialDataT<TVar> &data,
                                            REAL weight,
                                            TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                                            TPZBndCondT<TVar> &bc)
{
  return;
}
template<class TVar>
void SolutionProjection<TVar>::GetSolDimensions(uint64_t &u_len,
                                                uint64_t &du_row,
                                                uint64_t &du_col) const
{
  u_len=1;
  du_row=3;
  du_col=1;
}


template<class TVar>
int SolutionProjection<TVar>::VariableIndex(const std::string &name) const{
	if(!strcmp("Solution",name.c_str())) return ESolution;
  if(!strcmp("Derivative",name.c_str())) return EDerivative;
	return TPZMaterial::VariableIndex(name);
}

template<class TVar>
int SolutionProjection<TVar>::NSolutionVariables(int var) const{
	if(var == ESolution) return 1;
  if (var == EDerivative) {
    return fDim;
  }
	
  return TPZMaterial::NSolutionVariables(var);
}

template<class TVar>
void SolutionProjection<TVar>::Solution(const TPZMaterialDataT<TVar> &data,
                                        int var, TPZVec<TVar> &solOut)
{
  const auto &sol = data.sol[0];
  const auto &dsol = data.dsol[0];
	if (var == ESolution){
    solOut.Resize(sol.size());
    for (int i=0; i<sol.size(); i++) {
      solOut[i] = sol[i]/fScale;
    }
		return;
	}
  if (var == EDerivative) {
    solOut.Resize(fDim);
    for (int i=0; i<fDim; i++) {
      solOut[i] = dsol.GetVal(i,0)/fScale;
    }
    return;
  }
}

template<class TVar>
int SolutionProjection<TVar>::ClassId() const{
  return Hash("SolutionProjection") ^ TBase::ClassId() << 1;
}


template class wgma::materials::SolutionProjection<STATE>;
template class wgma::materials::SolutionProjection<CSTATE>;
