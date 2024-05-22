#include "materials/solutionprojection.hpp"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"

using namespace wgma::materials;

template<class TVar>
SolutionProjection<TVar>::SolutionProjection(int id, int dim, int soldim) :
  TPZRegisterClassId(&SolutionProjection::ClassId),
  TBase(id), fDim(dim), fSolDim(soldim)
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
  if(sz!=nc || sz != fSolDim){
    std::ostringstream sout;
    sout<<__PRETTY_FUNCTION__
	<<"Incompatible dimensions!\n"
	<<"sz:"<<sz<<" nc: "<<nc<<" sol dim: "<<fSolDim;
    PZError<<sout.str()<<std::endl;
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
  u_len=fSolDim;
  if(fSolDim==1){//scalar field
    du_row=3;
  }else{//hcurl field
    du_row = fDim == 1 ? 1 : 2*fDim - 3;
  }
  du_col=1;
}


template<class TVar>
int SolutionProjection<TVar>::VariableIndex(const std::string &name) const{
	if(!strcmp("Solution",name.c_str()) || !strcmp("Solution_real",name.c_str())){
    return ESolReal;
  }
  if(!strcmp("Solution_imag",name.c_str())) {return ESolImag;}
  if(!strcmp("Solution_abs",name.c_str())) {return ESolAbs;}
  if(!strcmp("Derivative",name.c_str())) return EDerivative;
	return TPZMaterial::VariableIndex(name);
}

template<class TVar>
int SolutionProjection<TVar>::NSolutionVariables(int var) const{
	if(var == ESolReal || var == ESolAbs || var == ESolImag) {return fSolDim;}
  if (var == EDerivative) {
    if(fSolDim==1) return fDim;
    return fDim == 1 ? 1: 2*fDim - 3;
  }
	
  return TPZMaterial::NSolutionVariables(var);
}

template<class TVar>
void SolutionProjection<TVar>::Solution(const TPZMaterialDataT<TVar> &data,
                                        int var, TPZVec<TVar> &solOut)
{
  const auto &sol = data.sol[0];
	if(var == ESolReal || var == ESolAbs || var == ESolImag){
    solOut.Resize(sol.size());
    if(var == ESolReal){
      for (int i=0; i<sol.size(); i++) {
        solOut[i] = std::real(sol[i])/fScale;
      }
    }else if (var == ESolAbs){
      for (int i=0; i<sol.size(); i++) {
        solOut[i] = std::abs(sol[i])/fScale;
      }
    }else{
      for (int i=0; i<sol.size(); i++) {
        solOut[i] = std::imag(sol[i])/fScale;
      }
    }
		return;
	}
  if (var == EDerivative && fSolDim==1) {
    const auto &dsol = data.dsol[0];
    solOut.Resize(fDim);
    for (int i=0; i<fDim; i++) {
      solOut[i] = std::real(dsol.GetVal(i,0))/fScale;
    }
    return;
  }
  else if (var == EDerivative && fSolDim==3) {
    const int curldim = fDim == 1 ? 1: 2*fDim - 3;
    const auto &curlsol = data.curlsol[0];
    solOut.Resize(curldim);
    for (int i=0; i<curldim; i++) {
      solOut[i] = std::real(curlsol[i])/fScale;
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
