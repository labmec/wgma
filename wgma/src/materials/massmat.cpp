#include "materials/massmat.hpp"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"

using namespace wgma::materials;

template<class TVar>
MassMatrix<TVar>::MassMatrix(int id, int dim, int nstate) :
    TPZRegisterClassId(&MassMatrix::ClassId),
    TBase(id), fDim(dim), fNStateVars(nstate), fSol(nstate,0.)
{
}
template<class TVar>
MassMatrix<TVar>::MassMatrix(int id, int dim, int nstate,
                                       const TPZVec<TVar> &sol) :
    TPZRegisterClassId(&MassMatrix::ClassId),
    TBase(id), fDim(dim), fNStateVars(nstate), fSol(sol)
{
    if(fSol.size()!=fNStateVars){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nn state variables: "<<nstate;
        PZError<<"\nn solutions: "<<fSol.size();
        PZError<<"\nAborting...\n";
        DebugStop();
    }
}

template<class TVar>
TPZMaterial * MassMatrix<TVar>::NewMaterial() const{
	return new MassMatrix(*this);
}

template<class TVar>
void MassMatrix<TVar>::Contribute(const TPZMaterialDataT<TVar> &data,
                                       REAL weight,
                                       TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef){
	
	const int nshape = data.phi.Rows();
    const int nvars = fNStateVars;
    TPZManVector<TVar,10> solLoc(fSol);
    if(this->HasForcingFunction()){
        this->fForcingFunction(data.x,solLoc);
    }
    const auto &phi = data.phi;
	for(int i = 0; i < nshape; i++){
		for(int j = 0; j < nshape; j++){
            const STATE phiIphiJ = phi.GetVal(i,0) * phi.GetVal(j,0);
			for(int ivi = 0; ivi < nvars; ivi++){
                const int posI = nvars*i+ivi;
                const int posJ = nvars*j+ivi;
                ek(posI, posJ) += weight*fScale*phiIphiJ;
			}//ivi
		}//for j
		for(int ivi = 0; ivi < nvars; ivi++){
			const int posI = nvars*i+ivi;
            ef(posI,0) += weight*fScale*phi.GetVal(i,0);
		}//ivi
	}//for i
}

template<class TVar>
void MassMatrix<TVar>::ContributeBC(const TPZMaterialDataT<TVar> &data,
                                         REAL weight,
                                         TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                                         TPZBndCondT<TVar> &bc)
{
    return;
}
template<class TVar>
void MassMatrix<TVar>::GetSolDimensions(uint64_t &u_len,
                                             uint64_t &du_row,
                                             uint64_t &du_col) const
{
    u_len=1;
    du_row=3;
    du_col=1;
}


template<class TVar>
int MassMatrix<TVar>::VariableIndex(const std::string &name) const{
	if(!strcmp("Solution",name.c_str())) return ESolution;
    if(!strcmp("Derivative",name.c_str())) return EDerivative;
	return TPZMaterial::VariableIndex(name);
}

template<class TVar>
int MassMatrix<TVar>::NSolutionVariables(int var) const{
	if(var == ESolution) return 1;
    if (var == EDerivative) {
        return fDim;
    }
	
    return TPZMaterial::NSolutionVariables(var);
}

template<class TVar>
void MassMatrix<TVar>::Solution(const TPZMaterialDataT<TVar> &data,
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
int MassMatrix<TVar>::ClassId() const{
    return Hash("MassMatrix") ^ TBase::ClassId() << 1;
}


template class wgma::materials::MassMatrix<STATE>;
template class wgma::materials::MassMatrix<CSTATE>;
