#include "materials/acousticmodesbase.hpp"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"

using namespace std::complex_literals;
using namespace wgma::materials;

AcousticModesBase::AcousticModesBase(int id, 
                             const STATE mu,
                             const STATE lambda,
                             const STATE rho,
                             const REAL scale) :
    fMu(mu), fRho(rho), fLambda(lambda),
    fScaleFactor(scale)
{
}

void
AcousticModesBase::ComputeLxy(const TPZFMatrix<STATE> &dphi,
                          const TPZVec<REAL> &xvec,
                          TPZFMatrix<CSTATE> &lxy)
{
    const int nphi = dphi.Cols();
    lxy.Redim(6,3*nphi);
    for(int i = 0; i < nphi; i++){
        const auto x = xvec[0];
        const auto y = xvec[1];
        const auto dphix = dphi.Get(0,i);
        const auto dphiy = dphi.Get(1,i);
        //(-x dudy + y dudx) * tau
        const auto lambda = (-x*dphiy + y*dphix)*0.5;
        //dudx
        lxy.Put(0,3*i+0, dphix);
        //dvdy
        lxy.Put(1,3*i+1, dphiy);
        //dudy
        lxy.Put(3,3*i+0, dphiy);
        //dvdx
        lxy.Put(3,3*i+1, dphix);
        
        lxy.Put(4,3*i+0, lambda);
        //dwdx
        lxy.Put(4,3*i+2, dphix);

        lxy.Put(5,3*i+1, lambda);
        //dwdy
        lxy.Put(5,3*i+2, dphiy);
    }
    
}
void
AcousticModesBase::ComputeLz(const TPZFMatrix<STATE> &phi,
                         TPZFMatrix<CSTATE> &lz)
{
    const int nphi = phi.Rows();
    lz.Redim(6,3*nphi);
    for(int i = 0; i < nphi; i++){
        const auto phival = phi.GetVal(i,0);
        lz.Put(2,3*i+2,phival);
        lz.Put(4,3*i+0,phival);
        lz.Put(5,3*i+1,phival);
    }
}

void
AcousticModesBase::ComputeC(TPZFMatrix<CSTATE> &Cmat)
{
    Cmat.Redim(6,6);
    const REAL scale = this->fScaleFactor;
    const REAL lambda = GetLambda()*scale*scale;
    const REAL mu = GetMu()*scale*scale;
    
    for(int i = 0; i < 3; i++){
        Cmat.Put(i,i,2*mu);
        for(int j = 0; j < 3; j++){
            Cmat(i,j) += lambda;
        }
    }
    for(int i = 3; i < 6; i++){
        Cmat.Put(i,i,1*mu);
    }
}

int AcousticModesBase::VariableIndex(const std::string &name) const
{
    if( strcmp(name.c_str(), "u_real") == 0) return 0;
    if( strcmp(name.c_str(), "u_abs") == 0) return 1;
    DebugStop();
    return 1;
}

int AcousticModesBase::NSolutionVariables(int var) const
{
    switch (var) {
        case 0: //u_real
            return 3;
        case 1://u_abs
            return 3;
        default:
            DebugStop();
            break;
    }
    return 1;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void AcousticModesBase::Solution(
    const TPZMaterialDataT<CSTATE> &data,
    int var,
    TPZVec<CSTATE> &solout)
{
    const auto & uvec = data.sol[0];
    switch (var) {
    case 0:{//u_real
        for (int i = 0; i < uvec.size(); ++i) {
            solout[i] = std::real(uvec[i]);
        }
        break;
    }
    case 1:{//u_abs
        for (int i = 0; i < uvec.size(); ++i) {
            solout[i] = std::abs(uvec[i]);
        }
        break;
    }
    default:
        DebugStop();
        break;
    }
}