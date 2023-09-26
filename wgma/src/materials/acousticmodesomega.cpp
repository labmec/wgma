#include "materials/acousticmodesomega.hpp"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"

using namespace std::complex_literals;
using namespace wgma::materials;

AcousticModesOmega::AcousticModesOmega() : TBase()
{
    SetMatrixA();
}

AcousticModesOmega::AcousticModesOmega(int id) : TBase(id){
    SetMatrixA();
}

AcousticModesOmega::AcousticModesOmega(int id, 
                             const STATE mu,
                             const STATE lambda,
                             const STATE rho,
                             const CSTATE beta,
                             const REAL scale) :
    TBase(id), AcousticModesBase(id,mu,lambda,rho,scale),
    fBeta(beta)
{
    SetMatrixA();
}


AcousticModesOmega* AcousticModesOmega::NewMaterial() const{
    return new AcousticModesOmega();
}

void AcousticModesOmega::SetMatrixA()
{
    TPZMatGeneralisedEigenVal::SetMatrixA();
    fCurrentContribute = [this](const auto &arg1, auto arg2, auto &arg3,
                                auto &arg4){
        this->ContributeA(arg1,arg2,arg3,arg4);
    };
    fCurrentContributeBC = [this](const auto &arg1, auto arg2, auto &arg3,
                                  auto &arg4, auto &arg5){
        this->ContributeBCA(arg1,arg2,arg3,arg4,arg5);
    };
}

void AcousticModesOmega::SetMatrixB()
{
    TPZMatGeneralisedEigenVal::SetMatrixB();
    fCurrentContribute = [this](const auto &arg1, auto arg2, auto &arg3,
                                auto &arg4){
        this->ContributeB(arg1,arg2,arg3,arg4);
    };
    fCurrentContributeBC = [this](const auto &arg1, auto arg2, auto &arg3,
                                  auto &arg4, auto &arg5){
        this->ContributeBCB(arg1,arg2,arg3,arg4,arg5);
    };
}

void
AcousticModesOmega::Contribute(
    const TPZMaterialDataT<CSTATE> &data,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)
{
    fCurrentContribute(data,weight,ek,ef);
}

void AcousticModesOmega::ContributeBC(
    const TPZMaterialDataT<CSTATE> &data,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    TPZBndCondT<CSTATE> &bc)
{
    fCurrentContributeBC(data,weight,ek,ef,bc);
}



//we will preallocate memory considering at most 30 shape functions
constexpr int nphimax{30};

void
AcousticModesOmega::ContributeA(
    const TPZMaterialDataT<CSTATE> &data,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)
{

    const auto &phi = data.phi;
    const int nphi = phi.Rows();
    TPZFNMatrix<3*nphimax,REAL> grad_phi(3, nphi, 0.);
    {
        const TPZFMatrix<REAL> &dphidaxes = data.dphix;
        TPZAxesTools<REAL>::Axes2XYZ(dphidaxes, grad_phi,
                                     data.axes);
    }

    //L matrix
    TPZFNMatrix<6*3*nphimax,CSTATE> lxy(6,3*nphi,0), lz(6,3*nphi,0);
    
    ComputeLxy(grad_phi,data.x,lxy);
    ComputeLz(phi,lz);

    //C matrix in voigt notation
    TPZFNMatrix<6*6,CSTATE> Cmat(6,6,0);

    ComputeC(Cmat);
    

    const REAL scale = this->fScaleFactor;
    const CSTATE beta = GetBeta()*scale;

    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/
    TPZFNMatrix<3*nphimax*3*nphimax,CSTATE> tmp;
    Cmat.Multiply(lxy,tmp);
    //lxy^T C lxy
    ek.AddContribution(0,0,lxy,true,tmp,false,weight);
    //lz^T C lxy
    ek.AddContribution(0,0,lz,true,tmp,false,  1i*beta*weight);
    Cmat.Multiply(lz,tmp);
    //lxy^T C lz
    ek.AddContribution(0,0,lxy,true,tmp,false, -1i*beta*weight);
    //lz^T C lz
    ek.AddContribution(0,0,lz,true,tmp,false, beta*beta*weight);
    
    
}

void
AcousticModesOmega::ContributeB(
    const TPZMaterialDataT<CSTATE> &data,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)
{
 
    const auto &phi = data.phi;
    const int nphi = phi.Rows();


    TPZFNMatrix<3*nphimax,CSTATE> phivec(3,3*nphi,0);
    for(int i = 0; i < nphi; i++){
        const auto phival = phi.Get(i,0);
        phivec.Put(0,3*i+0,phival);
        phivec.Put(1,3*i+1,phival);
        phivec.Put(2,3*i+2,phival);
    }
    const REAL scale = this->fScaleFactor;
    const REAL rho = this->GetRho()*scale*scale*scale;
    ek.AddContribution(0,0,phivec,true,phivec,false,rho*weight);
}



void AcousticModesOmega::ContributeBCA(
    const TPZMaterialDataT<CSTATE> &data,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    TPZBndCondT<CSTATE> &bc)
{
    const TPZFMatrix<REAL> &phi = data.phi;
    const int nphi  = phi.Rows();
        
    
    const auto& BIG = TPZMaterial::fBigNumber;
    
    const CSTATE v1 = bc.Val1()(0,0);
    const CSTATE v2 = bc.Val2()[0];
    constexpr STATE tol = std::numeric_limits<STATE>::epsilon();
    if(std::abs(v2) > tol){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThis method supports only homogeneous boundary conditions.\n";
        std::cout<<"Stopping now..."<<std::endl;
        DebugStop();
    }
    switch ( bc.Type() )
    {
    case 0:{
        const int nfuncs = nphi*3;
        for(int i = 0 ; i<nfuncs ; i++){
            const int iphi = i/3;
            for(int j=0;j<nfuncs;j++){
                const int jphi = i/3;
                const STATE stiff = phi(iphi,0) * phi(jphi,0) * BIG ;
                ek(i,j) += stiff*weight;
            }
        }
        break;
    }
    case 1:
        ///PMC condition just adds zero to both matrices. nothing to do here....
        break;
    case 2:
        /// periodic conditions are treated at a mesh level
        break;
    default:
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThis module supports only dirichlet and neumann boundary conditions.\n";
        PZError<<"Stopping now..."<<std::endl;
        DebugStop();
        break;
    }
}