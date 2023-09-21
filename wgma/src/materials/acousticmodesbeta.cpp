#include "materials/acousticmodesbeta.hpp"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"

using namespace std::complex_literals;
using namespace wgma::materials;

AcousticModesBeta::AcousticModesBeta() : TBase()
{
    SetMatrixK();
}

AcousticModesBeta::AcousticModesBeta(int id) : TBase(id){
    SetMatrixK();
}

AcousticModesBeta::AcousticModesBeta(int id, 
                             const STATE mu,
                             const STATE lambda,
                             const STATE rho,
                             const STATE freq,
                             const REAL scale) :
    TBase(id), AcousticModesBase(id,mu,lambda,rho,scale),
    fFreq(freq)
{
    SetMatrixK();

    if (freq <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative wavelength. Aborting..\n";
        DebugStop();
    }
}


AcousticModesBeta* AcousticModesBeta::NewMaterial() const{
    return new AcousticModesBeta();
}

void AcousticModesBeta::SetMatrixK()
{
    TPZMatQuadraticEigenVal::SetMatrixK();
    fCurrentContribute = [this](const auto &arg1, auto arg2, auto &arg3,
                                auto &arg4){
        this->ContributeK(arg1,arg2,arg3,arg4);
    };
    fCurrentContributeBC = [this](const auto &arg1, auto arg2, auto &arg3,
                                  auto &arg4, auto &arg5){
        this->ContributeBCK(arg1,arg2,arg3,arg4,arg5);
    };
}

void AcousticModesBeta::SetMatrixL()
{
    TPZMatQuadraticEigenVal::SetMatrixL();
    fCurrentContribute = [this](const auto &arg1, auto arg2, auto &arg3,
                                auto &arg4){
        this->ContributeL(arg1,arg2,arg3,arg4);
    };
    fCurrentContributeBC = [this](const auto &arg1, auto arg2, auto &arg3,
                                  auto &arg4, auto &arg5){
        this->ContributeBCL(arg1,arg2,arg3,arg4,arg5);
    };
}


void AcousticModesBeta::SetMatrixM()
{
    TPZMatQuadraticEigenVal::SetMatrixM();
    fCurrentContribute = [this](const auto &arg1, auto arg2, auto &arg3,
                                auto &arg4){
        this->ContributeM(arg1,arg2,arg3,arg4);
    };
    fCurrentContributeBC = [this](const auto &arg1, auto arg2, auto &arg3,
                                  auto &arg4, auto &arg5){
        this->ContributeBCM(arg1,arg2,arg3,arg4,arg5);
    };
}

void
AcousticModesBeta::Contribute(
    const TPZMaterialDataT<CSTATE> &data,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)
{
    fCurrentContribute(data,weight,ek,ef);
}

void AcousticModesBeta::ContributeBC(
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
AcousticModesBeta::ContributeK(
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

    TPZFNMatrix<3*nphimax,CSTATE> phivec(3,3*nphi,0);

    //L matrix
    TPZFNMatrix<6*3*nphimax,CSTATE> lxy(6,3*nphi,0);
    ComputeLxy(grad_phi,data.x,lxy);

    for(int i = 0; i < nphi; i++){
        const auto phival = phi.Get(i,0);
        phivec.Put(0,3*i+0,phival);
        phivec.Put(1,3*i+1,phival);
        phivec.Put(2,3*i+2,phival);
    }

    //C matrix in voigt notation
    TPZFNMatrix<6*6,CSTATE> Cmat(6,6,0);

    ComputeC(Cmat);
    

    const REAL scale = this->fScaleFactor;
    const REAL freq = GetFrequency();
    const REAL rho = GetRho();
    const REAL omega = 2*M_PI*freq;
    const REAL omega2rho = (omega*scale)*(omega*scale)*(rho*scale);
    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/
    TPZFNMatrix<3*nphimax*3*nphimax,CSTATE> tmp;

    //lxy^T C lxy
    Cmat.Multiply(lxy,tmp);
    ek.AddContribution(0,0,lxy,true,tmp,false,weight);
    //-w^2 rho phi_i phi_j
    ek.AddContribution(0,0,phivec,true,phivec,false,-omega2rho*weight);
}

void
AcousticModesBeta::ContributeL(
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

    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/
    TPZFNMatrix<3*nphimax*3*nphimax,CSTATE> tmp;

    //lz^T C lxy
    Cmat.Multiply(lxy,tmp);
    ek.AddContribution(0,0,lz,true,tmp,false,   weight);
    //lxy^T C lz
    Cmat.Multiply(lz,tmp);
    ek.AddContribution(0,0,lxy,true,tmp,false, -weight);
}

void
AcousticModesBeta::ContributeM(
    const TPZMaterialDataT<CSTATE> &data,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)
{
    
    const auto &phi = data.phi;
    const int nphi = phi.Rows();

    //L matrix
    TPZFNMatrix<6*3*nphimax,CSTATE> lz(6,3*nphi,0);

    ComputeLz(phi,lz);
    
    //C matrix in voigt notation
    TPZFNMatrix<6*6,CSTATE> Cmat(6,6,0);
    ComputeC(Cmat);

    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/
    TPZFNMatrix<3*nphimax*3*nphimax,CSTATE> tmp;

    //lz^T C lz
    Cmat.Multiply(lz,tmp);
    ek.AddContribution(0,0,lz,true,tmp,false,-weight);
}


void AcousticModesBeta::ContributeBCK(
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