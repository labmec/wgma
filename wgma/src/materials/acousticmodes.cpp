#include "materials/acousticmodes.hpp"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"

using namespace std::complex_literals;
using namespace wgma::materials;

AcousticModes::AcousticModes() : TBase()
{
    SetMatrixK();
}

AcousticModes::AcousticModes(int id) : TBase(id){
    SetMatrixK();
}

AcousticModes::AcousticModes(int id, 
                             const STATE mu,
                             const STATE lambda,
                             const STATE rho,
                             const STATE freq,
                             const REAL scale) :
    TBase(id), fMu(mu), fRho(rho), fLambda(lambda),
    fFreq(freq),fScaleFactor(scale)
{
    SetMatrixK();

    if (freq <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative wavelength. Aborting..\n";
        DebugStop();
    }
}


AcousticModes* AcousticModes::NewMaterial() const{
    return new AcousticModes();
}

void AcousticModes::SetMatrixK()
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

void AcousticModes::SetMatrixL()
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


void AcousticModes::SetMatrixM()
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
AcousticModes::Contribute(
    const TPZMaterialDataT<CSTATE> &data,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)
{
    fCurrentContribute(data,weight,ek,ef);
}

void AcousticModes::ContributeBC(
    const TPZMaterialDataT<CSTATE> &data,
    REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
    TPZBndCondT<CSTATE> &bc)
{
    fCurrentContributeBC(data,weight,ek,ef,bc);
}

void
AcousticModes::ComputeLxy(const TPZFMatrix<STATE> &dphi,
                          const TPZVec<REAL> &x,
                          TPZFMatrix<CSTATE> &lxy)
{
    const int nphi = dphi.Cols();
    lxy.Redim(6,3*nphi);
    for(int i = 0; i < nphi; i++){
        const auto dphix = dphi.Get(0,i);
        const auto dphiy = dphi.Get(1,i);
        //dudx
        lxy.Put(0,3*i+0, dphix);
        //dvdy
        lxy.Put(1,3*i+1, dphiy);
        //dudy
        lxy.Put(3,3*i+0, dphiy);
        //dvdx
        lxy.Put(3,3*i+1, dphix);
        //dwdx
        lxy.Put(4,3*i+2, dphix);
        //dwdy
        lxy.Put(5,3*i+2, dphiy);
    }
    
}
void
AcousticModes::ComputeLz(const TPZFMatrix<STATE> &phi,
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
AcousticModes::ComputeC(TPZFMatrix<CSTATE> &Cmat)
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

//we will preallocate memory considering at most 30 shape functions
constexpr int nphimax{30};

void
AcousticModes::ContributeK(
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

    TPZFNMatrix<3*nphimax,CSTATE> phivec(3*nphi,3,0);

    //L matrix
    TPZFNMatrix<6*3*nphimax,CSTATE> lxy(6,3*nphi,0);
    ComputeLxy(grad_phi,data.x,lxy);

    for(int i = 0; i < nphi; i++){
        phivec.Put(3*i+0,0,phi.Get(i,0));
        phivec.Put(3*i+1,1,phi.Get(i,0));
        phivec.Put(3*i+2,2,phi.Get(i,0));
    }

    //C matrix in voigt notation
    TPZFNMatrix<6*6,CSTATE> Cmat(6,6,0);

    ComputeC(Cmat);
    

    const REAL scale = this->fScaleFactor;
    const REAL freq = GetFrequency();
    const REAL rho = GetRho();
    const REAL omega = 2*M_PI*freq;
    const REAL omega2rho = omega*omega*rho*scale*scale*scale;
    /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/
    TPZFNMatrix<3*nphimax*3*nphimax,CSTATE> tmp;

    //lxy^T C lxy
    TPZFNMatrix<6*3*nphimax,CSTATE> lxyt(3*nphi,6,0);
    lxy.Transpose(&lxyt);
    lxyt.Multiply(Cmat,tmp,0);
    ek.AddContribution(0,0,tmp,false,lxy,false,weight);
    //-w^2 rho phi_i phi_j
    ek.AddContribution(0,0,phivec,false,phivec,true,-omega2rho*weight);
}

void
AcousticModes::ContributeL(
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

    //lxy^T C lz
    TPZFNMatrix<6*3*nphimax,CSTATE> lxyt(3*nphi,6,0);
    lxy.Transpose(&lxyt);
    lxyt.Multiply(Cmat,tmp,0);
    ek.AddContribution(0,0,tmp,false,lz,false,  weight);
    //lz^T C lxy
    TPZFNMatrix<6*3*nphimax,CSTATE> lzt(3*nphi,6,0);
    lz.Transpose(&lzt);
    lzt.Multiply(Cmat,tmp,0);
    ek.AddContribution(0,0,tmp,false,lxy,false, -weight);
}

void
AcousticModes::ContributeM(
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
    TPZFNMatrix<6*3*nphimax,CSTATE> lzt(3*nphi,6,0);
    lz.Transpose(&lzt);
    lzt.Multiply(Cmat,tmp,0);
    ek.AddContribution(0,0,tmp,false,lz,false,-weight);
}


void AcousticModes::ContributeBCK(
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

int AcousticModes::VariableIndex(const std::string &name) const
{
    if( strcmp(name.c_str(), "u_real") == 0) return 0;
    if( strcmp(name.c_str(), "u_abs") == 0) return 1;
    DebugStop();
    return 1;
}

int AcousticModes::NSolutionVariables(int var) const
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
void AcousticModes::Solution(
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