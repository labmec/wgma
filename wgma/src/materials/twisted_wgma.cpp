#include "materials/twisted_wgma.hpp"
#include <TPZMaterialDataT.h>
#include <cmath>
using namespace std::complex_literals;

namespace wgma::materials{
  //! Gets the equivalent permeability of the material
  void TwistedWgma::GetPermeability([[maybe_unused]] const TPZVec<REAL> &x,
                                    TPZFMatrix<CSTATE> &ur) const
  {
    TPZFNMatrix<9,CSTATE> t, tt, tmp;
    TPZAnisoWgma::GetPermeability(x,ur);
    TransformationMatrix(t,x);
    t.Transpose(&tt);
    ur.Multiply(tt,tmp);
    t.Multiply(tmp,ur);
  }
  //! Gets the equivalent permittivity of the material
  void TwistedWgma::GetPermittivity([[maybe_unused]] const TPZVec<REAL> &x,
                                    TPZFMatrix<CSTATE> &er) const
  {

    
    TPZFNMatrix<9,CSTATE> t, tt, tmp;
    TPZAnisoWgma::GetPermittivity(x,er);
    TransformationMatrix(t,x);
    t.Transpose(&tt);
    er.Multiply(tt,tmp);
    t.Multiply(tmp,er);
  }

  void TwistedWgma::TransformationMatrix(TPZFMatrix<CSTATE> &mat,
                                         const TPZVec<REAL> &x) const
  {

    /**
       let us denote by x,y,z the cartesian coordinates,
       by u,v,w the helicoidal coordinates and by
       a the torsion parameter

       we have
       J(u,v,w) =
       |  cos(aw)    sin(aw)    av cos(aw) - au sin(aw) |
       | -sin(aw)    cos(aw)   -au cos(aw) - av sin(aw) |
       |    0           0                  1            |

       and

       J(u,v,w)^-1 = 

       |  cos(aw)   -sin(aw)    -av |
       |  sin(aw)    cos(aw)     au |
       |    0           0         1 |


       setting w = z = 0

       J(u,v)^-1 =

       |  1   0    -av |
       |  0   1     au |
       |  0   0      1 |
       
       and

       T(u,v) = J^T J/det(J) =
       |    1           0          av       |
       |    0           1         -au       |
       |   av          -au   1+a*a(u*u+v*v) |
    */

    const auto &u = x[0];
    const auto &v = x[1];
    const auto &a = this->GetAlpha();
    mat.Redim(3,3);
    mat.Identity();
    mat.Put(0,2, -a*v);
    mat.Put(1,2, a*u);
  }

  void TwistedWgma::Solution(
      const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
      int var,TPZVec<CSTATE> &solout)
  {
    TPZFNMatrix<3,CSTATE> sol(3,1), tsol;
    sol.Put(0,0,datavec[fHCurlMeshIndex].sol[0][0]);
    sol.Put(1,0,datavec[fHCurlMeshIndex].sol[0][1]);
    sol.Put(2,0,datavec[fH1MeshIndex].sol[0][0]);
    TPZFNMatrix<9,CSTATE> mat, tmat;
    TransformationMatrix(mat,datavec[0].x);
    mat.Transpose(&tmat);
    tmat.Multiply(sol, tsol);

    switch (var) {
    case 0:{//et_real
      for (int i = 0; i < 2; ++i) {
        solout[i] = std::real(tsol.Get(i,0));
      }
      break;
    }
    case 1:{//ez_real
      solout[0] = std::real(tsol.Get(2,0));
      break;
    }
    case 2:{//et_abs
      for (int i = 0; i < 2; ++i) {
        solout[i] = std::abs(tsol.Get(i,0));
      }
      break;
    }
    case 3:{//ez_abs
      solout[0] = std::abs(tsol.Get(2,0));
      break;
    }
    case 4:{//material
      TPZFNMatrix<9,CSTATE> er(3,3,0.);
      GetPermittivity(datavec[0].x, er);
      solout[0] = er.Get(0,0);
      solout[1] = er.Get(1,1);
      break;
    }
    case 5:{//pOrder
      solout.Resize(1);
      solout[0] = datavec[fH1MeshIndex].p;
      break;
    }
    case 6:{//pOrder
      solout.Resize(1);
      solout[0] = datavec[fHCurlMeshIndex].p;
      break;
    }
    default:
      DebugStop();
      break;
    }
  }

#ifdef CUSTOMPML


  void TwistedWgmaPML::ComputeCylindricalPMLMat(TPZFMatrix<CSTATE> &mat,
                                                const TPZVec<REAL> &x_hel) const
  {
    const auto &u=x_hel[0];
    const auto &v=x_hel[1];
    const STATE r = sqrt(u*u+v*v);
    const CSTATE imag = 1i;
    const CSTATE alpha_pml = this->fAlphaMaxR;
    const CSTATE alpha = this->GetAlpha();
    const auto &rmin = this->fPmlBeginR;
    const auto &dr = this->fDR;
    const CSTATE diffr = r-rmin;
    const CSTATE sr = 1.- imag*alpha_pml * diffr*diffr/(dr*dr);
    //integral of sr(r') frodm rmin to r
    const CSTATE rt = r - imag*alpha_pml*diffr*diffr*diffr/(3*dr*dr);
    const CSTATE art = alpha*rt;
    const STATE phi = std::atan2(v,u);

    /*
      first we compute
                                           T1
                   ( r sr/rt    0                          0              )
     Tpml = R(phi) ( 0        rt/(r sr)             -alpha rt/sr          ) R(-phi)
                   ( 0      -alpha rt/sr  r*(1+alpha*alpha*rt*rt)/(rt*sr) )
     */
    mat.Redim(3,3);
    mat.Put(0,0, sr*r/rt);
    mat.Put(1,1, rt/(sr*r));
    mat.Put(1,2,-alpha* rt/ sr);
    mat.Put(2,1,-alpha* rt/ sr);
    mat.Put(2,2,r*(1.+alpha*alpha*rt*rt)/(rt*sr));
  }
  
// #define CUSTOMPML2  
  void TwistedWgmaPML::GetPermeability([[maybe_unused]] const TPZVec<REAL> &x_hel,
                                       TPZFMatrix<CSTATE> &urmat) const
  {
    
    TPZFNMatrix<9,CSTATE> t1, t2;

    auto RotationMatrix = [](TPZFMatrix<CSTATE> &mat, STATE theta){
      mat.Redim(3,3);
      mat.Put(0,0, std::cos(theta));
      mat.Put(0,1,-std::sin(theta));
      mat.Put(1,0, std::sin(theta));
      mat.Put(1,1, std::cos(theta));
      mat.Put(2,2,1);
    };

    this->ComputeCylindricalPMLMat(t1,x_hel);

    const auto &u = x_hel[0];
    const auto &v = x_hel[1];
    const STATE phi = std::atan2(v,u);
    RotationMatrix(t2,phi);

    TPZAnisoWgma::GetPermeability(x_hel, urmat);
    const CSTATE ur = urmat.Get(0,0);
    urmat.Zero();
    
    t2.Multiply(t1,urmat);
    RotationMatrix(t2,-phi);
    urmat.Multiply(t2,t1);
    //now t1 is Tpml
    urmat.Identity();
    t1.SolveDirect(urmat,ELU);
    urmat *= ur;
  }
  //! Gets the equivalent permittivity of the material
  void TwistedWgmaPML::GetPermittivity([[maybe_unused]] const TPZVec<REAL> &x_hel,
                                       TPZFMatrix<CSTATE> &ermat) const
  {

    TPZFNMatrix<9,CSTATE> t1, t2;

    auto RotationMatrix = [](TPZFMatrix<CSTATE> &mat, STATE theta){
      mat.Redim(3,3);
      mat.Put(0,0, std::cos(theta));
      mat.Put(0,1,-std::sin(theta));
      mat.Put(1,0, std::sin(theta));
      mat.Put(1,1, std::cos(theta));
      mat.Put(2,2,1);
    };

    this->ComputeCylindricalPMLMat(t1,x_hel);

    const auto &u = x_hel[0];
    const auto &v = x_hel[1];
    const STATE phi = std::atan2(v,u);
    RotationMatrix(t2,phi);

    TPZAnisoWgma::GetPermittivity(x_hel, ermat);
    const CSTATE er = ermat.Get(0,0);
    ermat.Zero();    
    t2.Multiply(t1,ermat);
    RotationMatrix(t2,-phi);
    ermat.Multiply(t2,t1);
    //now t1 is Tpml
    ermat.Identity();
    t1.SolveDirect(ermat,ELU);
    ermat *= er;
    
  }

#endif
  
};


#ifndef CUSTOMPML
template class TPZCombinedSpacesCylindricalPML<wgma::materials::TwistedWgma>;
#endif
