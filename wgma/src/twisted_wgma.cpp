#include "twisted_wgma.hpp"
#include <TPZMaterialDataT.h>
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
    mat.Put(1,2,  a*u);
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


  void TwistedWgmaPML::GetPermeability([[maybe_unused]] const TPZVec<REAL> &x,
                                       TPZFMatrix<CSTATE> &ur) const
  {
    TPZAnisoWgma::GetPermeability(x,ur);
    CSTATE sr{1}, sz{1};
    const auto r = sqrt(x[0]*x[0]+x[1]*x[1]);
    const auto z = x[2];
    ComputeSParameters(r,z,sr,sz);
    const auto imagsr = sr.imag();
    const auto sx = CSTATE(1. + 1i*imagsr*x[0]/r);
    const auto sy = CSTATE(1. + 1i*imagsr*x[1]/r);
    const auto dets = sx*sy*sz;
    TPZFNMatrix<9,CSTATE> smat(3,3,0.), tmp(3,3,0.);
    smat.PutVal(0,0,sx);
    smat.PutVal(1,1,sy);
    smat.PutVal(2,2,sz);
    smat.Multiply(ur,tmp);
    tmp.Multiply(smat,ur);
    ur *= dets;

    TPZFNMatrix<9,CSTATE> t, tt;
    TransformationMatrix(t,x);
    t.Transpose(&tt);
    ur.Multiply(tt,tmp);
    t.Multiply(tmp,ur);
  }
  //! Gets the equivalent permittivity of the material
  void TwistedWgmaPML::GetPermittivity([[maybe_unused]] const TPZVec<REAL> &x,
                                       TPZFMatrix<CSTATE> &er) const
  {

    
    TPZAnisoWgma::GetPermittivity(x,er);
    CSTATE sr{1}, sz{1};
    const auto r = sqrt(x[0]*x[0]+x[1]*x[1]);
    const auto z = x[2];
    ComputeSParameters(r,z,sr,sz);
    const auto imagsr = sr.imag();
    const auto sx = CSTATE(1. + 1i*imagsr*x[0]/r);
    const auto sy = CSTATE(1. + 1i*imagsr*x[1]/r);
    const auto dets = sx*sy*sz;
    TPZFNMatrix<9,CSTATE> smat(3,3,0.), tmp(3,3,0.);
    smat.PutVal(0,0,sx);
    smat.PutVal(1,1,sy);
    smat.PutVal(2,2,sz);
    smat.Multiply(er,tmp);
    tmp.Multiply(smat,er);
    er *= dets;

    TPZFNMatrix<9,CSTATE> t, tt;
    TransformationMatrix(t,x);
    t.Transpose(&tt);
    er.Multiply(tt,tmp);
    t.Multiply(tmp,er);
    
  }
  
};