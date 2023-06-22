#include "twisted_wgma.hpp"

namespace wgma::materials{
  //! Gets the equivalent permeability of the material
  void TwistedWgma::GetPermeability([[maybe_unused]] const TPZVec<REAL> &x,
                                    TPZFMatrix<CSTATE> &ur) const
  {
    TPZFNMatrix<9,CSTATE> t, tt, urtmp, tmp;
    TPZAnisoWgma::GetPermeability(x,urtmp);
    TransformationMatrix(t,x);
    t.Transpose(&tt);
    urtmp.Multiply(tt,tmp);
    t.Multiply(tmp,ur);
  }
  //! Gets the equivalent permittivity of the material
  void TwistedWgma::GetPermittivity([[maybe_unused]] const TPZVec<REAL> &x,
                                    TPZFMatrix<CSTATE> &er) const
  {

    
    TPZFNMatrix<9,CSTATE> t, tt, ertmp, tmp;
    TPZAnisoWgma::GetPermittivity(x,ertmp);
    TransformationMatrix(t,x);
    t.Transpose(&tt);
    ertmp.Multiply(tt,tmp);
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
};