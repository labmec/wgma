#include "materials/scatteredfield.hpp"
#include "TPZMaterialDataT.h"

using namespace std::complex_literals;


using namespace wgma::materials;

ScatteredField::ScatteredField(int id, const CSTATE er,
                               const CSTATE ur, const STATE wl,
                               const REAL scale)
  : TBase(id),m_wl(wl), m_scale(scale)
{
  SetPermeability(ur);
  SetPermittivity(er);
  SetBackgroundPermittivity(er);
}
ScatteredField::ScatteredField(int id, const TPZFMatrix<CSTATE>& er,
                               const TPZFMatrix<CSTATE> &ur,
                               const STATE wl,
                               const REAL scale)
  : TBase(id),m_wl(wl), m_scale(scale)
{
  SetPermeability(ur);
  SetPermittivity(er);
  SetBackgroundPermittivity(er);
}


void ScatteredField::GetPermittivity(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &er) const
{
  er = m_er;
}

void ScatteredField::GetPermeability(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &ur) const
{
  ur = m_ur;
}

void ScatteredField::SetPermeability(CSTATE ur)
{
  m_ur.Redim(3,3);
  m_ur.PutVal(0,0,ur);
  m_ur.PutVal(1,1,ur);
  m_ur.PutVal(2,2,ur);
}


void ScatteredField::SetPermeability(const TPZFMatrix<CSTATE>& ur)
{
  if(ur.Rows()!=3 || ur.Cols()!= 3){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nSize of ur != 3. Aborting...\n";
    DebugStop();
  }
  m_ur = ur;
}
void ScatteredField::SetPermittivity(CSTATE er)
{
  m_er.Redim(3,3);
  m_er.PutVal(0,0,er);
  m_er.PutVal(1,1,er);
  m_er.PutVal(2,2,er);
}

void ScatteredField::SetPermittivity(const TPZFMatrix<CSTATE>&er)
{
  if(er.Rows()!=3 || er.Cols()!= 3){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nSize of er != 3. Aborting...\n";
    DebugStop();
  }
  m_er = er;
}

void ScatteredField::Contribute(const TPZMaterialDataT<CSTATE> &data,
                                REAL weight,
                                TPZFMatrix<CSTATE> &ek,
                                TPZFMatrix<CSTATE> &ef)
{
  //we now for sure that we wont have more than 200 shape functions
  constexpr int MEMSHAPE{3*200};
  TPZFNMatrix<9,CSTATE> er_mat,ur_inv_mat;
  GetPermittivity(data.x,er_mat);
  GetPermeability(data.x,ur_inv_mat);
  ur_inv_mat.Decompose(ELU);
  const int nshape = data.phi.Rows();
  const auto &phi_real = data.phi;
  const auto &curl_phi_real = data.curlphi;
  
  const STATE k0 = m_scale * 2*M_PI/m_wl;
  
  //making complex version of phi
  TPZFNMatrix<MEMSHAPE,CSTATE> phi(3,nshape);
  TPZFNMatrix<MEMSHAPE,CSTATE> curl_phi(3,nshape);

  //we cannot use phi_real_ptr because its actually transposed
  CSTATE *phi_ptr = phi.Elem();
  CSTATE *curl_phi_ptr = curl_phi.Elem();
  const STATE *curl_phi_real_ptr = curl_phi_real.Elem();
  const int sz = 3*nshape;
  for(int i = 0; i < sz; i++){
    //g instead of GetVal will be most likely inlined
    *phi_ptr++ = phi_real.g(i/3, i%3);
    *curl_phi_ptr++ = *curl_phi_real_ptr++;
  }

  constexpr int no_transp{0}, transp{1}, conj{2};
  TPZFNMatrix<MEMSHAPE,CSTATE> tmp;
  tmp = curl_phi;
  ur_inv_mat.Substitution(&tmp);
  ek.AddContribution(0, 0, curl_phi, conj, tmp, no_transp, weight);
  er_mat.Multiply(phi, tmp);
  ek.AddContribution(0, 0, phi, conj, tmp, no_transp, -k0*k0*weight);

  if(m_has_sol == false){return;}
  //now for the rhs we use the background field
  CSTATE const *solvec = &(data.sol[0][0]);
  CSTATE const *curlvec = &(data.curlsol[0][0]);
  TPZFNMatrix<3,CSTATE> sol(3,1,0),curl(3,1,0);

  CSTATE *sol_ptr = sol.Elem();
  CSTATE *curl_ptr = curl.Elem();
  for(int i = 0; i < 3; i++){
    *sol_ptr++ = *solvec++;
    *curl_ptr++ = *curlvec++;
  }

  TPZFNMatrix<9,CSTATE> er_background_mat;
  GetBackgroundPermittivity(er_background_mat);

  constexpr int sign{1};
  er_mat -= er_background_mat;
  er_mat.Multiply(sol, tmp);
  ef.AddContribution(0, 0, phi, conj, tmp, no_transp, sign*k0*k0*weight);

  // constexpr int sign{1};
  // tmp = curl;
  // ur_inv_mat.Substitution(&tmp);
  // ef.AddContribution(0, 0, curl_phi, conj, tmp, no_transp, -sign*weight);
  // er_mat.Multiply(sol, tmp);
  // ef.AddContribution(0, 0, phi, conj, tmp, no_transp, sign*k0*k0*weight);
  
  // tmp = sol;
  // ur_inv_mat.Substitution(&tmp);
  // ef.AddContribution(0, 0, phi, true, tmp, false, -sign*k0*k0*weight);
  // er_mat.Multiply(sol, tmp);
  // ef.AddContribution(0, 0, phi, true, tmp, false, sign*k0*k0*weight);
}


int ScatteredField::ClassId() const {
  return Hash("ScatteredField") ^
    TBase::ClassId() << 1;


}

//! Variable index of a given solution
int ScatteredField::VariableIndex(const std::string &name) const
{
  if( strcmp(name.c_str(), "Field_real") == 0) return 0;
  if( strcmp(name.c_str(), "Field_imag") == 0) return 1;
  if( strcmp(name.c_str(), "Field_abs") == 0) return 2;
  if( strcmp(name.c_str(), "Deriv_real") == 0) return 3;
  if( strcmp(name.c_str(), "Deriv_imag") == 0) return 4;
  if( strcmp(name.c_str(), "Deriv_abs") == 0) return 5;
  if( strcmp(name.c_str(), "Material") == 0) return 6;
  if( strcmp(name.c_str(), "Permittivity") == 0) return 7;
  return TPZMaterial::VariableIndex(name);
}
//! Number of variables associated with a given solution
int ScatteredField::NSolutionVariables(int var) const
{
  switch(var){
  case 0: //field (real part)
  case 1: //field (imag val)
  case 2: //field (abs val)
  case 3://deriv (real part)
  case 4://deriv (imag val)
  case 5://deriv (abs val)
  case 7://permittivity
    return this->Dimension();
  case 6:
    return 1;
  default:
    return TPZMaterial::NSolutionVariables(var);
  }
}
//! Computes the solution at an integration point
void ScatteredField::Solution(const TPZMaterialDataT<CSTATE> &data,
                              int var, TPZVec<CSTATE> &solout)
{

  TPZFNMatrix<9,CSTATE> er(3,3,0.);

  GetPermittivity(data.x, er);
  TPZManVector<CSTATE,3> epsvec = {er.g(0,0), er.g(1,1),er.g(2,2)};
  
  const TPZVec<CSTATE> &sol = data.sol[0];
  const TPZVec<CSTATE>&curlsol = data.curlsol[0];
  
  if(var == 6){
    solout[0] = this->Id();
    return;
  }
  const auto op = [var](auto &val){
    switch(var){
    case 0:
    case 3:
      return std::real(val);
    case 1:
    case 4:
      return std::imag(val);
    case 2:
    case 5:
    case 6:
    case 7:
      return std::abs(val);
    default:
      DebugStop();
      return std::real(val);
    }
  };

  const TPZVec<CSTATE> &val = var < 3 ? sol :
    (var == 7 ? epsvec : curlsol);

  for(auto x = 0; x < 3; x++){
    solout[x] = op(val[x]);
  }
}

#include <Electromagnetics/TPZCartesianPML.h>
template class TPZSingleSpaceCartesianPML<ScatteredField>;
