#include "materials/planewavesolutions.hpp"

#include <TPZMaterialDataT.h>

#include <numeric>
using namespace wgma::materials;
using namespace std::complex_literals;


PlaneWaveSolutions::PlaneWaveSolutions(int id, STATE wl, CSTATE ref_index,
                                       REAL lx, REAL ly, int max_k) :
  TPZRegisterClassId(&PlaneWaveSolutions::ClassId),
  TBase(id), m_wl(wl), m_ref_index(ref_index), m_max_k(max_k),
  m_lx(lx), m_ly(ly)
{
  ComputeBeta();
  const int sz = (2*m_max_k+1)*(2*m_max_k+1);
  SetNumLoadCases(2*sz);
}

int PlaneWaveSolutions::MinimumNumberofLoadCases() const{
  const int sz = (2*m_max_k+1)*(2*m_max_k+1);
  return 2*sz;
}

void PlaneWaveSolutions::ComputeBeta()
{
  
  const CSTATE wavenumber{2*M_PI*m_ref_index/m_wl};
  const CSTATE k2{wavenumber*wavenumber};

  const int max_k = m_max_k;
  const int sz = (max_k*2+1)*(max_k*2+1);
  //temporary vectors to keep the code simple
  TPZVec<STATE> kx(sz,0), ky(sz,0);
  TPZVec<CSTATE> beta(sz,0);

  auto is_propagating = [](CSTATE b){
    return std::abs(b.real()) > std::abs(b.imag());
  };
  
  int count{0};
  for(int i = -max_k; i <= max_k; i++){
    for(int j = -max_k; j <= max_k; j++){
      kx[count] = i*2*M_PI/m_lx;
      ky[count] = j*2*M_PI/m_ly;
      const auto kt2 = kx[count]*kx[count]+ky[count]*ky[count];
      CSTATE betaval = std::sqrt((CSTATE)(k2 - kt2));
      //we ensure that beta has negative imag part for evanescent modes
      if(!is_propagating(betaval) && betaval.imag() > 0){
        betaval = -betaval;
      }
      beta[count] = betaval;
      count++;
    }
  }

  //now we sort beta as we want and keep the indices

  // initialize original index locations
  std::vector<size_t> idx(beta.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
                   [&beta, &is_propagating](size_t i1, size_t i2) {
                     if (is_propagating(beta[i1])){
                       //propagating beta
                       if (is_propagating(beta[i2])){
                         const auto b1 = beta[i1].real();
                         const auto b2 = beta[i2].real();
                         return b1 > b2;
                       }else{
                         return true;
                       }
                     }else{
                       //evanescent beta
                       if (is_propagating(beta[i2])){
                         return false;
                       }else{
                         const auto b1 = beta[i1].imag();
                         const auto b2 = beta[i2].imag();
                         return b1 > b2;
                       }

                     }
                   });

  //now we sort each array according to idx
  m_kx.resize(sz);
  m_ky.resize(sz);
  m_beta.resize(sz);
  for(auto i = 0; i < sz; i++){
    m_kx[i] = kx[idx[i]];
  }
  for(auto i = 0; i < sz; i++){
    m_ky[i] = ky[idx[i]];
  }
  for(auto i = 0; i < sz; i++){
    m_beta[i] = beta[idx[i]];
  } 
  
}

void PlaneWaveSolutions::GetBeta(TPZVec<CSTATE> &beta)
{
  const int nsol = (2*m_max_k+1)*(2*m_max_k+1);
  const int nsolvec = 2*nsol;
  beta.resize(nsolvec);
  for(int i = 0; i < nsol; i++){
    beta[2*i] = m_beta[i];
    beta[2*i+1] = m_beta[i];
  }
}

void
PlaneWaveSolutions::Contribute(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
                               REAL weight,TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)
{
  const int nsol = (2*m_max_k+1)*(2*m_max_k+1);

  const auto &phi_hcurl_real = datavec[m_hcurl_index].phi;
  const auto &phi_h1_real = datavec[m_h1_index].phi;
  const int nhcurl  = phi_hcurl_real.Rows();
  const int nh1  = phi_h1_real.Rows();
  //making complex version of phi hcurl
  TPZFNMatrix<200,CSTATE> phi_hcurl(2,nhcurl,0.);
  for(int i = 0; i < nhcurl; i++){
    for(int x = 0; x < 2; x++){
      phi_hcurl.PutVal(x,i,phi_hcurl_real.GetVal(i,x));
    }
  }
  //making complex version of phi h1
  TPZFNMatrix<100,CSTATE> phi_h1(1,nh1,0.);
  for(int i = 0; i < nh1; i++){
    phi_h1.PutVal(0,i,phi_h1_real.GetVal(i,0));
  }

  /*****************ACTUAL COMPUTATION OF CONTRIBUTION****************/
  const int firsthcurl = m_hcurl_index * nh1;
  const int firsth1 = m_h1_index * nhcurl;

  constexpr int transp{1};
  constexpr int no_transp{0};
  ek.AddContribution(firsthcurl,firsthcurl,phi_hcurl,transp,phi_hcurl,no_transp, weight);
  ek.AddContribution(firsth1,firsth1,phi_h1,transp,phi_h1,no_transp, weight);
  //now we must compute the actual modes, both Et and Ez
  const int nsolvec = 2*nsol;
  TPZFNMatrix<2000,CSTATE> sol_et(2,nsolvec,0.);
  TPZFNMatrix<1000,CSTATE> sol_ez(1,nsolvec,0.);
  const auto x = datavec[0].x[0];
  const auto y = datavec[0].x[1];
  for(int im = 0; im < nsol; im++){
    const CSTATE kx = m_kx[im];
    const CSTATE ky = m_ky[im];
    const CSTATE beta = m_beta[im];

    //solutions with et in the x direction
    const auto eval = std::exp(-1i*(kx*x+ky*y));
    sol_et.PutVal(0,2*im,eval);
    sol_ez.PutVal(0,2*im,-eval*kx/beta);
    //solutions with et in the y direction
    sol_et.PutVal(1,2*im+1,eval);
    sol_ez.PutVal(0,2*im+1,-eval*ky/beta);
  }
  ef.AddContribution(firsthcurl,0,phi_hcurl,transp,sol_et,no_transp,weight);
  ef.AddContribution(firsth1,0,phi_h1,transp,sol_ez,no_transp,weight);
}


int
PlaneWaveSolutions::VariableIndex(const std::string &name) const
{
  if( strcmp(name.c_str(), "Et_real") == 0) return 0;
  if( strcmp(name.c_str(), "Ez_real") == 0) return 1;
  if( strcmp(name.c_str(), "Et_abs") == 0) return 2;
  if( strcmp(name.c_str(), "Ez_abs") == 0) return 3;
  DebugStop();
  return 1;
}

int
PlaneWaveSolutions::NSolutionVariables(int var) const
{
  switch (var) {
  case 0: //Et_real
    return 3;
  case 1://Ez_real
    return 1;
  case 2: //Et_abs
    return 3;
  case 3://Ez_abs
    return 1;
  default:
    DebugStop();
    break;
  }
  return 1;
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void
PlaneWaveSolutions::Solution(
  const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
  int var,
  TPZVec<CSTATE> &solout)
{
  TPZManVector<CSTATE,3> et(3,0.);
  TPZManVector<CSTATE,1> ez(1,0.);

  const auto idx = this->fPostProcIndex;
  et = datavec[ m_hcurl_index ].sol[idx];
  ez = datavec[ m_h1_index ].sol[idx];

  switch (var) {
  case 0:{//et_real
    for (int i = 0; i < et.size(); ++i) {
      et[i] = std::real(et[i]);
    }
    solout = et;
    break;
  }
  case 1:{//ez_real
    for (int i = 0; i < ez.size(); ++i) {
      ez[i] = std::real(ez[i]);
    }
    solout = ez;
    break;
  }
  case 2:{//et_abs
    for (int i = 0; i < et.size(); ++i) {
      et[i] = std::abs(et[i]);
    }
    solout = et;
    break;
  }
  case 3:{//ez_abs
    for (int i = 0; i < ez.size(); ++i) {
      ez[i] = std::abs(ez[i]);
    }
    solout = ez;
    break;
  }
  default:
    DebugStop();
    break;
  }
}

PlaneWaveSolutions* PlaneWaveSolutions::NewMaterial() const{
  return new PlaneWaveSolutions();
}

int PlaneWaveSolutions::ClassId() const{
  return Hash("PlaneWaveSolutions") ^ TBase::ClassId() << 1;
}