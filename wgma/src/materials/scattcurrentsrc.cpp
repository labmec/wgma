#include "materials/scattcurrentsrc.hpp"
#include <TPZMaterialDataT.h>

using namespace std::complex_literals;
using namespace wgma::materials;

//! Unique identifier for serialization purposes
int CurrentSource3D::ClassId() const
{
  return Hash("CurrentSource3D");
}

void CurrentSource3D::Write(TPZStream &buf, int withclassid) const
{
  buf.Write(j);
}
  //! Read from stream(serialization method)
void CurrentSource3D::Read(TPZStream &buf, void *context)
{
  j.Read(buf, context);
}

void CurrentSource3D::Print(std::ostream &out) const
{
  j.Print("j",out);
  out << "x:";
  for(auto xi : x){
    out << ' '<< xi;
  }
  out <<'\n';
}


void ScattCurrentSrc::FillDataRequirements(
  TPZMaterialData &data) const
{
  data.fNeedsNormal = true;
}

//! Contribution to the integration point
void ScattCurrentSrc::Contribute(const TPZMaterialDataT<CSTATE> &data,
                                     REAL weight, TPZFMatrix<CSTATE> &ef)
{
  //index of integration point
  const int gp_index = data.intGlobPtIndex;
  const auto &mem_item = this->MemItem(gp_index);
  
  const int nshape = data.phi.Rows();
  const auto &phi_real = data.phi;
  TPZFNMatrix<3000,CSTATE> phi(3,nshape);

  //we cannot use phi_real_ptr because its actually transposed
  CSTATE *phi_ptr = phi.Elem();
  const int sz = 3*nshape;
  for(int i = 0; i < sz; i++){
    //g instead of GetVal will be most likely inlined
    *phi_ptr++ = phi_real.g(i/3, i%3);
  }

  const auto &j = mem_item.j;
  constexpr int conjtransp{2}, notransp{0};
  ef.AddContribution(0, 0, phi, conjtransp, j, notransp, weight);
}

ScattCurrentSrc * ScattCurrentSrc::NewMaterial() const
{
  return new ScattCurrentSrc(*this);
}
//! Unique identifier for serialization purposes
int ScattCurrentSrc::ClassId() const
{
  return
    Hash("ScattCurrentSrc")
    ^
    TPZScattering::ClassId() << 1
    ^
    TPZMatWithMem<CurrentSource3D>::ClassId() << 2;
}
//! Write to stream(serialization method)
void ScattCurrentSrc::Write(TPZStream &buf, int withclassid) const
{
  TPZScattering::Write(buf,withclassid);
  TPZMatWithMem<CurrentSource3D>::Write(buf, withclassid);
}
//! Read from stream(serialization method)
void ScattCurrentSrc::Read(TPZStream &buf, void *context)
{
  TPZScattering::Read(buf,context);
  TPZMatWithMem<CurrentSource3D>::Read(buf, context);
}

#include "Electromagnetics/TPZCartesianPML.h"
template class TPZSingleSpaceCartesianPML<ScattCurrentSrc>;