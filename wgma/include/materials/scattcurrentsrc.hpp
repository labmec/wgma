/**

 * @file volcurrentsource.hpp
 * @brief Header file for class ScattCurrentSource, that implements a volumetric current source
 for electromagnetic wave propagation problems.
 The material will simply compute F\cdotE^*, where
 F is the source term stored in memory and E^* the conjugated basis function.\n
*/

#ifndef SCATTCURRENTSRC_H
#define SCATTCURRENTSRC_H

#include <Electromagnetics/TPZScattering.h>
#include <TPZMatWithMem.h>


namespace wgma::materials{

  //! Data to be stored at each integration point
  struct CurrentSource3D{
    TPZFNMatrix<3,CSTATE> j = {{0},{0},{0}};//< source term at point
    TPZManVector<REAL,3> x = {0,0,0};//<for debugging purposes
    //! Unique identifier for serialization purposes
    [[nodiscard]] int ClassId() const;
    //! Write to stream(serialization method)
    void Write(TPZStream &buf, int withclassid) const;
    //! Read from stream(serialization method)
    void Read(TPZStream &buf, void *context);
    //! Print contents(debugging method)
    void Print(std::ostream &out) const;
  };

  inline std::ostream& operator<<( std::ostream& out, const CurrentSource3D& t ){
    t.Print(out);
    return out;
  }
  //! Implements a volumetric current source for electromagnetic wave propagation problems
  class ScattCurrentSrc : public TPZScattering,
                               public TPZMatWithMem<CurrentSource3D>{
  public:
    //! All constructors from base class shall be available
    using TPZScattering::TPZScattering;
    //! Contribution to the MATRIX ONLY at the integration point
    void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override {
      TPZScattering::Contribute(data,weight,ek,ef);
    }
    
    //! Contribution to the rhs at the integration point
    void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                    TPZFMatrix<CSTATE> &ef) override;
  
    void FillDataRequirements(TPZMaterialData &data) const override;
    //! Returns the integrable dimension of the material
    int Dimension() const override {return 3;}
    //! Creates a copy of this instance
    ScattCurrentSrc * NewMaterial() const override;
    //! Returns name of the class
    std::string Name() const override { return "ScattCurrentSrc"; }
    //! Unique identifier for serialization purposes
    [[nodiscard]] int ClassId() const override;
    //! Write to stream(serialization method)
    void Write(TPZStream &buf, int withclassid) const override;
    //! Read from stream(serialization method)
    void Read(TPZStream &buf, void *context) override;
  };
};
#endif
