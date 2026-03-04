/*
 * @file scatteredfield.hpp
 * @brief Contains the ScatteredField class which computes the scattered field
 for a given background field, which must be loaded into the mesh
*/

#ifndef WGMASCATTEREDFIELD_H
#define WGMASCATTEREDFIELD_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMaterialData.h"
#include "pzvec.h"

namespace wgma::materials{
  /**
   * @brief Computes the periodic modes associated with the propagating
   and evanescent solutions of a given homogeneous medium with periodic BCs
  */
  class ScatteredField :
    public TPZMatBase<CSTATE,TPZMatSingleSpaceT<CSTATE>>{
    using TBase = TPZMatBase<CSTATE,TPZMatSingleSpaceT<CSTATE>>;	
  public:
    /**
     * @brief Class constructor 
     @param[in] id Material identifier.
     @param[in] er Relative permittivity.
     @param[in] ur Relative permeability.
     @param[in] scale Scale for geometric domain.
     @note the `scale` param might help with floating point arithmetics on really small domains.
    */
    ScatteredField(int id, const CSTATE er,
                   const CSTATE ur, const STATE wl,
                   const REAL scale);

    /**
     * @brief Class constructor 
     @param[in] id Material identifier.
     @param[in] er Relative permittivity.
     @param[in] ur Relative permeability.
     @param[in] scale Scale for geometric domain.
     @note the `scale` param might help with floating point arithmetics on really small domains.
    */
    ScatteredField(int id, const TPZFMatrix<CSTATE>& er, const TPZFMatrix<CSTATE>& ur, const STATE lambda,
                   const REAL scale = 1.);

    //! Sets the solution as a required data
    void FillDataRequirements(TPZMaterialData &data) const override{
      TBase::FillDataRequirements(data);
      data.fNeedsSol = m_has_sol;
    }

    /**
       @name ParamMethods
       @{
    */
    //! Sets the wavelength being analysed
    void SetWavelength(STATE wl) {m_wl = wl;}
    //! Gets the current wavelength
    [[nodiscard]] inline STATE GetWavelength() const{ return m_wl;}
    //! Sets the permeability of the material
    void SetPermeability(CSTATE ur);
    //! Sets the permeability of the material
    void SetPermeability(const TPZFMatrix<CSTATE> &ur);
    //! Gets the permeability of the material
    inline virtual void GetPermeability([[maybe_unused]] const TPZVec<REAL> &x,
                                        TPZFMatrix<CSTATE> &ur) const
    {ur = m_ur;}
    //! Sets the permittivity of the material
    void SetPermittivity(CSTATE er);
    //! Sets the permittivity of the material
    void SetPermittivity(const TPZFMatrix<CSTATE> &er);
    //! Gets the permittivity of the material
    inline virtual void GetPermittivity([[maybe_unused]] const TPZVec<REAL> &x,
                                        TPZFMatrix<CSTATE> &er) const
    {er = m_er;}

    //! Sets the permittivity of the material
    void SetBackgroundPermittivity(CSTATE er){
        m_back_er.Identity();
        m_back_er *= er;
    }
    //! Sets the permittivity of the material
    void SetBackgroundPermittivity(const TPZFMatrix<CSTATE> &er){
        m_back_er = er;
    }
    //! Gets the permittivity of the material
    void GetBackgroundPermittivity(TPZFMatrix<CSTATE> &er) const{
        er = m_back_er;
    }
      
    /**@}*/

    std::string Name() const override { return "ScatteredField"; }
	
    /** @brief Solution indices of post-processing */
    int Dimension() const  override { return this->m_dim; }

    int NStateVariables() const override { return 1; }

  
    void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
            TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override;
  
    void ContributeBC(const TPZMaterialDataT<CSTATE> &data, REAL weight,
              TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
              TPZBndCondT<CSTATE> &bc) override {}
    /** @brief To create another material of the same type */
    ScatteredField * NewMaterial() const override
    {
      return new ScatteredField(*this);
    }
    /**
       @name SolutionMethods
       @{*/
    /** @brief Variable index of a given solution.
      Possibilities are:
      -Field_real
      -Field_imag
      -Field_abs
      -Deriv_real
      -Deriv_imag
      -Deriv_abs
      -Material
      -Permittivity
    */
    int VariableIndex(const std::string &name) const override;
    //! Number of variables associated with a given solution
    int NSolutionVariables(int var) const override;
    //! Computes the solution at an integration point
    void Solution(const TPZMaterialDataT<CSTATE> &data,
                  int var, TPZVec<CSTATE> &solout) override;

    //! Gets dimensions for solution variable
    void GetSolDimensions(uint64_t &u_len,
                          uint64_t &du_row,
                          uint64_t &du_col) const override{
      u_len = m_dim;
      du_row = m_dim;
      du_col = 1;
    }
    //! Sets a multiplicative scale factor for .vtk fields
    void SetScaleVTK(REAL scale){m_scale_post = scale;}
    //! Gets the multiplicative scale factor for .vtk fields
    REAL GetScaleVTK() const {return m_scale_post ;}
    //! Sets whether solution should be computed in this specific instance
    void SetComputeSol(bool val){m_has_sol = val;}
    /**@}*/
    virtual int ClassId() const override;
  protected:
    //! problem dimension
    constexpr static int m_dim{3};
    ScatteredField() = default;
    //! Relative magnetic permeability (xx, yy, zz)
    TPZFNMatrix<9,CSTATE> m_ur{{1.,0,0},{0,1,0},{0,0,1}};
    //! Relative electric permittivity (xx, yy, zz)
    TPZFNMatrix<9,CSTATE> m_er{{1.,0,0},{0,1,0},{0,0,1}};
    //! Relative background electric permittivity (xx, yy, zz)
    TPZFNMatrix<9,CSTATE> m_back_er{{1.,0,0},{0,1,0},{0,0,1}};
    //! Wavelength being analysed
    STATE m_wl{0};
    //! Scale factor for the domain (helps with floating point arithmetic on small domains)
    const REAL m_scale{1.};
    //! Scale factor for post processing vtk files (ONLY)
    REAL m_scale_post{1.};
    //! Whether it has a background field solution or not (ex: PML domain)
    bool m_has_sol{true};
  };
};
#endif
