/*
 * @file planewavesolutions.hpp
 * @brief Contains the PlaneWaveSolutions class which projects the plane wave
 solutions of a periodic homogeneous waveguide in an hcurl x h1 approx space
*/

#ifndef TPZPLANEWAVESOLUTIONS_H
#define TPZPLANEWAVESOLUTIONS_H

#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatLoadCases.h"
#include "Electromagnetics/TPZWgma.h"//for hcurl and h1 indexes
#include "pzvec.h"

namespace wgma::materials{
  /**
   * @brief Computes the periodic modes associated with the propagating
   and evanescent solutions of a given homogeneous medium with periodic BCs
  */
  class PlaneWaveSolutions :
        public TPZMatBase<CSTATE,TPZMatCombinedSpacesT<CSTATE>,
        TPZMatLoadCases<CSTATE>>{
      using TBase = TPZMatBase<CSTATE,TPZMatCombinedSpacesT<CSTATE>,
                               TPZMatLoadCases<CSTATE>>;	
  public:
    //! Default constructor (should not be called)
    PlaneWaveSolutions() {DebugStop();}
    /**
     * @brief Class constructor 
     * @param id material id
     * @param wavelength operational wavelength
     * @param ref_index refractive index
     * @param lx domain's length in x direction
     * @param ly domain's length in y direction
     * @param max_k maximum value of k for which we compute the periodic modes
     */
    PlaneWaveSolutions(int id, STATE wavelength, CSTATE ref_index,
               REAL lx, REAL ly, int max_k);

    //! Sets new wavelength and recompute beta values
    inline void SetWavelength(STATE wl){
      m_wl = wl;
      ComputeBeta();
    }

    //! Returns the minimum number of load cases (rhs columns).
    [[nodiscard]] int MinimumNumberofLoadCases() const override;

    //! Returns DUPLICATED beta values (Ex, Ey, Ex, Ey)
    void GetBeta(TPZVec<CSTATE> &beta);

    std::string Name() const override { return "PlaneWaveSolutions"; }
	
    /** @brief Solution indices of post-processing */
    int Dimension() const  override { return this->m_dim; }

    int NStateVariables() const override { return 1; }

    int IntegrationRuleOrder(const TPZVec<int> &elpmax_ord) const override{
        return TBase::IntegrationRuleOrder(elpmax_ord)+4;
    }
  
    void Contribute(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
            TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override;
  
    void ContributeBC(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
              TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
              TPZBndCondT<CSTATE> &bc) override {}
    /** @brief To create another material of the same type */
    PlaneWaveSolutions * NewMaterial() const override;
    /**
       @name SolutionMethods
       @{*/
    /** @brief Variable index of a given solution.
      Possibilities are:
      -Et_real
      -Ez_real
      -Et_abs
      -Ez_abs
    */
    int VariableIndex(const std::string &name) const override;
    //! Number of variables associated with a given solution
    int NSolutionVariables(int var) const override;
    //! Computes the solution at an integration point
    void Solution(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
            int var, TPZVec<CSTATE> &solout) override;

    /**@}*/
    virtual int ClassId() const override;
  protected:
    //Updates the beta,kx and ky vector after a new wavelength has been set
    void ComputeBeta();
    //! problem dimension
    constexpr static int m_dim{2};
    //! maximum value of k taken in periodic solutions
    int m_max_k{-1};
    //! domain length in x direction
    STATE m_lx{0};
    //! domain length in y direction
    STATE m_ly{0};
    //! operational wavelength
    STATE m_wl{0};
    //! refractive index of homogeneous medium
    CSTATE m_ref_index{0};
    ///all data structures belowhave dim nsol/2 
    //! vector of kx (transverse wavenumber in x direction)
    TPZVec<STATE> m_kx;
    //! vector of ky (transverse wavenumber in y direction)
    TPZVec<STATE> m_ky;
    //! vector of beta
    TPZVec<CSTATE> m_beta;

    //! hcurl mesh index
    constexpr static int m_hcurl_index = TPZWgma::HCurlIndex();
    //! h1 mesh index
    constexpr static int m_h1_index = TPZWgma::H1Index();
  
  };
};
#endif
