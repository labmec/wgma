/*
 * @file solutionprojection.hpp
 * @brief Contains the SolutionProjection class which projects a solution in an approx space
 */

#ifndef TPZSOLPROJ_H
#define TPZSOLPROJ_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"

namespace wgma::materials{
  /**
   * @brief Projects a given solution in an approx space
   */
  template<class TVar=STATE>
  class SolutionProjection :
    public TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>>{
    using TBase = TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>>;	
  public:
    //! Default constructor
    SolutionProjection() = default;
    /**
     * @brief Class constructor 
     * @param id material id
     * @param dim problem dimension
     * @param soldim dimension of solution (1 for h1, 3 for hcurl)
     */
    SolutionProjection(int id, int dim, int soldim);

    /** 
	 * \brief Fill the TPZMaterialData parameter.
     * The needed parameters for the calculations of the 
     * Contribute methods should be set.By default, 
     * all requirements are considered to be false.
     * \param The TPZMaterialData instance to be filled
     */
    virtual void FillDataRequirements(TPZMaterialData &data) const override;

    std::string Name() const override { return "SolutionProjection"; }
	
    /** @brief Solution indices of post-processing */
    enum ESolutionVars { ENone = 0, ESolReal = 1 , ESolImag = 2, ESolAbs = 3,
                         EDerivative = 3};
	
    /**
     * @brief Set a scale factor for the stiffness matrix and right hand side
     * the default value of the scale factor is 1
     */
    void SetScaleFactor(REAL scale)
    {fScale = scale;}
    
    [[nodiscard]] REAL ScaleFactor() const
    {return fScale;}

    /**
     * @brief Sets dimension of solution (1 for scalar, 3 for vector)
     */
    void SetSolDimension(REAL scale)
    {fScale = scale;}
    
    [[nodiscard]] REAL SolDimension() const
    {return fScale;}

    int Dimension() const  override { return this->fDim; }
    
    /** @brief Sets problem dimension */
    virtual void SetDimension(int dim) { this->fDim = dim; }
	
    int NStateVariables() const override { return 1; }
    
    void Contribute(const TPZMaterialDataT<TVar> &data,
                    REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override;
    /** @brief To create another material of the same type */
    TPZMaterial * NewMaterial() const override;
	
    /** @brief It returns the variable index associated with the name */
    int VariableIndex(const std::string &name) const override;
	
    int NSolutionVariables(int var) const override;

    void Solution(const TPZMaterialDataT<TVar> &data,
                  int var, TPZVec<TVar> &solOut) override;
    
    void GetSolDimensions(uint64_t &u_len,
                          uint64_t &du_row,
                          uint64_t &du_col) const override;
    
    virtual int ClassId() const override;
  protected:
    /** @brief Problem dimension */
    int fDim;
    /** @brief Scale factor applied to the stiffness matrix and right hand side */
    REAL fScale{1};
    /** @brief Dimension of solution*/
    int fSolDim{-1};
  };


  extern template class SolutionProjection<STATE>;
  extern template class SolutionProjection<CSTATE>;
};
#endif
