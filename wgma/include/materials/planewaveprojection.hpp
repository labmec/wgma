/*
 * @file planewaveprojection.hpp
 * @brief Contains the PlanewaveProjection class which projects a linearly polarised
 plane wave solution onto a Hcurl 2d mesh
 Useful for analysis of meta surfaces
 
*/

#ifndef WAVEPROJ_H
#define WAVEPROJ_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"

namespace wgma::materials{
  /**
   * @brief Projects a plane wave solution onto a HCurl 2d mesh.
   Angle with respect to x-axis is set by theta
   */
  template<class TVar=STATE>
  class PlanewaveProjection :
    public TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>>{
    using TBase = TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>>;  
  public:
    //! Default constructor
    PlanewaveProjection() = default;
    /**
     * @brief Class constructor 
     * @param id material id
     * @param theta angle in radians w.r.t. x axis
     */
    PlanewaveProjection(int id, REAL theta);

    /** 
     * \brief Fill the TPZMaterialData parameter.
     * The needed parameters for the calculations of the 
     * Contribute methods should be set.By default, 
     * all requirements are considered to be false.
     * \param The TPZMaterialData instance to be filled
     */
    void FillDataRequirements(TPZMaterialData &data) const override final{}

    std::string Name() const override final{ return "PlanewaveProjection"; }
  
    /** @brief Solution indices of post-processing */
    enum ESolutionVars { ENone = 0, ESol = 1, EDerivative = 2};
  
    /**
     * @brief Set a scale factor for the stiffness matrix and right hand side
     * the default value of the scale factor is 1
     */
    void SetScaleFactor(REAL scale)
    {m_scale = scale;}
    
    [[nodiscard]] REAL ScaleFactor() const
    {return m_scale;}

    /**@brief Sets angle with respect to x-axis(in radian)*/
    void SetTheta(REAL theta)
    {m_theta = theta;}
    
    [[nodiscard]] REAL Theta() const
    {return m_theta;}

    int Dimension() const  override final{ return this->fDim; }
  
    int NStateVariables() const override final{ return 1; }
    
    void Contribute(const TPZMaterialDataT<TVar> &data,
                    REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override final;

    void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override final{}
    /** @brief To create another material of the same type */
    TPZMaterial * NewMaterial() const override final;
  
    /** @brief It returns the variable index associated with the name */
    int VariableIndex(const std::string &name) const override final;
  
    int NSolutionVariables(int var) const override final;

    void Solution(const TPZMaterialDataT<TVar> &data,
                  int var, TPZVec<TVar> &solOut) override final;
    
    void GetSolDimensions(uint64_t &u_len,
                          uint64_t &du_row,
                          uint64_t &du_col) const override final
    {
      u_len = 3; du_row = 1; du_col = 1;
    }
    
    virtual int ClassId() const override final;
  protected:
    /** @brief Problem dimension */
    static constexpr int fDim{2};
    /** @brief Angle of solution with respect to x-axis(in radians)*/
    REAL m_theta{0};
    /** @brief Scale factor applied to the stiffness matrix and right hand side */
    REAL m_scale{1};
  };


  extern template class PlanewaveProjection<STATE>;
  extern template class PlanewaveProjection<CSTATE>;
};
#endif
