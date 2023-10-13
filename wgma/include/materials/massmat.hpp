/**
 * @file MassMatrix.h
 * @brief Contains the MassMatrix class which implements a mass matrix
 */

#ifndef TPZMASSMAT_H
#define TPZMASSMAT_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"

namespace wgma::materials{
/**
 * @brief Computes mass matrix phi_i phi_j. rhs is phi_i*phi_i
 */
template<class TVar=STATE>
class MassMatrix :
    public TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>>{
	using TBase = TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>>;	
public:
    //! Default constructor
    MassMatrix() = default;
    /**
	 * @brief Class constructor 
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 */
	MassMatrix(int id, int dim, int nstate=1);
	/**
	 * @brief Class constructor setting a default constant solution
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 * @param sol constant solution. Ignored if a forcing function is set.
	 */
	MassMatrix(int id, int dim, int nstate, const TPZVec<TVar> &sol);

    std::string Name() const override { return "MassMatrix"; }
	
	/** @brief Solution indices of post-processing */
	enum ESolutionVars { ENone = 0, ESolution = 1 , EDerivative = 2};
	
    /**
     * @brief Set a scale factor for the stiffness matrix and right hand side
     * the default value of the scale factor is 1
     */
    void SetScaleFactor(REAL scale)
    {fScale = scale;}
    
    [[nodiscard]] REAL ScaleFactor() const
    {return fScale;}

	int Dimension() const  override { return this->fDim; }
    
    /** @brief Sets problem dimension */
    virtual void SetDimension(int dim) { this->fDim = dim; }
	
    int NStateVariables() const override { return this->fNStateVars; }
    void SetNStateVariables(int nstate) { this->fNStateVars = nstate; }
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
	
	/** @brief Number of state variables */
	int fNStateVars{1};
	
	/** @brief Constant solution vector. Ignored if forcing function is set. */
	TPZVec<TVar> fSol;
    
    /** @brief Scale factor applied to the stiffness matrix and right hand side */
    REAL fScale{1};
};


extern template class MassMatrix<STATE>;
extern template class MassMatrix<CSTATE>;
};
#endif
