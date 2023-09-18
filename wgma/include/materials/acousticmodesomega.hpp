/**

 * @file AcousticModesOmega.h
 * @brief Header file for class AcousticModesOmega.\n
 */

#ifndef ACOUSTICMODESBETA_H
#define ACOUSTICMODESBETA_H


#include <materials/acousticmodesbase.hpp>
#include <TPZMatBase.h>
#include <TPZMatSingleSpace.h>
#include <TPZMatGeneralisedEigenVal.h>

namespace wgma::materials{
    /**
     * @ingroup material
     * @brief This class implements the weak statement for the modal analysis of acoustic modes.
    */
    class  AcousticModesOmega :
        public TPZMatBase<CSTATE,
                          TPZMatSingleSpaceT<CSTATE>,
                          AcousticModesBase,
                          TPZMatGeneralisedEigenVal>
    {
        using TBase = TPZMatBase<CSTATE,
                                 TPZMatSingleSpaceT<CSTATE>,
                                 AcousticModesBase,
                                 TPZMatGeneralisedEigenVal>;
    public:

        /**
           @brief Constructor taking a few material parameters
           @param[in] id Material identifier.
           @param[in] mu Lame's parameter
           @param[in] lambda Lame's parameter
           @param[in] rho Density of the material
           @param[in] beta Propagation constant
           @param[in] scale Scale for geometric domain.
           @note the `scale` param might help with floating point arithmetics on really small domains.
        */
        AcousticModesOmega(int id, const STATE mu,
                          const STATE lambda, const STATE rho,
                          const CSTATE beta,
                          const REAL scale = 1.);
    
        explicit AcousticModesOmega(int id);

        AcousticModesOmega * NewMaterial() const override;
    
        std::string Name() const override { return "AcousticModesOmega"; }
    
        /** @brief Returns the integrable dimension of the material */
        int Dimension() const override {return AcousticModesBase::Dimension();}

        [[nodiscard]] int NStateVariables() const override{
            return AcousticModesBase::NStateVariables();
        }
        /**
           @name ContributeMethods
           @{
        */
        void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                        TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override;
    
        void ContributeBC(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                          TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                          TPZBndCondT<CSTATE> &bc) override;
        /**@}*/
        /**
           @name SolutionMethods
           @{*/
        /** @brief Variable index of a given solution.
            Available variables are:
            u_real
            u_abs
        */
        
        int VariableIndex(const std::string &name) const override
        { return AcousticModesBase::VariableIndex(name);}
        //! Number of variables associated with a given solution
        int NSolutionVariables(int var) const override
        { return AcousticModesBase::NSolutionVariables(var);}
        //! Computes the solution at an integration point
        void Solution(const TPZMaterialDataT<CSTATE> &data,
                      int var, TPZVec<CSTATE> &solout) override
        { AcousticModesBase::Solution(data,var,solout);}
        /**@}*/

        /** @name Generalised EVP Methods */
        /** @{*/
        //! Set the material to compute matrix A
        void SetMatrixA() override;
        //! Set the material to compute matrix B
        void SetMatrixB() override;
        /**@{*/
        /**
           @name BetaMethods
           @{*/
        //! Get the propagation constant
        inline const CSTATE &GetBeta() const
        {return fBeta;}
        //! Set the frequency (Hz)
        inline void SetBeta(const CSTATE &beta)
        { fBeta = beta;}
        /**@}*/
    protected:
        /** @name InternalGeneralisedTypes */
        /** @{ */
        using TContributeType =
            std::function<
            void (const TPZMaterialDataT<CSTATE> &data,
                  REAL weight,TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)>;
        using TContributeBCType =
            std::function<
            void (const TPZMaterialDataT<CSTATE> &data,
                  REAL weight,
                  TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                  TPZBndCondT<CSTATE> &bc)>;
        /** @} */
        //! Propagation constant
        CSTATE fBeta{0};
        /** @name InternalGeneralisedAttributes */
        /** @{ */
        //! Pointer to the current Contribute function
        TContributeType fCurrentContribute{nullptr};
        //! Pointer to the current ContributeBC function
        TContributeBCType fCurrentContributeBC{nullptr};
        /** @} */
        AcousticModesOmega();//< Default constructor

        /** @name InternalContributeMethods */
        /** @{*/
        //! Contribution of the A Matrix
        void ContributeA(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                         TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef);
        //! Boundary contribution of the A Matrix
        void ContributeBCA(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                           TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                           TPZBndCondT<CSTATE> &bc);
        //! Contribution of the B Matrix
        void ContributeB(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                         TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef);
        //! Boundary contribution of the B Matrix
        void ContributeBCB(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                           TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                           TPZBndCondT<CSTATE> &bc) {}
        /** @{ */
    };
};
#endif
