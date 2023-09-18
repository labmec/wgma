/**

 * @file AcousticModesBeta.h
 * @brief Header file for class AcousticModesBeta.\n
 */

#ifndef ACOUSTICMODESBETA_H
#define ACOUSTICMODESBETA_H


#include <materials/acousticmodesbase.hpp>
#include <TPZMatBase.h>
#include <TPZMatSingleSpace.h>
#include <TPZMatQuadraticEigenVal.h>

namespace wgma::materials{
    /**
     * @ingroup material
     * @brief This class implements the weak statement for the modal analysis of acoustic modes
     in waveguides supporting full anisotropy.
    */
    class  AcousticModesBeta :
        public TPZMatBase<CSTATE,
                          TPZMatSingleSpaceT<CSTATE>,
                          AcousticModesBase,
                          TPZMatQuadraticEigenVal>
    {
        using TBase = TPZMatBase<CSTATE,
                                 TPZMatSingleSpaceT<CSTATE>,
                                 AcousticModesBase,
                                 TPZMatQuadraticEigenVal>;
    public:

        /**
           @brief Constructor taking a few material parameters
           @param[in] id Material identifier.
           @param[in] mu Lame's parameter
           @param[in] lambda Lame's parameter
           @param[in] rho Density of the material
           @param[in] freq Frequency
           @param[in] scale Scale for geometric domain.
           @note the `scale` param might help with floating point arithmetics on really small domains.
        */
        AcousticModesBeta(int id, const STATE mu,
                          const STATE lambda, const STATE rho,
                          const STATE freq,
                          const REAL scale = 1.);
    
        explicit AcousticModesBeta(int id);

        AcousticModesBeta * NewMaterial() const override;
    
        std::string Name() const override { return "AcousticModesBeta"; }
    
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

        /** @name Quadratic EVP Methods */
        /** @{*/
        //! Set the material to compute matrix K
        void SetMatrixK() override;
        //! Set the material to compute matrix L
        void SetMatrixL() override;
        //! Set the material to compute matrix M
        void SetMatrixM() override;
        /**@{*/
        /**
           @name FrequencyMethods
           @{*/
        //! Get the frequency (Hz)
        inline const STATE &GetFrequency() const
        {return fFreq;}
        //! Set the frequency (Hz)
        inline void SetFrequency(const STATE &freq)
        { fFreq = freq;}
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
        //! Frequency being analysed (hz)
        STATE fFreq{11e+9};
        /** @name InternalGeneralisedAttributes */
        /** @{ */
        //! Pointer to the current Contribute function
        TContributeType fCurrentContribute{nullptr};
        //! Pointer to the current ContributeBC function
        TContributeBCType fCurrentContributeBC{nullptr};
        /** @} */
        AcousticModesBeta();//< Default constructor

        /** @name InternalContributeMethods */
        /** @{*/
        //! Contribution of the K Matrix
        void ContributeK(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                         TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef);
        //! Boundary contribution of the K Matrix
        void ContributeBCK(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                           TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                           TPZBndCondT<CSTATE> &bc);
        //! Contribution of the L Matrix
        void ContributeL(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                         TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef);
        //! Boundary contribution of the L Matrix
        void ContributeBCL(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                           TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                           TPZBndCondT<CSTATE> &bc) {}
        //! Contribution of the K Matrix
        void ContributeM(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                         TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef);
        //! Boundary contribution of the M Matrix
        void ContributeBCM(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                           TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                           TPZBndCondT<CSTATE> &bc) {}
        /** @{ */
    };
};
#endif
