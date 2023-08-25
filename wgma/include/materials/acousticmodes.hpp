/**

 * @file AcousticModes.h
 * @brief Header file for class AcousticModes.\n
 */

#ifndef ACOUSTICMODES_H
#define ACOUSTICMODES_H


#include <TPZMatBase.h>
#include <TPZMatSingleSpace.h>
#include <TPZMatQuadraticEigenVal.h>

namespace wgma::materials{
    /**
     * @ingroup material
     * @brief This class implements the weak statement for the modal analysis of acoustic modes
     in waveguides supporting full anisotropy.
    */
    class  AcousticModes :
        public TPZMatBase<CSTATE,
                          TPZMatSingleSpaceT<CSTATE>,
                          TPZMatQuadraticEigenVal>
    {
        using TBase = TPZMatBase<CSTATE,
                                 TPZMatSingleSpaceT<CSTATE>,
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
        AcousticModes(int id, const STATE mu,
                      const STATE lambda, const STATE rho,
                      const STATE freq,
                      const REAL scale = 1.);
    
        explicit AcousticModes(int id);

        AcousticModes * NewMaterial() const override;
    
        std::string Name() const override { return "AcousticModes"; }
    
        /** @brief Returns the integrable dimension of the material */
        int Dimension() const override {return 2;}

        [[nodiscard]] int NStateVariables() const override{return 3;}

        /**
           @name ParamMethods
           @{
        */
        //! Sets the wavelength being analysed
        inline void SetFrequency(STATE f){fFreq = f;} 
        inline void SetMu(STATE mu){fMu = mu;}
        inline void SetLambda(STATE lambda){fLambda = lambda;}
        inline void SetRho(STATE rho){fRho = rho;}
        //! Gets the current wavelength
        [[nodiscard]] inline STATE GetFrequency() const{ return fFreq;}
        [[nodiscard]] inline STATE GetMu() const{return fMu;}
        [[nodiscard]] inline STATE GetLambda() const{return fLambda;}
        [[nodiscard]] inline STATE GetRho() const{return fRho;}
        /**@}*/
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
        int VariableIndex(const std::string &name) const override;
        //! Number of variables associated with a given solution
        int NSolutionVariables(int var) const override;
        //! Computes the solution at an integration point
        void Solution(const TPZMaterialDataT<CSTATE> &data,
                      int var, TPZVec<CSTATE> &solout) override;
        //! Get the propagation constant used for post processing the solution
        inline const CSTATE &GetKz() const
        {return fKz;}
        //! Set the propagation constant used for post processing the solution
        inline void SetKz(const CSTATE &kz)
        { fKz = kz;}
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
    protected:
        void ComputeLxy(const TPZFMatrix<STATE> &dphi,
                        const TPZVec<REAL> &x,
                        TPZFMatrix<CSTATE> &lxy);
        void ComputeLz(const TPZFMatrix<STATE> &phi,
                       TPZFMatrix<CSTATE> &lz);
        void ComputeC(TPZFMatrix<CSTATE> &cmat);
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
        //! Lame's parameter
        STATE fLambda{1};
        //! Lame's parameter
        STATE fMu{1};
        //! Density
        STATE fRho{1};
        //! Frequency being analysed (hz)
        STATE fFreq{11e+9};
        //! Scale factor for the domain (helps with floating point arithmetic on small domains)
        const REAL fScaleFactor{1.};
        //!Fixes a propagation constant for printing the solution
        CSTATE fKz{1};

        /** @name InternalGeneralisedAttributes */
        /** @{ */
        //! Pointer to the current Contribute function
        TContributeType fCurrentContribute{nullptr};
        //! Pointer to the current ContributeBC function
        TContributeBCType fCurrentContributeBC{nullptr};
        /** @} */
        AcousticModes();//< Default constructor

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
