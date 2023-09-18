/**

 * @file AcousticModesBase.h
 * @brief Header file for class AcousticModesBase.\n
 */

#ifndef ACOUSTICMODESBASE_H
#define ACOUSTICMODESBASE_H

#include <pzreal.h>
#include <TPZMaterial.h>
template<class T>
class TPZMaterialDataT;
template<class T>
class TPZVec;
template<class T>
class TPZFMatrix;


namespace wgma::materials{

    //forward declare bc auxiliary class
    class AcousticModesBaseBC;
    /**
     * @ingroup material
     * @brief This class the common interface for materials that compute acoustic modes
     in waveguides
    */
    class  AcousticModesBase
    {
    public:
        using TInterfaceBC = AcousticModesBaseBC;

        /**
           @brief Constructor taking a few material parameters
           @param[in] id Material identifier.
           @param[in] mu Lame's parameter
           @param[in] lambda Lame's parameter
           @param[in] rho Density of the material
           @param[in] scale Scale for geometric domain.
           @note the `scale` param might help with floating point arithmetics on really small domains.
        */
        AcousticModesBase(int id, const STATE mu,
                          const STATE lambda, const STATE rho,
                          const REAL scale = 1.);
    
        explicit AcousticModesBase(int id);
    
        /** @brief Returns the integrable dimension of the material */
        int Dimension() const {return 2;}

        [[nodiscard]] int NStateVariables() const{return 3;}

        /**
           @name ParamMethods
           @{
        */
        inline void SetMu(STATE mu){fMu = mu;}
        inline void SetLambda(STATE lambda){fLambda = lambda;}
        inline void SetRho(STATE rho){fRho = rho;}
        [[nodiscard]] inline STATE GetMu() const{return fMu;}
        [[nodiscard]] inline STATE GetLambda() const{return fLambda;}
        [[nodiscard]] inline STATE GetRho() const{return fRho;}
        /**@}*/
        /**
           @name SolutionMethods
           @{*/
        /** @brief Variable index of a given solution.
            Available variables are:
            u_real
            u_abs
        */
        int VariableIndex(const std::string &name) const;
        //! Number of variables associated with a given solution
        int NSolutionVariables(int var) const;
        //! Computes the solution at an integration point
        void Solution(const TPZMaterialDataT<CSTATE> &data,
                      int var, TPZVec<CSTATE> &solout);
        /**@}*/
        
    protected:
        void ComputeLxy(const TPZFMatrix<STATE> &dphi,
                        const TPZVec<REAL> &x,
                        TPZFMatrix<CSTATE> &lxy);
        void ComputeLz(const TPZFMatrix<STATE> &phi,
                       TPZFMatrix<CSTATE> &lz);
        void ComputeC(TPZFMatrix<CSTATE> &cmat);
        
        //! Lame's parameter
        STATE fLambda{1};
        //! Lame's parameter
        STATE fMu{1};
        //! Density
        STATE fRho{1};
        //! Scale factor for the domain (helps with floating point arithmetic on small domains)
        const REAL fScaleFactor{1.};
        AcousticModesBase() = default;//< Default constructor
    };

    
    class AcousticModesBaseBC : public AcousticModesBase{
    protected:
        AcousticModesBase *fMatAcousticModesBase{nullptr};
        void SetMaterialImpl(TPZMaterial* mat){
            fMatAcousticModesBase = dynamic_cast<AcousticModesBase*>(mat);
        }
    };
};
#endif
