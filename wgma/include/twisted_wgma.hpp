#ifndef _TWISTED_WGMA_H_
#define _TWISTED_WGMA_H_

#include <Electromagnetics/TPZAnisoWgma.h>
#include <Electromagnetics/TPZCylindricalPML.h>

// #define CUSTOMPML

namespace wgma::materials{
    class TwistedWgma : public TPZAnisoWgma{
    public:
        //! Creates a material based on an anisotropic mat
        TwistedWgma(const TPZAnisoWgma &cp) : TPZAnisoWgma(cp), m_alpha{1}{}
        //! Gets the equivalent permeability of the material
        void GetPermeability([[maybe_unused]] const TPZVec<REAL> &x,
                             TPZFMatrix<CSTATE> &ur) const override;
        //! Gets the equivalent permittivity of the material
        void GetPermittivity([[maybe_unused]] const TPZVec<REAL> &x,
                             TPZFMatrix<CSTATE> &er) const override;

        //! Sets torsion parameter
        void SetAlpha(const STATE alpha) {m_alpha = alpha;}
        //! Gets torsion parameter
        [[nodiscard]] STATE GetAlpha() const {return m_alpha;}

        [[nodiscard]] int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override{
            return TPZAnisoWgma::IntegrationRuleOrder(elPMaxOrder)+2;
        }
        
        void Solution(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
                      int var, TPZVec<CSTATE> &solout) override;
    protected:
        void TransformationMatrix(TPZFMatrix<CSTATE> &mat, const TPZVec<REAL> &x) const;
        STATE m_alpha{1};
    };

#ifdef CUSTOMPML
    class TwistedWgmaPML : public TPZCombinedSpacesCylindricalPML<TwistedWgma>{
    public:
        using TPZCombinedSpacesCylindricalPML<TwistedWgma>::TPZCombinedSpacesCylindricalPML;
        //! Gets the permeability of the material
        void GetPermeability(const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &ur) const override;
        //! Gets the permittivity of the material
        void GetPermittivity(const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &er) const override;
    };
#else
    using TwistedWgmaPML=TPZCombinedSpacesCylindricalPML<TwistedWgma>;
#endif
};
#endif /* _TWISTED_WGMA_H_ */
