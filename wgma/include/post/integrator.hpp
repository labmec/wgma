#ifndef _INTEGRATOR_HPP_
#define _INTEGRATOR_HPP_

#include <pzcmesh.h>
#include <tpzautopointer.h>

#include "eldata.hpp"

template<class T>
class TPZMaterialDataT;

//! Contains a set of routines for manipulating FEM solutions
namespace wgma::post{
  /** @brief Base class, integrates a given FEM solution.
      Actual contribution at the integration points is defined
      in the derived classes*/
  class Integrator{
  public:
    Integrator() = default;
    virtual ~Integrator() = default;
    Integrator(const Integrator &) = default;
    Integrator(Integrator &&) = default;
    Integrator& operator=(const Integrator &) = default;
    Integrator& operator=(Integrator &&) = default;
    //! Set number of threads to be used in the computations
    void SetNThreads(int nt){m_nthreads = nt;}
    //! Get number of threads
    int NThreads() const {return m_nthreads;}
  protected:
    Integrator(TPZAutoPointer<TPZCompMesh> mesh,
               std::set<int> matids,
               int nThreads);
    
    Integrator(TPZAutoPointer<TPZCompMesh> mesh,
               TPZVec<TPZCompEl*> elvec,
               int nThreads) : m_cmesh(mesh),
                               m_elvec(elvec),
                               m_nthreads(nThreads){}
    
    //! Gets access to the mesh
    TPZAutoPointer<TPZCompMesh> Mesh() {return m_cmesh;}
    /** 
        @brief Computes an integral over given elements.
        The integral to be performed is provided by the virtual method Compute.
        @param[in] elvec elements in which the quantity will be integrated.
        @param[in] data needed data for computing the solution at the integration point
*/
    void Integrate(const TPZVec<TPZCompEl*> &elvec);
    //! Creates instance of element data
    virtual ElData* CreateElData() {return new ElData;}
    //! Initialises element data
    virtual void InitData(TPZCompEl* el, ElData& data) = 0;
    //! Gets element data at an integration point
    virtual void IntPointData(TPZCompEl* el, ElData& data, TPZVec<REAL> &x) = 0;
    //! Allows for inspecting element data after integrating
    virtual void PostProcessData(ElData& data) {}
    //! Override this method with the desired calculation
    virtual void Compute(const ElData &data, REAL weight, int thread) = 0;
    //! Computational mesh
    TPZAutoPointer<TPZCompMesh> m_cmesh{nullptr};
    //! List of elements for which the solution will be computed
    TPZVec<TPZCompEl*> m_elvec;
    //!Number of threads
    int m_nthreads{4};
  };

  //! Integrator for comp meshes containing one approx space
  class SingleSpaceIntegrator : public Integrator{
  public:
    //! Constructors
    using Integrator::Integrator;
    //! Initialises element data
    void InitData(TPZCompEl* el, ElData& data) override;
    //! Gets element data at an integration point
    void IntPointData(TPZCompEl* el, ElData& data, TPZVec<REAL> &x) override;
  };

  //! Integrator for comp meshes containing multiple approx spaces
  class MultiphysicsIntegrator : public Integrator{
  public:
    //! Constructors
    using Integrator::Integrator;
    //! Initialises element data
    void InitData(TPZCompEl* el, ElData& data) override;
    //! Gets element data at an integration point
    void IntPointData(TPZCompEl* el, ElData& data, TPZVec<REAL> &x) override;
  };
};

#endif /* _INTEGRATOR_HPP_ */
