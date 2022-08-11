#ifndef _SCATTERING_HPP_
#define _SCATTERING_HPP_

#include <cmeshtools.hpp>
#include <variant>

#include <pzcmesh.h>
#include <TPZMatrixSolver.h>
#include <TPZLinearAnalysis.h>

using namespace std::string_view_literals;


namespace wgma::scattering{
  /**
     @brief  Class responsible for managing the scattering analysis of planar waveguides
  */
  class Analysis{
  public:
    /**
       @brief Creates the analysis module based on a given computational mesh
       as returned by cmeshtools::CMeshScattering2D
       @param [in] n_threads Number of threads to be used in the analysis
       @param [in] meshvec Vector containing the computational meshes
       @param [in] reorder_eqs whether the equations are reordered for optimising bandwidth
       @param [in] filter_bound whether to impose homogeneous dirichlet BCs by removing the equations
    */
    Analysis(const TPZAutoPointer<TPZCompMesh> cmesh,
                       const int n_threads, const bool reorder_eqs=true,
                       const bool filter_bound=true);
    //! Sets a custom linear solver to be copied to underlying TPZAnalysis(advanced)
    void SetSolver(const TPZMatrixSolver<CSTATE> &solv);
    /**
       @brief Gets a copy of the linear solver for easier configuration (advanced)
       @note A call to Analysis::SetSolver must be made afterwards.
    */
    TPZMatrixSolver<CSTATE> & GetSolver() const;
    /**
       @brief Run the analysis
     */
    void Run();

    
  protected:
    //! H1 mesh
    TPZAutoPointer<TPZCompMesh> m_cmesh{nullptr};
    //! Analysis instance
    TPZAutoPointer<TPZLinearAnalysis> m_an{nullptr};
    //! Whether the matrices have been assembled already
    bool m_assembled{false};
    //! Number of H1 dofs
    int m_n_dofs{-1};
    //! Whether the equations have been filtered
    bool m_filter_bound{false};
  };


  using srcfunc = std::function<void (const TPZVec<REAL> &loc,
                                      CSTATE &val,
                                      TPZVec<CSTATE> &deriv)>;
  /**
     @brief Struct using for using prescribed (analytical) sources in
     scattering analysis.
   */
  struct Source1D{
    std::set<int> id;//< material identifiers for all regions with a prescribed source
    srcfunc func;//function
  };

  /**
     @brief Struct using for using solutions from modal analysis as sources in
     scattering analysis.
   */
  struct SourceWgma{
    std::set<int> id;//< material identifiers for all regions with a prescribed source
    TPZAutoPointer<TPZCompMesh> modal_cmesh;//function
  };

  /**
     @brief Creates the computational mesh used for the scattering analysis of planar waveguides.
     The mesh will be a H1 conforming approximation space and it will approximate
     either Ex or Hx, depending on whether TE/TM modes are approximated.
     The source can be provided by two different means:
     either by 
     - providing a computational mesh and a material identifier where
     the solution will be evaluated
     - providing an analytical source via wgma::scattering::source2D
     @note The computational domain is in the xy-plane, even though physically it
     corresponds to the yz-plane.
     @param[in] gmesh geometrical mesh
     @param[in] mode whether to solve for TE or TM modes
     @param[in] pOrder polynomial order
     @param[in] data information regarding domain's regions
     @param[in] source contains the function that will excite the waveguide
     @param[in] lambda operational wavelength
     @param[in] scale geometric scaling (characteristic length) for better floating point precision
  */
  TPZAutoPointer<TPZCompMesh>
  CMeshScattering2D(TPZAutoPointer<TPZGeoMesh> gmesh,
                    const planarwg::mode mode, int pOrder,
                    cmeshtools::PhysicalData &data,
                    std::variant<
                    wgma::scattering::Source1D,
                    wgma::scattering::SourceWgma> source,
                    const STATE lambda, const REAL scale);

  void
  SetPropagationConstant(TPZAutoPointer<TPZCompMesh> cmesh,
                         const CSTATE beta);
};

#endif /* _SCATTERING_HPP_ */
