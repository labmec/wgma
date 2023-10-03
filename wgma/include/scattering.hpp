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
  class Analysis : private TPZLinearAnalysis {
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
             const bool filter_bound=true,
             const bool is_sym=true);
    //! Sets a custom linear solver to be copied to underlying TPZAnalysis(advanced)
    using TPZLinearAnalysis::SetSolver;
    /**
       @brief Gets a copy of the linear solver for easier configuration (advanced)
       @note A call to Analysis::SetSolver must be made afterwards.
    */
    TPZMatrixSolver<CSTATE> & GetSolver(){
      return TPZLinearAnalysis::MatrixSolver<CSTATE>();
    }
    /**
       @brief Assembles the algebraic system
     */
    void Assemble();
    /**
       @brief Assembles the rhs of the algebraic system
       @param[in] identifiers of the materials to be assembled (source materials)
     */
    void AssembleRhs(std::set<int> matids);
    /**
       @brief Solves the algebraic system
    */
    void Solve();
    /**
       @brief Assembles and solves the algebraic system
     */
    void Run();

   [[nodiscard]] TPZAutoPointer<TPZCompMesh> GetMesh() {return m_cmesh;}

    void SetSymMatrix(bool sym){m_sym = sym;}

    using TPZAnalysis::StructMatrix;

    using TPZAnalysis::LoadSolution;

    using TPZAnalysis::Solution;

    using TPZLinearAnalysis::Rhs;

    TPZAutoPointer<TPZMatrixSolver<CSTATE>> BuildBlockPrecond(const TPZVec<int64_t> &eqgraph,
                                                              const TPZVec<int64_t> &graphindex,
                                                              const bool overlap);
    using TPZAnalysis::BuildPreconditioner;
  protected:
    //! H1 mesh
    TPZAutoPointer<TPZCompMesh> m_cmesh{nullptr};
    //! Whether the matrices have been assembled already
    bool m_assembled{false};
    //! Number of H1 dofs
    int m_n_dofs{-1};
    //! Whether the equations have been filtered
    bool m_filter_bound{false};
    //! Whether to use a symmetric matrix storage format
    bool m_sym{true};
  };


  using srcfunc1d = std::function<void (const TPZVec<REAL> &loc,
                                      CSTATE &val,
                                      TPZVec<CSTATE> &deriv)>;
  /**
     @brief Struct using for using prescribed (analytical) sources in
     scattering analysis in two dimensions.
   */
  struct Source1D{
    std::set<int> id;//< material identifiers for all regions with a prescribed source
    srcfunc1d func;//function
  };

  using srcfunc2d = std::function<void (const TPZVec<REAL> &loc,
                                        TPZVec<CSTATE> &val)>;
  /**
     @brief Struct using for using prescribed (analytical) sources in
     scattering analysis in three dimensions.
   */
  struct Source2D{
    std::set<int> id;//< material identifiers for all regions with a prescribed source
    srcfunc2d func;//function
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
     For loading a source, see scattering::LoadSource and scattering::SetPropagationConstant. 
     @note The computational domain is in the xy-plane, even though physically it
     corresponds to the yz-plane.
     @param[in] gmesh geometrical mesh
     @param[in] mode whether to solve for TE or TM modes
     @param[in] pOrder polynomial order
     @param[in] data information regarding domain's regions
     @param[in] source_ids contains the ids of the excitation source regions
     @param[in] lambda operational wavelength
     @param[in] scale geometric scaling (characteristic length) for better floating point precision
  */
  TPZAutoPointer<TPZCompMesh>
  CMeshScattering2D(TPZAutoPointer<TPZGeoMesh> gmesh,
                    const planarwg::mode mode, int pOrder,
                    cmeshtools::PhysicalData &data,
                    const std::set<int> source_ids,
                    const STATE lambda, const REAL scale);

  /**
     @brief Creates the computational mesh used for the scattering analysis in 3D.
Usually the source will be the result of a previously computed modal analysis.
     The mesh will be a Hcurl conforming approximation space and it will approximate
     the electric field.
     For loading a source, see scattering::LoadSource and scattering::SetPropagationConstant. 
     @param[in] gmesh geometrical mesh
     @param[in] pOrder polynomial order
     @param[in] data information regarding domain's regions
     @param[in] source_ids contains the ids of the excitation source regions
     @param[in] lambda operational wavelength
     @param[in] scale geometric scaling (characteristic length) for better floating point precision
  */
  TPZAutoPointer<TPZCompMesh>
  CMeshScattering3D(TPZAutoPointer<TPZGeoMesh> gmesh,
                    int pOrder,
                    cmeshtools::PhysicalData &data,
                    const std::set<int> source_ids,
                    const STATE lambda, const REAL scale);

  /**
     @brief Set the propagation constant value for the source of the scattering analysis.
     @tparam[in] type of the material corresponding to the source
     @param[in] cmesh computational mesh of the scattering problem
     @param[in] beta propagation constant
  */
  void
  SetPropagationConstant(TPZAutoPointer<TPZCompMesh> cmesh,
                         const CSTATE beta);

  /**
     @brief Loads a source for the scattering analysis in two dimensions.
     The source can be either an analytical source or
     a source from the modal analysis.

     The source can be provided by two different means:
     either by 
     - providing a computational mesh and a material identifier where
     the solution will be evaluated
     - providing an analytical source via wgma::scattering::source1D

     The value of the propagation constant beta is set by SetPropagationConstant.

     @param[in] scatt_cmesh computational mesh of the scattering problem
     @param[in] source source for the scattering problem
     @param[in] coeff coefficient for modal analysis source (amplitude)
  */
  void
  LoadSource1D(TPZAutoPointer<TPZCompMesh> cmesh,
             std::variant<
             wgma::scattering::Source1D,
             wgma::scattering::SourceWgma> source,
             CSTATE coeff=1.);

  /**
     @brief Loads a source for the scattering analysis in three dimensions.
     The source can be either an analytical source or
     a source from the modal analysis.

     The source can be provided by two different means:
     either by 
     - providing a computational mesh and a material identifier where
     the solution will be evaluated
     - providing an analytical source via wgma::scattering::source2D

     The value of the propagation constant beta is set by SetPropagationConstant.

     @param[in] scatt_cmesh computational mesh of the scattering problem
     @param[in] source source for the scattering problem
     @param[in] coeff coefficient for modal analysis source (amplitude)
  */
  void
  LoadSource2D(TPZAutoPointer<TPZCompMesh> cmesh,
             std::variant<
             wgma::scattering::Source2D,
             wgma::scattering::SourceWgma> source,
             CSTATE src=1.);
};

#endif /* _SCATTERING_HPP_ */
