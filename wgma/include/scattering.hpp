#ifndef _SCATTERING_HPP_
#define _SCATTERING_HPP_

#include <pzcmesh.h>
#include <TPZMatrixSolver.h>
#include <TPZLinearAnalysis.h>

namespace wgma{
  /**
     @brief  Class responsible for managing the scattering analysis of planar waveguides
  */
  class ScatteringAnalysis{
  public:
    /**
       @brief Creates the analysis module based on a given computational mesh
       as returned by cmeshtools::CMeshScattering2D
       @param [in] n_threads Number of threads to be used in the analysis
       @param [in] meshvec Vector containing the computational meshes
       @param [in] reorder_eqs whether the equations are reordered for optimising bandwidth
       @param [in] filter_bound whether to impose homogeneous dirichlet BCs by removing the equations
    */
    ScatteringAnalysis(const TPZAutoPointer<TPZCompMesh> cmesh,
                       const int n_threads, const bool reorder_eqs=true,
                       const bool filter_bound=true);
    //! Sets a custom linear solver to be copied to underlying TPZAnalysis(advanced)
    void SetSolver(const TPZMatrixSolver<CSTATE> &solv);
    /**
       @brief Gets a copy of the linear solver for easier configuration (advanced)
       @note A call to ScatteringAnalysis::SetSolver must be made afterwards.
    */
    TPZMatrixSolver<CSTATE> & GetSolver() const;
    /**
       @brief Run the analysis
     */
    void Run();
    /**
       @brief Post-process the solution in .vtk format
       @param[in] filename without extension in which to export the solution
       @param[in] vtk_res resolution of the .vtk file (number of element subdivisions)
     */
    void PostProcess(std::string filename, const int vtk_res = 0);

    
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
};

#endif /* _SCATTERING_HPP_ */
