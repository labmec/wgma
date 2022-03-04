#ifndef _WGANALYSIS_HPP_
#define _WGANALYSIS_HPP_

#include <pzcmesh.h>
#include <TPZEigenSolver.h>
#include <TPZEigenAnalysis.h>

namespace wgma{
  /**
     @brief  Class responsible for managing the modal analysis of waveguides
     @note The eigensolver still has to be configured.
  */
  class WGAnalysis{
  public:
    /**
       @brief Creates the analysis module based on a given set
       of computational meshes as returned by cmeshtools::CMeshWgma2D
       @param [in] n_threads Number of threads to be used in the analysis
       @param [in] meshvec Vector containing the computational meshes
       @param [in] reorder_eqs whether the equations are reordered for optimising bandwidth
       @param [in] filter_bound whether to impose homogeneous dirichlet BCs by removing the equations
    */
    WGAnalysis(const TPZVec<TPZAutoPointer<TPZCompMesh>> &meshvec,
               const int n_threads, const bool reorder_eqs=true,
               const bool filter_bound=true);
    //! Sets a custom eigensolver to be copied to underlying TPZAnalysis
    void SetSolver(const TPZEigenSolver<CSTATE> &solv);
    /**
       @brief Gets a copy of the eigensolver for easier configuration
       @note A call to WGAnalysis::SetSolver must be made afterwards.
    */
    TPZEigenSolver<CSTATE> & GetSolver() const;
    /**
       @brief Run the analysis
       @param [in] compute_eigenvectors whether to compute eigenvectors (or just eigenvalues)
     */
    void Run(bool compute_eigenvectors);
    /*
      @brief Gets calculated eigenvalues
     */
    TPZVec<CSTATE> GetEigenvalues() const;

    /**
       @brief Post-process the solution in .vtk format
       @param[in] filename without extension in which to export the solution
       @param[in] vtk_res resolution of the .vtk file (number of element subdivisions)
       @param[in] print_real_part whether to print the real part (or magnitude) of the field
     */
    void PostProcess(std::string filename, const int vtk_res = 0, const bool print_real_part=true);

    /**
       @brief Export eigenvalues to in a csv format and append it to a file.
       The following values are exported:
       neq , nel , h , p , lambda , nev, real(w1), imag(w1) , ... , real(wnev) , imag(wnev)
       where:
       - neq    = number of equations
       - nev    = number of elements
       - h      = radius of the biggest element in the mesh
       - p      = default polynomial order of the hcurl mesh
       - lambda = operational wavelength
       - nev    = number of eigenvalues
       @param [in] filename name of the file 
       @param [in] lambda operational wavelength
     */
    void WriteToCsv(std::string filename, STATE lambda);
      
  protected:
    //! Combined computational mesh (hcurl and h1)
    TPZAutoPointer<TPZCompMesh> m_cmesh_mf{nullptr};
    //! H1 mesh
    TPZAutoPointer<TPZCompMesh> m_cmesh_h1{nullptr};
    //! Hcurl mesh
    TPZAutoPointer<TPZCompMesh> m_cmesh_hcurl{nullptr};
    //! Analysis instance
    TPZAutoPointer<TPZEigenAnalysis> m_an{nullptr};
    //! Calculated eigenvalues
    TPZVec<CSTATE> m_evalues;
    //! Calculated eigenvectors
    TPZFMatrix<CSTATE> m_evectors;
    //! Whether the matrices have been assembled already
    bool m_assembled{false};
    //! Total number of dofs
    int m_n_dofs_mf{-1};
    //! Number of H1 dofs
    int m_n_dofs_h1{-1};
    //! Number of HCurl dofs
    int m_n_dofs_hcurl{-1};
    //! Whether the equations have been filtered
    bool m_filter_bound{false};
  };
};

#endif /* _WGANALYSIS_HPP_ */
