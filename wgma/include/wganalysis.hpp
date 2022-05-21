#ifndef _WGANALYSIS_HPP_
#define _WGANALYSIS_HPP_

#include <cmeshtools.hpp>

#include <pzcmesh.h>
#include <TPZEigenSolver.h>
#include <TPZEigenAnalysis.h>

namespace wgma::wganalysis{
  /**
     @brief  Abstract base class responsible for managing the modal analysis of waveguides
  */
  class Wgma{
  public:
    //! Default constructor
    Wgma() = default;
    //! Copy destructor
    Wgma(const Wgma &) = default;
    //! Move constructor
    Wgma(Wgma &&) = default;
    //! Copy assignment operator
    Wgma& operator=(const Wgma &) = default;
    //! Move assignment operator
    Wgma& operator=(Wgma &&) = default;
    //! Default destructor
    virtual ~Wgma();
    //default dtor in the cpp file to ensure correct creation of vtable (https://stackoverflow.com/a/57504289/1870801)
    
    //! Sets a custom eigensolver to be copied to underlying TPZWgma2D
    void SetSolver(const TPZEigenSolver<CSTATE> &solv);
    /**
       @brief Gets a copy of the eigensolver for easier configuration
       @note A call to Wgma2D::SetSolver must be made afterwards.
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
    virtual void WriteToCsv(std::string filename, STATE lambda) = 0;
  protected:
    //! Perform necessary adjustments on the eigensolver
    virtual void AdjustSolver(TPZEigenSolver<CSTATE> *solv) {}
    //! Wgma2D instance
    TPZAutoPointer<TPZEigenAnalysis> m_an{nullptr};
    //! Calculated eigenvalues
    TPZVec<CSTATE> m_evalues;
    //! Calculated eigenvectors
    TPZFMatrix<CSTATE> m_evectors;
    //! Whether the matrices have been assembled already
    bool m_assembled{false};
    //! Whether the equations have been filtered
    bool m_filter_bound{false};
    //! Indices of bound connects (do not count as additional eqs)
    std::set<int64_t> m_bound_cons;
  };

  /**
     @brief  Class responsible for managing the modal analysis of waveguides with a 2D cross section
     @note The eigensolver still has to be configured.
  */
  class Wgma2D : public Wgma{
  public:
    /**
       @brief Creates the analysis module based on a given set
       of computational meshes as returned by cmeshtools::CMeshWgma2D
       @param [in] n_threads Number of threads to be used in the analysis
       @param [in] meshvec Vector containing the computational meshes
       @param [in] reorder_eqs whether the equations are reordered for optimising bandwidth
       @param [in] filter_bound whether to impose homogeneous dirichlet BCs by removing the equations
    */
    Wgma2D(const TPZVec<TPZAutoPointer<TPZCompMesh>> &meshvec,
               const int n_threads, const bool reorder_eqs=true,
               const bool filter_bound=true);
    
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
    void WriteToCsv(std::string filename, STATE lambda) override;

    /** @brief Counts active equations per approximation space for the 2D waveguide modal analysis.*/
    void CountActiveEqs(int &neq, int&nh1, int &nhcurl);
  private:
    void AdjustSolver(TPZEigenSolver<CSTATE> *solv) override;
    //! Combined computational mesh (hcurl and h1)
    TPZAutoPointer<TPZCompMesh> m_cmesh_mf{nullptr};
    //! H1 mesh
    TPZAutoPointer<TPZCompMesh> m_cmesh_h1{nullptr};
    //! Hcurl mesh
    TPZAutoPointer<TPZCompMesh> m_cmesh_hcurl{nullptr};
    //! Total number of dofs
    int m_n_dofs_mf{-1};
    //! Number of H1 dofs
    int m_n_dofs_h1{-1};
    //! Number of HCurl dofs
    int m_n_dofs_hcurl{-1};
  };

  /**
     @brief Performs modal analysis of periodic planar waveguides.
     Since it is a non-linear eigenvalue in beta it is supposed
     to be solved iteratively.
     See pcwg example for usage.
   */
  class WgmaPeriodic2D : public Wgma{
  public:
    /**
       @brief Creates the analysis module based on a given set
       of computational meshes as returned by cmeshtools::CMeshWgma2D
       @param [in] n_threads Number of threads to be used in the analysis
       @param [in] meshvec Vector containing the computational meshes
       @param [in] reorder_eqs whether the equations are reordered for optimising bandwidth
       @param [in] filter_bound whether to impose homogeneous dirichlet BCs by removing the equations
    */
    WgmaPeriodic2D(TPZAutoPointer<TPZCompMesh> cmesh,
                   const int n_threads, const bool reorder_eqs=true,
                   const bool filter_bound=true);

    

    //! Sets propagation constant to be used in the material
    void SetBeta(const CSTATE beta);
    
    /**
       @brief Post-process the solution in .vtk format
       @param[in] filename without extension in which to export the solution
       @param[in] vtk_res resolution of the .vtk file (number of element subdivisions)
     */
    void PostProcess(std::string filename, const int vtk_res = 0);

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
    void WriteToCsv(std::string filename, STATE lambda) override;

    /** @brief Counts active equations per approximation space for the 2D waveguide modal analysis.*/
    void CountActiveEqs(int &neq);
  private:
    void AdjustSolver(TPZEigenSolver<CSTATE> *solv) override;
    //! Computational mesh
    TPZAutoPointer<TPZCompMesh> m_cmesh{nullptr};
    //! Total number of dofs
    int m_n_dof{-1};
    //! Current value of propagation constant
    CSTATE m_beta{0};
  };

  /**
     @brief Creates the computational meshes used for approximating the waveguide EVP in two dimensions.
     Three meshes will be created: one for the H1 approximation space, one for the
     HCurl approximation space and one multiphysics mesh combining both spaces.
     @param[in] gmesh geometrical mesh
     @param[in] pOrder polynomial order
     @param[in] data information regarding domain's regions
     @param[in] lambda operational wavelength
     @param[in] scale geometric scaling (characteristic length) for better floating point precision
  */
  TPZVec<TPZAutoPointer<TPZCompMesh>>
  CMeshWgma2D(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
              cmeshtools::PhysicalData &data,
              const STATE lambda, const REAL &scale);

  //creates computational mesh for modal analysis of periodic planar waveguides
  TPZAutoPointer<TPZCompMesh>
  CMeshWgmaPeriodic(TPZAutoPointer<TPZGeoMesh> gmesh,
                    const wgma::planarwg::mode mode, int pOrder,
                    wgma::cmeshtools::PhysicalData &data,
                    std::map<int64_t, int64_t> periodic_els, const STATE lambda,
                    const REAL scale);
};

#endif /* _WGANALYSIS_HPP_ */
