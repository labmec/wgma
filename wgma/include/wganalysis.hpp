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
  class Wgma : private TPZEigenAnalysis{
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
    
    //! Sets a custom eigensolver to be copied to this instance
    using TPZEigenAnalysis::SetSolver;
    /**
       @brief Gets a copy of the eigensolver for easier configuration
       @note A call to Wgma2D::SetSolver must be made afterwards.
    */
    TPZEigenSolver<CSTATE> & GetSolver(){
      return TPZEigenAnalysis::EigenSolver<CSTATE>();
    }
    /**
       @brief Assemble both FEM matrices
    */
    using TPZEigenAnalysis::Assemble;
    /**
       @brief Assemble one of the FEM matrices (usefull for nlin problems) 
     */
    void Assemble(TPZEigenAnalysis::Mat mat);
    /**
       @brief Run the analysis
       @param [in] compute_eigenvectors whether to compute eigenvectors (or just eigenvalues)
     */
    void Run(bool compute_eigenvectors);
    /**
       @brief Solve the eigen problem
       @param [in] compute_eigenvectors whether to compute eigenvectors (or just eigenvalues)
     */
    void Solve(bool compute_eigenvectors);

    /*
      @brief Loads the isol-th solution (eigenvector) in the computational mesh.
      This method is specially useful if the solution from the modal analysis
      will be used in another FEM scheme (as a source, for instance).
      @note Derived classes will implement LoadSolutionInternal.
     */
    void LoadSolution(const int isol);

    //! Loads all computed eigenvectors into the mesh
    void LoadAllSolutions();
    /*
      @brief Gets calculated eigenvalues
     */
    using TPZEigenAnalysis::GetEigenvalues;
    //! Sets eigenvalues
    using TPZEigenAnalysis::SetEigenvalues;
    /*
      @brief Gets calculated eigenvectors
     */
    using TPZEigenAnalysis::GetEigenvectors;
    //! Sets eigenvectors
    using TPZEigenAnalysis::SetEigenvectors;

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

    [[nodiscard]] virtual TPZAutoPointer<TPZCompMesh> GetMesh() = 0;
  protected:

    virtual void LoadSolutionInternal(const int isol, const int nsol) = 0;
    using TPZEigenAnalysis::LoadSolution;
    using TPZAnalysis::SetCompMeshInit;
    using TPZAnalysis::SetStructuralMatrix;
    using TPZAnalysis::StructMatrix;
    //! Perform necessary adjustments on the eigensolver
    virtual void AdjustSolver(TPZEigenSolver<CSTATE> *solv) {}
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

    [[nodiscard]] TPZAutoPointer<TPZCompMesh> GetMesh() override{return m_cmesh_mf;}
  private:
    void LoadSolutionInternal(const int isol, const int ncols) override;
    
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
     @brief Base class for performing modal analysis of planar waveguides.
     See slab_disc example for usage.
   */
  class WgmaPlanar : public Wgma{
  public:
    /**
       @brief Creates the analysis module based on a given set
       of computational meshes as returned by cmeshtools::CMeshWgma2D
       @param [in] n_threads Number of threads to be used in the analysis
       @param [in] meshvec Vector containing the computational meshes
       @param [in] reorder_eqs whether the equations are reordered for optimising bandwidth
       @param [in] filter_bound whether to impose homogeneous dirichlet BCs by removing the equations
    */
    WgmaPlanar(TPZAutoPointer<TPZCompMesh> cmesh,
               const int n_threads, const bool reorder_eqs=true,
               const bool filter_bound=true);
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

    [[nodiscard]] TPZAutoPointer<TPZCompMesh> GetMesh() override{return m_cmesh;}
  protected:
    void LoadSolutionInternal(const int isol, const int nsol) override;
    //! Computational mesh
    TPZAutoPointer<TPZCompMesh> m_cmesh{nullptr};
    //! Total number of dofs
    int m_n_dof{-1};
  };

    /**
     @brief Performs modal analysis of periodic planar waveguides.
     Since it is a non-linear eigenvalue in beta it is supposed
     to be solved iteratively.
     See pcwg example for usage.
     @note It will always search just for one eigenvalue.
   */
  class WgmaPeriodic2D : public WgmaPlanar{
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
                   const bool filter_bound=true) :
      WgmaPlanar(cmesh,n_threads,reorder_eqs,filter_bound) {}

    

    //! Sets propagation constant to be used in the material
    void SetBeta(const CSTATE beta);
  private:
    void AdjustSolver(TPZEigenSolver<CSTATE> *solv) override;
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

  /**
     @brief Creates the computational meshes used for approximating the waveguide EVP in one dimension, for planar waveguides.
     @param[in] gmesh geometrical mesh
     @param[in] mode whether to solve for TE or TM modes
     @param[in] pOrder polynomial order
     @param[in] data information regarding domain's regions
     @param[in] lambda operational wavelength
     @param[in] scale geometric scaling (characteristic length) for better floating point precision
  */
  TPZAutoPointer<TPZCompMesh>
  CMeshWgma1D(TPZAutoPointer<TPZGeoMesh> gmesh,
              wgma::planarwg::mode mode, int pOrder,
              cmeshtools::PhysicalData &data,
              const STATE lambda, const REAL scale);

  //creates computational mesh for modal analysis of periodic planar waveguides
  TPZAutoPointer<TPZCompMesh>
  CMeshWgmaPeriodic(TPZAutoPointer<TPZGeoMesh> gmesh,
                    const wgma::planarwg::mode mode, int pOrder,
                    wgma::cmeshtools::PhysicalData &data,
                    std::map<int64_t, int64_t> periodic_els, const STATE lambda,
                    const REAL scale);
};

#endif /* _WGANALYSIS_HPP_ */
