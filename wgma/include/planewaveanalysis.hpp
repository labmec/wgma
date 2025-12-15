#ifndef _PWANALYSIS_HPP_
#define _PWANALYSIS_HPP_

#include <cmeshtools.hpp>

#include <pzcmesh.h>
#include <TPZMatrixSolver.h>
#include <TPZLinearAnalysis.h>

namespace wgma::planewaveanalysis{
  //temporary variable (set to false in .cpp)
  extern bool using_tbb_mat;

  class Analysis : public TPZLinearAnalysis {
  public:
    /**
       @brief Creates the analysis module based on a given computational mesh
       as returned by cmeshtools::CMeshScattering2D
       @param [in] n_threads Number of threads to be used in the analysis
       @param [in] meshvec Vector containing the computational meshes
       @param [in] reorder_eqs whether the equations are reordered for optimising bandwidth
       @param [in] filter_bound whether to impose homogeneous dirichlet BCs by removing the equations
    */
    Analysis(const TPZVec<TPZAutoPointer<TPZCompMesh>> &meshvec,
             const int n_threads, const bool reorder_eqs=true,
             const bool filter_bound=true);

    void
    CountActiveEqs(int64_t &neq,int64_t &nH1Equations, int64_t &nHCurlEquations);
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
    void Assemble() override;
    /**
       @brief Assembles the rhs of the algebraic system
       @param[in] identifiers of the materials to be assembled (source materials)
     */
    void AssembleRhs(std::set<int> matids);
    /**
       @brief Solves the algebraic system
    */
    void Solve() override;
    /**
       @brief Assembles and solves the algebraic system
     */
    void Run();

    [[nodiscard]] const TPZAutoPointer<TPZCompMesh> GetMesh(){return m_cmesh_mf;}

    [[nodiscard]] const TPZAutoPointer<TPZCompMesh> GetHCurlMesh(){return m_cmesh_hcurl;}
    [[nodiscard]] const TPZAutoPointer<TPZCompMesh> GetH1Mesh(){return m_cmesh_h1;}

    [[nodiscard]] TPZVec<CSTATE> GetEigenvalues() {return m_eigenvalues;}

    
    
    void LoadSolution() override;
    
    using TPZAnalysis::StructMatrix;

    using TPZAnalysis::LoadSolution;

    using TPZAnalysis::Solution;

    using TPZLinearAnalysis::Rhs;

  protected:
    //! Combined computational mesh (hcurl and h1)
    TPZAutoPointer<TPZCompMesh> m_cmesh_mf{nullptr};
    //! H1 mesh
    TPZAutoPointer<TPZCompMesh> m_cmesh_h1{nullptr};
    //! Hcurl mesh
    TPZAutoPointer<TPZCompMesh> m_cmesh_hcurl{nullptr};
    //! Analytical eigenvalues
    TPZVec<CSTATE> m_eigenvalues;
    //! Total number of dofs
    int64_t m_n_dofs_mf{-1};
    //! Number of H1 dofs
    int64_t m_n_dofs_h1{-1};
    //! Number of HCurl dofs
    int64_t m_n_dofs_hcurl{-1};
    //! Whether the matrices have been assembled already
    bool m_assembled{false};
    //! Whether the equations have been filtered
    bool m_filter_bound{false};
    //! Indices of bound connects (do not remember why it is needed)
    std::set<int64_t> m_bound_cons;
  };
  

  /**
     @brief Split materials between volumetric/pml materials and find BC neighbours
     @param [in] gmesh geometrical mesh
     @param [in/out] data information regarding domain's reginos
     @param [out] volmats mat ids of volumetric materials
     @param [out] pmlmats mat ids of pml materials
   */
  void SetupModalAnalysisMaterials(TPZAutoPointer<TPZGeoMesh> &gmesh, cmeshtools::PhysicalData& data,
                                   std::set<int>&volmats, std::set<int>&pmlmats);

  /**
     @brief Creates the computational meshes used for projecting plane wave
     solutions for the periodic homogeneous waveguide EVP in two dimensions.
     Three meshes will be created: one for the H1 approximation space, one for the
     HCurl approximation space and one multiphysics mesh combining both spaces.
     @param[in] gmesh geometrical mesh
     @param[in] pOrder polynomial order
     @param[in] data information regarding domain's regions
     @param[in] periodic_els map
     @param[in] lambda operational wavelength
     @param[in] lx domain length in x direction
     @param[in] ly domain length in y direction
     @param[in] max_k maximum integer value for computing solutions
     @param[in] scale geometric scaling (characteristic length) for better floating point precision
  */
  TPZVec<TPZAutoPointer<TPZCompMesh>>
  CMeshPlaneWave2D(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
                   cmeshtools::PhysicalData &data,
                   const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &el_map,
                   const STATE lambda,
                   const STATE lx,
                   const STATE ly,
                   const int max_k,
                   const REAL &scale,
                   bool verbose);

  //!aux function used in CMeshPlaneWave2D
  TPZAutoPointer<TPZCompMesh>
  CreateAtomicWgma2D(TPZAutoPointer<TPZGeoMesh> gmesh,
                     bool isH1,
                     int p,
                     const std::set<int> &volmats,
                     const std::set<int> &pmlmats,
                                                 
                     const std::vector<wgma::bc::data> &bcmats,
                     const std::vector<std::pair<int,int>> &probevec);
};

#endif /* _PWANALYSIS_HPP_ */
