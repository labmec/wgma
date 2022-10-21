/**
   @file slepcepshandler.hpp
   This contains a handler for the SLEPc EPS Module
*/
#ifndef WGMASLEPCHANDLER_H
#define WGMASLEPCHANDLER_H

#include <TPZEigenSolver.h>

namespace wgma::slepc{
  /*
    The following are simple wrappers for avoiding exposing SLEPc types
  */
  /**
     @brief Eigenvalue problem type (generalised, hermitian, etc)
  */
  enum class EPSProblemType{
    EPS_HEP=1,///< Hermitian eigenvalue problem
    EPS_GHEP,///< Generalized Hermitian eigenvalue problem
    EPS_NHEP,///< Non-Hermitian eigenvalue problem
    EPS_GNHEP,///< Generalized non-Hermitian eigenvalue problem
    EPS_PGNHEP,///< Generalized non-Hermitian eigenvalue problem with positive semi-definite B
    EPS_GHIEP,///< Generalized Hermitian-indefinite eigenvalue problem
    EPS_NOTSET///<Default value
  };
  
  /**
     @brief How to evaluate convergence
  */
  enum class EPSConv{
    EPS_CONV_ABS,///< Sets the absolute convergence test 
    EPS_CONV_REL,///< Sets the convergence test relative to the eigenvalue 
    EPS_CONV_NORM///< Sets the convergence test relative to the matrix norms 
  };

  /**
     @brief Type of eigensolver
  */
  enum class EPSType{
    POWER, ///< Power
    SUBSPACE,///< Subspace iteration
    ARNOLDI,///< Arnoldi
    LANCZOS,///< Lanczos
    KRYLOVSCHUR,///< Krylov-Schur
    GD,///< Generalised Davidson
    JD///< Jacobi-Davidson
  };
  /**
     @brief Which eigenvalues are to be sought
  */
  enum class EPSWhich{
    EPS_LARGEST_MAGNITUDE,///< Largest eigenvalues in magnitude
    EPS_SMALLEST_MAGNITUDE,///< Smallest eigenvalues in magnitude
    EPS_LARGEST_REAL,///< Largest real parts
    EPS_SMALLEST_REAL,///< Smallest real parts
    EPS_LARGEST_IMAGINARY,///< Largest imaginary parts
    EPS_SMALLEST_IMAGINARY,///< Smallest imaginary parts
    EPS_TARGET_MAGNITUDE,///< Eigenvalues closest to the target (in magnitude)
    EPS_TARGET_REAL,///< Eigenvalues with real part closest to target
    EPS_TARGET_IMAGINARY///< Eigenvalues with imaginary part closest to target
  };

  


  /**
     @brief Linear solver to be used.
     In the shift-and-invert method, a matrix is inverted.
     This class selects the solver to be used.
     @note Suggested: KSSPPREONLY
  */
  enum class KSPSolver{
    RICHARDSON,///< Check PETSC manual
    CHEBYSHEV,///< Check PETSC manual
    CG,///< Check PETSC manual
    GROPPCG,///< Check PETSC manual
    PIPECG,///< Check PETSC manual
    PIPECGRR,///< Check PETSC manual
    PIPELCG,///< Check PETSC manual
    PIPEPRCG,///< Check PETSC manual
    PIPECG2,///< Check PETSC manual
    CGNE,///< Check PETSC manual
    NASH,///< Check PETSC manual
    STCG,///< Check PETSC manual
    GLTR,///< Check PETSC manual
    FCG,///< Check PETSC manual
    PIPEFCG,///< Check PETSC manual
    GMRES,///< Check PETSC manual
    PIPEFGMRES,///< Check PETSC manual
    FGMRES,///< Check PETSC manual
    LGMRES,///< Check PETSC manual
    DGMRES,///< Check PETSC manual
    PGMRES,///< Check PETSC manual
    TCQMR,///< Check PETSC manual
    BCGS,///< Check PETSC manual
    IBCGS,///< Check PETSC manual
    FBCGS,///< Check PETSC manual
    FBCGSR,///< Check PETSC manual
    BCGSL,///< Check PETSC manual
    PIPEBCGS,///< Check PETSC manual
    CGS,///< Check PETSC manual
    TFQMR,///< Check PETSC manual
    CR,///< Check PETSC manual
    PIPECR,///< Check PETSC manual
    LSQR,///< Check PETSC manual
    PREONLY,///< Check PETSC manual
    QCG,///< Check PETSC manual
    BICG,///< Check PETSC manual
    MINRES,///< Check PETSC manual
    SYMMLQ,///< Check PETSC manual
    LCD,///< Check PETSC manual
    PYTHON,///< Check PETSC manual
    GCR,///< Check PETSC manual
    PIPEGCR,///< Check PETSC manual
    TSIRM,///< Check PETSC manual
    CGLS,///< Check PETSC manual
    FETIDP,///< Check PETSC manual
    HPDDM///< Check PETSC manual
  };
  /**
     @brief Preconditioner to be used.
     @note Its choice is thightly related to the KSPSolver.*/
  enum class Precond{
    NONE,///< Uses just KSP Solver
    JACOBI,///< Check PETSc manual
    SOR,///< Check PETSc manual
    LU,///< Check PETSc manual
    SHELL,///< Check PETSc manual
    BJACOBI,///< Check PETSc manual
    MG,///< Check PETSc manual
    EISENSTAT,///< Check PETSc manual
    ILU,///< Check PETSc manual
    ICC,///< Check PETSc manual
    ASM,///< Check PETSc manual
    GASM,///< Check PETSc manual
    KSP,///< Check PETSc manual
    REDUNDANT,///< Check PETSc manual
    CHOLESKY,///< Check PETSc manual
    PBJACOBI,///< Check PETSc manual
    VPBJACOBI,///< Check PETSc manual
    SVD,///< Check PETSc manual
    BDDC,///< Check PETSc manual
    KACZMARZ,///< Check PETSc manual
    TELESCOPE,///< Check PETSc manual
    PATCH,///< Check PETSc manual
    LMVM,///< Check PETSc manual
    HMG,///< Check PETSc manual
    DEFLATION///< Check PETSc manual
  };

  
  /**
     @brief Handler for EPS Module of SLEPc.
     This class is mainly designed thinking about solving a generalised
     eigenvalue problem using the shift and invert spectral transform.
     Consult SLEPc manual for advanced usage. 
     WARNING: Currently it only supports NeoPZ TPZSpStructMatrix<T> type,
     where std::is_same_v<PetscScalar,T> == true.
     @note All configurations must be made before associating the EPSHandler instance
     with the TPZEigenAnalysis instance.
  */
  template<class TVar>
  class EPSHandler : public TPZEigenSolver<TVar> {
  public:
    //! Initializes SLEPc
    static void InitSLEPc();
    //! Should be called *once* before ending the execution of the program
    static void FinalizeSLEPc();
    //! If TVar!=PetscScalar, it will not compile
    EPSHandler();
      
    EPSHandler(const EPSHandler &) = default;
    EPSHandler(EPSHandler &&) = default;
    EPSHandler& operator=(const EPSHandler &) = default;
    EPSHandler& operator=(EPSHandler &&) = default;
    ~EPSHandler() = default;
    //! This function does NOT clone the instance. It merely passes a pointer to it.
    EPSHandler * Clone() const override;
    
    int SolveEigenProblem(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors) override;
    int SolveEigenProblem(TPZVec<CTVar> &w) override;

    int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w,
                                     TPZFMatrix<CTVar> &eigenVectors) override;

    int SolveGeneralisedEigenProblem(TPZVec<CTVar> &w) override;

    void SetNEigenpairs(int n) override;

    //! Details the problem type (hermitian or not, generalised or not, etc)
    void SetProblemType(const EPSProblemType epsProblem);
    //! Get details about the problem type
    EPSProblemType GetProblemType() const;
    //! It is recommend to call SetProblemType instead
    void SetAsGeneralised(bool isGeneralised) override;
    //! Sets the portion of the spectrum in which evs are to be sought
    void SetWhichEigenpairs(const EPSWhich eps_which);
    //! Gets the portion of the spectrum in which evs are to be sought
    EPSWhich GetWhichEigenpairs() const;
    /**
       @brief Set options related to Krylov-schur algorithm
       @param lock Whether the locking variant is used
       @param restart Percentage of evs kept after restart. negative value for PETSC_DECIDE
       @note ignored if other algorithms are in use.
     */
    void SetKrylovOptions(const bool lock, const RTVar restart);
    /**
       @brief Set options related to Krylov-schur algorithm
       @param lock Whether the locking variant is used
       @param restart Percentage of evs kept after restart
       @note ignored if other algorithms are in use.
     */
    void GetKrylovOptions(bool &lock , RTVar &restart) const;
    /**
       @brief Sets tolerances for the eigenvalue solver
       @param tol tolerance of the eigensolver
       @param max_its Maximum iterations of the eigensolver
       @note Negative values for PETSC_DECIDE
    */
    void SetTolerances(const RTVar tol, const int max_its);
    /**
       @brief Gets tolerances for the eigenvalue solver
       @param tol tolerance of the eigensolver
       @param max_its Maximum iterations of the eigensolver
    */
    void GetTolerances(RTVar &tol, int &max_its) const;
    //! Sets convergence test for the eigensolver
    void SetConvergenceTest(const EPSConv test);
    //! Gets convergence test for the eigensolver
    EPSConv GetConvergenceTest() const;
    //! Whether or not to compute the residual explicitly
    void SetTrueResidual(const bool opt);
    //! Sets eigensolver algorithm
    void SetType(const EPSType type);
    //! Gets eigensolver algorithm
    EPSType GetType() const;
    /**
       @brief Sets the number of eigenvalues to compute and the dimension of the subspace. 
       @param nev number of eigenvalues to compute
       @param ncv maximum dimension of the subspace
       @param mpd maximum dimension for the projected problem
       @note Set ncv and mpd to negative values for PETSC_DECIDE
     */
    void SetEPSDimensions(const int nev, const int ncv, const int mpd);
    /**
       @brief Gets the number of eigenvalues to compute and the dimension of the subspace. 
       @param nev number of eigenvalues to compute
       @param ncv maximum dimension of the subspace
       @param mpd maximum dimension for the projected problem
     */
    void GetEPSDimensions(int &nev, int &ncv, int &mpd) const;
    //!Sets verbosity level
    void SetVerbose(bool verbose);

    //!Sets how to invert the matrix obtained by shift-and invert
    void SetLinearSolver(const KSPSolver solver);
    /**
       @brief Sets tolerances regarding the linear solver
       @param rtol Relative tolerance for the linear solver
       @param atol Absolute tolerance for the linear solver
       @param dtol Divergence condition for the linear solver
       @param max_its Maximum iterations of the linear solver
       @note Set negative values for PETSC_DECIDE
    */
    void SetLinearSolverTol(const RTVar rtol, const RTVar atol, const RTVar dtol, const int max_its);
    //!Gets tolerances regarding the linear solver
    void GetLinearSolverTol(RTVar &rtol, RTVar &atol, RTVar &dtol, int &max_its);
    //!Sets preconditioner to be used and tolerance for zero pivot
    void SetPrecond(const Precond precond,
                    RTVar zero = std::numeric_limits<RTVar>::epsilon());

    //! Returns nullptr since pardiso is not used in this solver
    TPZPardisoSolver<TVar> *GetPardisoControlA() override
    {return nullptr;}
    TPZPardisoSolver<TVar> *GetPardisoControlB() override
    {return nullptr;}
  private:
    //! Actual solver implementation with SLEPc calls
    int SolveImpl(TPZVec<CTVar> &w,TPZFMatrix<CTVar> &eigenVectors,
                  bool calcVectors);
    //! Checks whether the matrices are in a correct format for SLEPc
    bool CheckMatrixTypes();
    //! Controls verbosity level
    bool fVerbose = true;
    //! Problem Type
    EPSProblemType fProbType{EPSProblemType::EPS_NOTSET};
    //! Which part of the spectrum to solve for
    EPSWhich fWhich{EPSWhich::EPS_LARGEST_MAGNITUDE};
    //! Whether to use the locking variant if using Krylov algorithm
    bool fLocking{false};
    //! Percentage of eigenvectors kept after restart if using Krylov algorithm
    RTVar fRestart{0.5};
    //! Tolerance of the eigensolver
    RTVar fEpsTol{-1};
    //! Maximum iterations of the eigensolver
    RTVar fEpsMaxIts{-1};
    //! Convergence test for the eigensolver
    EPSConv fConvTest{EPSConv::EPS_CONV_REL};
    //! Whether to compute true residual explicitly
    bool fTrueResidual{false};
    //! Eigensolver algorithm
    EPSType fEpsType{EPSType::KRYLOVSCHUR};
    //! the maximum dimension of the subspace to be used by the subsolve 
    int fNcv{-1};
    //! the maximum dimension allowed for the projected problem
    int fMpd{-1};
    //! Linear solver to be used
    KSPSolver fKsp{KSPSolver::PREONLY};
    //! Relative tolerance for the linear solver
    RTVar fKspRtol{-1};
    //! Absolute tolerance for the linear solver
    RTVar fKspAtol{-1};
    //! Divergence condition for the linear solver
    RTVar fKspDtol{-1};
    //! Maximum iterations of the linear solver
    RTVar fKspMaxIts{-1};
    //! Preconditioner for the linear solver
    Precond fPc{Precond::LU};
    //! Zero pivot tolerance for the preconditioner
    RTVar fPcZero{std::numeric_limits<RTVar>::epsilon()};
    //! Whether SLEPc has been initialized
    static bool fSlepcInit;
  };
}
#endif //WGMASLEPCHANDLER_H
