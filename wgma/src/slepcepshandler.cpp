//
// Created by Francisco Teixeira Orlandini on 12/13/17.
//
#include "slepcepshandler.hpp"

#include "pzysmp.h"
#include "TPZSimpleTimer.h"

#ifdef WGMA_USING_SLEPC
#include <slepceps.h>
#include <petsctime.h>
#include <petscksp.h>
#endif


/*******************
 *    GENERAL       *
 *******************/
namespace wgma::slepc{

  std::string SolverName(EPSType);
  
#ifdef WGMA_USING_SLEPC
  ::EPSProblemType ConvertProblemType(EPSProblemType);
  EPSProblemType ConvertProblemType(::EPSProblemType);

  ::EPSConv ConvertConv(EPSConv);
  EPSConv ConvertConv(::EPSConv);

 ::EPSType ConvertType(EPSType);
  EPSType ConvertType(::EPSType);

  ::EPSWhich ConvertWhich(EPSWhich);
  EPSWhich ConvertWhich(::EPSWhich);

  KSPSolver ConvertKSP(::KSPType in);

  ::KSPType ConvertKSP(KSPSolver in);

  Precond ConvertPrecond(PCType in);
  
  ::PCType ConvertPrecond(Precond in);
  
#endif

  template<class TVar>
  bool EPSHandler<TVar>::fSlepcInit = false;
  
  template<class TVar>
  void EPSHandler<TVar>::InitSLEPc(){
    //initialize SLEPc
#ifdef WGMA_USING_SLEPC    
    SlepcInitialize((int *)0, (char ***)0, (const char*)0,(const char*)0 );
    fSlepcInit = true;
#else
    std::cerr <<"WARNING:The EPSHandler module is only available if the wgma"
              <<"was configured with SLEPc library. Aborting..."<<std::endl;
    exit(-1); 
#endif
    
  }

  template<class TVar>
  void EPSHandler<TVar>::FinalizeSLEPc(){
#ifdef WGMA_USING_SLEPC    
    SlepcFinalize();
#else
    std::cerr <<"WARNING:The EPSHandler module is only available if the wgma"
              <<"was configured with SLEPc library. Aborting..."<<std::endl;
    exit(-1); 
#endif
    
  }
  
  template<class TVar>
  EPSHandler<TVar>::EPSHandler(){
#ifdef WGMA_USING_SLEPC    
    static_assert(std::is_same_v<TVar,PetscScalar>,
                  "Cannot create EPSHandler<T> for T!=PetscScalar\n"
                  "Check your configuration");
#else
    std::cerr <<"WARNING:The EPSHandler module is only available if the wgma"
              <<"was configured with SLEPc library. Aborting..."<<std::endl;
    exit(-1); 
#endif
  }
  
  template<class TVar>
  EPSHandler<TVar>* EPSHandler<TVar>::Clone() const {
    return (EPSHandler*) new EPSHandler(*this);
  }

  template<class TVar>
  void EPSHandler<TVar>::SetNEigenpairs(int n) {
    std::cout<<__PRETTY_FUNCTION__
             <<"\nNote: setting ncv = 0 and mpd = 0.\n"
             <<"Call instead EPSHandler<T>::SetEPSDimensions";
    SetEPSDimensions(n, 0, 0);
  }
  
  template<class TVar>
  bool EPSHandler<TVar>::CheckMatrixTypes(){
    auto sparseA =
      TPZAutoPointerDynamicCast<TPZFYsmpMatrix<TVar>>(this->MatrixA());
    auto sparseB =
      TPZAutoPointerDynamicCast<TPZFYsmpMatrix<TVar>>(this->MatrixB());

    if(!sparseA ||
       (this->fIsGeneralised && !sparseB)){
      PZError<<__PRETTY_FUNCTION__
             <<"\nIncompatible matrix types. Use:"
             <<"\t\tTPZSpStructMatrix<CSTATE>"
             <<"for obtaining a"
             <<"\t\tTPZFYsmpMatrix<CSTATE>"
             <<"sparse matrix. Aborting...\n";
      return false;
    }
    return true;
  }
  
  template<class TVar>
  int EPSHandler<TVar>::SolveEigenProblem(TPZVec<CTVar> &w, TPZFMatrix<CTVar> &ev){
    return SolveImpl(w, ev, true);
  }
  template<class TVar>
  int EPSHandler<TVar>::SolveEigenProblem(TPZVec<CTVar> &w){
    TPZFMatrix<CTVar> ev;
    return SolveImpl(w, ev, false);
  }

  template<class TVar>
  int EPSHandler<TVar>::SolveGeneralisedEigenProblem(TPZVec<CTVar> &w,
                                               TPZFMatrix<CTVar> &ev){
    return SolveImpl(w, ev, true);
  }

  template<class TVar>
  int EPSHandler<TVar>::SolveGeneralisedEigenProblem(TPZVec<CTVar> &w){
    TPZFMatrix<CTVar> ev;
    return SolveImpl(w, ev, false);
  }
  
  template<class TVar>
  int EPSHandler<TVar>::SolveImpl(TPZVec <CTVar> &w, TPZFMatrix <CTVar>& eigenVectors,
                                  bool calcVectors){
#ifdef WGMA_USING_SLEPC
    CheckMatrixTypes();
    auto &pzA = dynamic_cast<TPZFYsmpMatrix<TVar>&>(this->MatrixA().operator*());
    auto &pzB = dynamic_cast<TPZFYsmpMatrix<TVar>&>(this->MatrixB().operator*());

    if(! fSlepcInit) {InitSLEPc();}
    
    auto CreatePetscMat = [](TPZFYsmpMatrix<TVar> &pzmat,
                             Mat &mat,
                             PetscInt *&ia,
                             PetscInt *&ja,
                             PetscScalar *&aa) ->int{
      const int nRows = pzmat.Rows();
      const int nCols = pzmat.Cols();
      TPZVec<int64_t> I, J;
      TPZVec<TVar> A;
      pzmat.GetData(I,J,A);
      
      PetscErrorCode ierr;
      ierr = PetscMalloc1(I.size(),&ia);CHKERRQ(ierr);
      ierr = PetscMalloc1(J.size(),&ja);CHKERRQ(ierr);
      ierr = PetscMalloc1(J.size(),&aa);CHKERRQ(ierr);

      for (int j = 0; j < I.size(); ++j) {
        ia[j]=I[j];
      }
      for (int j = 0; j < J.size(); ++j) {
        ja[j]=J[j];
        aa[j]=A[j];
      }
      ierr = MatCreateSeqAIJWithArrays(MPI_COMM_WORLD,nRows,nCols,ia,ja,
                                       aa,&mat);
      CHKERRQ(ierr);
      return 0;
    };
    Mat petscA{nullptr}, petscB{nullptr};
    PetscInt *iaP{nullptr}, *jaP{nullptr}, *ibP{nullptr}, *jbP{nullptr};
    PetscScalar *aaP{nullptr}, *abP{nullptr};
    {
      TPZSimpleTimer create("CreatePetscMat");
      std::cout<<"Creating PETSc Amat...";
      CreatePetscMat(pzA,petscA, iaP, jaP, aaP);
      std::cout<<"Created!"<<std::endl;
      //if it is not generalised, we just need petscB to be nullptr
      if(this->fIsGeneralised){
        std::cout<<"Creating PETSc Bmat...";
        CreatePetscMat(pzB,petscB, ibP, jbP, abP);
        std::cout<<"Created!"<<std::endl;
      }
    }
    /**
       SET UP PC, KSP, ST and EPS based on user options
     **/
    PetscErrorCode ierr;
    ::PC pc;
    ::KSP ksp;
    ::ST st;
    ::EPS eps;
    //PC settings
    {
      ierr = PCCreate(PETSC_COMM_WORLD, &pc);

      const PCType pc_type = ConvertPrecond(fPc);
      ierr = PCSetType(pc, pc_type);

      const STATE pc_zero = fPcZero > 0 ? fPcZero : PETSC_DECIDE;
      ierr = PCFactorSetZeroPivot(pc,pc_zero);
#ifdef PETSC_HAVE_MUMPS
      if(!strcmp(pc_type,"lu")||!strcmp(pc_type,"cholesky")){
        ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
      }
#endif
      CHKERRQ(ierr);
      //KSP settings
    
      ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
      ierr = KSPSetPC(ksp, pc);
      const ::KSPType ksp_type = ConvertKSP(fKsp);
      ierr = KSPSetType(ksp,ksp_type);
      const STATE ksp_rtol = fKspRtol > 0 ? fKspRtol : PETSC_DEFAULT;
      const STATE ksp_atol = fKspAtol > 0 ? fKspAtol : PETSC_DEFAULT;
      const STATE ksp_dtol = fKspDtol > 0 ? fKspDtol : PETSC_DEFAULT;
      const int ksp_max_ints = fKspMaxIts > 0 ? fKspMaxIts : PETSC_DEFAULT;
      ierr = KSPSetTolerances(ksp, ksp_rtol, ksp_atol, ksp_dtol, ksp_max_ints);
      CHKERRQ(ierr);
      //ST settings
      
      STCreate(PETSC_COMM_WORLD, &st);
      STSetKSP(st, ksp);
      const ::STType st_type = STSINVERT;
      ierr = STSetType(st, st_type);
      CHKERRQ(ierr);
      
      //EPS settings
      ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
      ierr = EPSSetST(eps, st);
      const ::EPSType eps_type = ConvertType(fEpsType);
      ierr = EPSSetType(eps, eps_type);
      const ::EPSProblemType eps_prob_type = ConvertProblemType(fProbType);
      ierr = EPSSetProblemType(eps, eps_prob_type);
      const STATE eps_tol = fEpsTol > 0 ? fEpsTol : PETSC_DEFAULT;
      const STATE eps_max_its = fEpsMaxIts > 0 ? fEpsMaxIts : PETSC_DEFAULT;
      ierr = EPSSetTolerances(eps, eps_tol, eps_max_its);
      const ::EPSConv eps_conv = ConvertConv(fConvTest);
      ierr = EPSSetConvergenceTest(eps, eps_conv);
      const ::EPSWhich eps_which = ConvertWhich(fWhich);
      ierr = EPSSetWhichEigenpairs(eps, eps_which);
      const auto eps_target = this->fTarget;
      ierr = EPSSetTarget(eps, eps_target);
      const PetscBool eps_true_res = fTrueResidual ? PETSC_TRUE : PETSC_FALSE;
      ierr = EPSSetTrueResidual(eps, eps_true_res);
      const PetscInt nev = this->fNEigenpairs;
      const PetscInt ncv = fNcv > 0 ? fNcv : PETSC_DEFAULT;
      const PetscInt mpd = fMpd > 0 ? fMpd : PETSC_DEFAULT;
      ierr = EPSSetDimensions(eps, nev, ncv, mpd);
      if(!strcmp(eps_type,EPSKRYLOVSCHUR)){
        const PetscBool locking = fLocking ? PETSC_TRUE : PETSC_FALSE;
        const PetscReal restart = fRestart > 0 ? fRestart : PETSC_DEFAULT;
        ierr = EPSKrylovSchurSetLocking(eps, locking);
        ierr = EPSKrylovSchurSetRestart(eps, restart);
      }
      CHKERRQ(ierr);
    }
    

    /**
       BEGINNING OF THE OPERATIONS
     */
    
    EPSSetOperators(eps, petscA, petscB);

    {
      TPZSimpleTimer setup("EPSSetUp");
      EPSSetUp(eps);
    }

    if(fVerbose){
      EPSView(eps,PETSC_VIEWER_STDOUT_WORLD);
      PetscViewerAndFormat *vf;
      PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);
      EPSMonitorSet(eps,
                    (PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,
                                        PetscReal*,PetscInt,void*))EPSMonitorFirst,vf,
                    (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
    }
    /*****************************
     *  SOLVE
     *****************************/
    {
      TPZSimpleTimer solver("EPSSolve");
      PetscLogDouble t1,t2;
      ierr = PetscTime(&t1);CHKERRQ(ierr);
      ierr = EPSSolve(eps);CHKERRQ(ierr);
      ierr = PetscTime(&t2);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Elapsed Time in EPSSolve: %f\n",t2-t1);CHKERRQ(ierr);
      
    }
    /*
      Optional: Get some information from the solver and display it
    */
    {
      PetscInt its,lits,maxit,nev;
      PetscReal tol;
      ::EPSType type;
      ::KSP ksp;
      ::ST st;
      ierr = EPSGetIterationNumber(eps, &its);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
      ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
      ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
      ierr = KSPGetTotalIterations(ksp, &lits);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %D\n",lits);CHKERRQ(ierr);
      ierr = EPSGetType(eps, &type);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
      ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n", nev);CHKERRQ(ierr);
      ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);
    }
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Display solution and clean up
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /*
      Show detailed info unless -terse option is given by user
    */
    ::EPSConv eps_conv_test;
    EPSErrorType eps_error_type = EPS_ERROR_RELATIVE;
    EPSGetConvergenceTest(eps, &eps_conv_test);
    switch(eps_conv_test){
    case EPS_CONV_ABS:
      eps_error_type = EPS_ERROR_ABSOLUTE;
      break;
    case EPS_CONV_REL:
      eps_error_type = EPS_ERROR_RELATIVE;
      break;
    case EPS_CONV_NORM:
      eps_error_type = EPS_ERROR_BACKWARD;
      break;
    case EPS_CONV_USER:
      PZError<<__PRETTY_FUNCTION__
             <<"\n EPS_CONV_USER not supported. Aborting...\n";
      DebugStop();
    }
    if (fVerbose) {
      ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
      ierr = EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      ierr = EPSErrorView(eps,eps_error_type,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    } else {
      ierr = EPSErrorView(eps,eps_error_type,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
    PetscInt nconv;
    EPSGetConverged(eps, &nconv);
    const PetscInt nEigen = this->fNEigenpairs > nconv ? nconv : this->fNEigenpairs;
    w.Resize(nEigen);

    //let us sort the eigenvalues
    TPZManVector<int,20> indices;
    {
      PetscScalar eigr{0}, eigi{0};
      for(int i = 0; i < nEigen; i++){
        EPSGetEigenvalue(eps,i, &eigr, &eigi);
        if constexpr(std::is_same_v<PetscScalar,STATE>){
          w[i] = eigr + 1i * eigi;
        }else{
          w[i] = eigr;
        }
        
      }
      this->SortEigenvalues(w,indices);
    }

    if(!calcVectors) return 0;
    
    eigenVectors.Resize(pzA.Rows(),nEigen);
    for (int i = 0; i < nEigen; ++i) {
      auto il = indices[i];
      if constexpr(std::is_same_v<PetscScalar,STATE>){
        Vec eigVecRe, eigVecIm;
        PetscScalar *eigVecReArray, *eigVecImArray;
        
        MatCreateVecs(petscA,&eigVecRe,nullptr);
        MatCreateVecs(petscA,&eigVecIm,nullptr);
        EPSGetEigenvector(eps,il,eigVecRe,eigVecIm);
        VecGetArray(eigVecRe,&eigVecReArray);
        VecGetArray(eigVecIm,&eigVecImArray);
        for (int j = 0; j < pzA.Rows(); ++j) {
          eigenVectors(j,i) = eigVecReArray[j] + 1i * eigVecImArray[j];
        }
        VecRestoreArray(eigVecRe,&eigVecReArray);
        VecRestoreArray(eigVecIm,&eigVecImArray);
      }else{
        Vec eigVec;
        PetscScalar  *eigVecArray;
        MatCreateVecs(petscA,&eigVec,nullptr);
        EPSGetEigenvector(eps,il,eigVec,nullptr);
        VecGetArray(eigVec,&eigVecArray);
        for (int j = 0; j < pzA.Rows(); ++j) {
          eigenVectors(j,i) = eigVecArray[j];
        }
        VecRestoreArray(eigVec,&eigVecArray);
      }
    }

    EPSDestroy(&eps);
    
    MatDestroy(&petscA);
    MatDestroy(&petscB);


    if(iaP) PetscFree(iaP);
    if(jaP) PetscFree(jaP);
    if(aaP) PetscFree(aaP);
    if(ibP) PetscFree(ibP);
    if(jbP) PetscFree(jbP);
    if(abP) PetscFree(abP);
    
    return 0;
#endif
    return -1;
  }

  template<class TVar>
  void EPSHandler<TVar>::SetProblemType(const EPSProblemType eps_problem)
  {
    fProbType = eps_problem;
    this->SetAsGeneralised(!(eps_problem == EPSProblemType::EPS_HEP ||
                             eps_problem == EPSProblemType::EPS_NHEP));
  }
  
  template<class TVar>
  EPSProblemType EPSHandler<TVar>::GetProblemType() const
  {
    return fProbType;
  }
  
  template<class TVar>
  void EPSHandler<TVar>::SetAsGeneralised(bool isGeneralised){
    TPZEigenSolver<TVar>::SetAsGeneralised(isGeneralised);
    if(fProbType == EPSProblemType::EPS_NOTSET){
      std::cout<<__PRETTY_FUNCTION__
               <<"\nWARNING: call SetProblemType instead"
               <<std::endl;
      if(isGeneralised) fProbType = EPSProblemType::EPS_GNHEP;
      else fProbType = EPSProblemType::EPS_NHEP;
    }
  }

  template<class TVar>
  void EPSHandler<TVar>::SetWhichEigenpairs (const EPSWhich which){
    fWhich = which;
  }

  template<class TVar>
  EPSWhich EPSHandler<TVar>::GetWhichEigenpairs () const {
    return fWhich;
  }

  template<class TVar>
  void EPSHandler<TVar>::SetKrylovOptions (const bool pLocking, const RTVar restart){
    if(fEpsType != EPSType::KRYLOVSCHUR){
      std::cout<<__PRETTY_FUNCTION__
               <<"\n Solver is not KRYLOVSCHUR"
               <<"\teps_type = "<<SolverName(fEpsType)
               <<"\nThese options may be meaningless.\n"; 
    }
    fLocking = pLocking;
    fRestart = restart;
  }

  template<class TVar>
  void EPSHandler<TVar>::GetKrylovOptions (bool &pLocking, RTVar &restart) const{
    if(fEpsType != EPSType::KRYLOVSCHUR){
      std::cout<<__PRETTY_FUNCTION__
               <<"\n Solver is not KRYLOVSCHUR"
               <<"\teps_type = "<<SolverName(fEpsType)
               <<"\nThese options may be meaningless.\n"; 
    }
    pLocking = fLocking;
    restart = fRestart;
  }

  template<class TVar>
  void EPSHandler<TVar>::SetTolerances(const RTVar t, const int m_its) {
    fEpsTol = t;
    fEpsMaxIts = m_its;
  }

  template<class TVar>
  void EPSHandler<TVar>::GetTolerances(RTVar &tol, int &max_its) const {
    tol = fEpsTol;
    max_its = fEpsMaxIts;
  }

  template<class TVar>
  void EPSHandler<TVar>::SetConvergenceTest(const EPSConv test) {
    fConvTest = test;
  }

  template<class TVar>
  EPSConv EPSHandler<TVar>::GetConvergenceTest() const {
    return fConvTest;
  }

  template<class TVar>
  void EPSHandler<TVar>::SetTrueResidual(const bool pOpt) {
    fTrueResidual = pOpt;
  }

  template<class TVar>
  void EPSHandler<TVar>::SetType(const EPSType type) {
    fEpsType = type;
  }
  template<class TVar>
  EPSType EPSHandler<TVar>::GetType() const {
    return fEpsType;
  }

  template<class TVar>
  void EPSHandler<TVar>::SetEPSDimensions(const int nev,
                                          const int nc,
                                          const int mp) {    
    this->fNEigenpairs = nev;
    fNcv = nc;
    fMpd = mp;
  }

  template<class TVar>
  void EPSHandler<TVar>::GetEPSDimensions(int &nev,
                                          int &ncv,
                                          int &mpd) const {
    nev = this->fNEigenpairs;
    ncv = fNcv;
    mpd = fMpd;
  }

  template<class TVar>
  void EPSHandler<TVar>::SetVerbose(bool fVerbose) {
    EPSHandler::fVerbose = fVerbose;
  }

  template<class TVar>
  void EPSHandler<TVar>::SetLinearSolver(const KSPSolver solver){
    fKsp = solver;
  }

  template<class TVar>
  void EPSHandler<TVar>::SetLinearSolverTol(const RTVar r, const RTVar a,
                                            const RTVar d, const int m){
    fKspRtol = r;
    fKspAtol = a;
    fKspDtol = d;
    fKspMaxIts = m;
  }

  template<class TVar>
  void EPSHandler<TVar>::GetLinearSolverTol(RTVar &rtol, RTVar &atol,
                                            RTVar &dtol, int &max_its){
    rtol = fKspRtol;
    atol = fKspAtol;
    dtol = fKspDtol;
    max_its = fKspMaxIts;
  }
  
  template<class TVar>
  void EPSHandler<TVar>::SetPrecond(const Precond pre, RTVar zero){
    fPc = pre;
    fPcZero = zero;
  }

  std::string SolverName(EPSType type){
    switch(type){
    case EPSType::POWER: return "POWER";
    case EPSType::SUBSPACE: return "SUBSPACE";
    case EPSType::ARNOLDI: return "ARNOLDI";
    case EPSType::LANCZOS: return "LANCZOS";
    case EPSType::KRYLOVSCHUR: return "KRYLOVSCHUR";
    case EPSType::GD: return "GD";
    case EPSType::JD: return "JD";
    }
  };
  
#ifdef WGMA_USING_SLEPC
  ::EPSProblemType ConvertProblemType(EPSProblemType in)
  {
    switch(in){
    case EPSProblemType::EPS_HEP: return EPS_HEP;
    case EPSProblemType::EPS_GHEP: return EPS_GHEP;
    case EPSProblemType::EPS_NHEP: return EPS_NHEP;
    case EPSProblemType::EPS_GNHEP: return EPS_GNHEP;
    case EPSProblemType::EPS_PGNHEP: return EPS_PGNHEP;
    case EPSProblemType::EPS_GHIEP: return EPS_GHIEP;
    default:
      unreachable();
    }
    DebugStop();
  }
  
  EPSProblemType ConvertProblemType(::EPSProblemType in)
  {
    switch(in){
    case EPS_HEP: return EPSProblemType::EPS_HEP;
    case EPS_GHEP: return EPSProblemType::EPS_GHEP;
    case EPS_NHEP: return EPSProblemType::EPS_NHEP;
    case EPS_GNHEP: return EPSProblemType::EPS_GNHEP;
    case EPS_PGNHEP: return EPSProblemType::EPS_PGNHEP;
    case EPS_GHIEP: return EPSProblemType::EPS_GHIEP;
    }
  }

  ::EPSConv ConvertConv(EPSConv in)
  {
    switch(in){
    case EPSConv::EPS_CONV_ABS: return EPS_CONV_ABS;
    case EPSConv::EPS_CONV_REL: return EPS_CONV_REL;
    case EPSConv::EPS_CONV_NORM: return EPS_CONV_NORM;
    }
  }
  EPSConv ConvertConv(::EPSConv in)
  {
    switch(in){
    case EPS_CONV_ABS: return EPSConv::EPS_CONV_ABS;
    case EPS_CONV_REL: return EPSConv::EPS_CONV_REL;
    case EPS_CONV_NORM: return EPSConv::EPS_CONV_NORM;
    default:
      unreachable();
    }
  }

  ::EPSType ConvertType(EPSType in)
  {
    switch(in){
    case EPSType::POWER: return EPSPOWER;
    case EPSType::SUBSPACE: return EPSSUBSPACE;
    case EPSType::ARNOLDI: return EPSARNOLDI;
    case EPSType::LANCZOS: return EPSLANCZOS;
    case EPSType::KRYLOVSCHUR: return EPSKRYLOVSCHUR;
    case EPSType::GD: return EPSGD;
    case EPSType::JD: return EPSJD;
    default:
      unreachable();
    }
  }
  
  EPSType ConvertType(::EPSType in)
  {
    
    if(!strcmp(in,EPSPOWER)) return EPSType::POWER;
    else if(!strcmp(in,EPSSUBSPACE)) return EPSType::SUBSPACE;
    else if(!strcmp(in,EPSARNOLDI)) return EPSType::ARNOLDI;
    else if(!strcmp(in,EPSLANCZOS)) return EPSType::LANCZOS;
    else if(!strcmp(in,EPSKRYLOVSCHUR)) return EPSType::KRYLOVSCHUR;
    else if(!strcmp(in,EPSGD)) return EPSType::GD;
    else if(!strcmp(in,EPSJD)) return EPSType::JD;
    else{
      PZError<<__PRETTY_FUNCTION__
             <<"\nUnsupported type.\nAborting...\n";
      DebugStop();
    }
    return EPSType::POWER;
  }

  ::EPSWhich ConvertWhich(EPSWhich in)
  {
    switch(in){
    case EPSWhich::EPS_LARGEST_MAGNITUDE: return EPS_LARGEST_MAGNITUDE;
    case EPSWhich::EPS_SMALLEST_MAGNITUDE: return EPS_SMALLEST_MAGNITUDE;
    case EPSWhich::EPS_LARGEST_REAL: return EPS_LARGEST_REAL;
    case EPSWhich::EPS_SMALLEST_REAL: return EPS_SMALLEST_REAL;
    case EPSWhich::EPS_LARGEST_IMAGINARY: return EPS_LARGEST_IMAGINARY;
    case EPSWhich::EPS_SMALLEST_IMAGINARY: return EPS_SMALLEST_IMAGINARY;
    case EPSWhich::EPS_TARGET_MAGNITUDE: return EPS_TARGET_MAGNITUDE;
    case EPSWhich::EPS_TARGET_REAL: return EPS_TARGET_REAL;
    case EPSWhich::EPS_TARGET_IMAGINARY: return EPS_TARGET_IMAGINARY;
    default:
      unreachable();
    }
  }
  
  EPSWhich ConvertWhich(::EPSWhich in)
  {
    switch(in){
    case EPS_LARGEST_MAGNITUDE: return EPSWhich::EPS_LARGEST_MAGNITUDE;
    case EPS_SMALLEST_MAGNITUDE: return EPSWhich::EPS_SMALLEST_MAGNITUDE;
    case EPS_LARGEST_REAL: return EPSWhich::EPS_LARGEST_REAL;
    case EPS_SMALLEST_REAL: return EPSWhich::EPS_SMALLEST_REAL;
    case EPS_LARGEST_IMAGINARY: return EPSWhich::EPS_LARGEST_IMAGINARY;
    case EPS_SMALLEST_IMAGINARY: return EPSWhich::EPS_SMALLEST_IMAGINARY;
    case EPS_TARGET_MAGNITUDE: return EPSWhich::EPS_TARGET_MAGNITUDE;
    case EPS_TARGET_REAL: return EPSWhich::EPS_TARGET_REAL;
    case EPS_TARGET_IMAGINARY: return EPSWhich::EPS_TARGET_IMAGINARY;
    default:
      unreachable();
    }
  }

  KSPSolver ConvertKSP(::KSPType in)
  {
    
    if(!strcmp(in,KSPRICHARDSON)) return KSPSolver::RICHARDSON;
    else if(!strcmp(in,KSPCHEBYSHEV)) return KSPSolver::CHEBYSHEV;
    else if(!strcmp(in,KSPCG)) return KSPSolver::CG;
    else if(!strcmp(in,KSPGROPPCG)) return KSPSolver::GROPPCG;
    else if(!strcmp(in,KSPPIPECG)) return KSPSolver::PIPECG;
    else if(!strcmp(in,KSPPIPECGRR)) return KSPSolver::PIPECGRR;
    else if(!strcmp(in,KSPPIPELCG)) return KSPSolver::PIPELCG;
    else if(!strcmp(in,KSPPIPEPRCG)) return KSPSolver::PIPEPRCG;
    else if(!strcmp(in,KSPPIPECG2)) return KSPSolver::PIPECG2;
    else if(!strcmp(in,KSPCGNE)) return KSPSolver::CGNE;
    else if(!strcmp(in,KSPNASH)) return KSPSolver::NASH;
    else if(!strcmp(in,KSPSTCG)) return KSPSolver::STCG;
    else if(!strcmp(in,KSPGLTR)) return KSPSolver::GLTR;
    else if(!strcmp(in,KSPFCG)) return KSPSolver::FCG;
    else if(!strcmp(in,KSPPIPEFCG)) return KSPSolver::PIPEFCG;
    else if(!strcmp(in,KSPGMRES)) return KSPSolver::GMRES;
    else if(!strcmp(in,KSPPIPEFGMRES)) return KSPSolver::PIPEFGMRES;
    else if(!strcmp(in,KSPFGMRES)) return KSPSolver::FGMRES;
    else if(!strcmp(in,KSPLGMRES)) return KSPSolver::LGMRES;
    else if(!strcmp(in,KSPDGMRES)) return KSPSolver::DGMRES;
    else if(!strcmp(in,KSPPGMRES)) return KSPSolver::PGMRES;
    else if(!strcmp(in,KSPTCQMR)) return KSPSolver::TCQMR;
    else if(!strcmp(in,KSPBCGS)) return KSPSolver::BCGS;
    else if(!strcmp(in,KSPIBCGS)) return KSPSolver::IBCGS;
    else if(!strcmp(in,KSPFBCGS)) return KSPSolver::FBCGS;
    else if(!strcmp(in,KSPFBCGSR)) return KSPSolver::FBCGSR;
    else if(!strcmp(in,KSPBCGSL)) return KSPSolver::BCGSL;
    else if(!strcmp(in,KSPPIPEBCGS)) return KSPSolver::PIPEBCGS;
    else if(!strcmp(in,KSPCGS)) return KSPSolver::CGS;
    else if(!strcmp(in,KSPTFQMR)) return KSPSolver::TFQMR;
    else if(!strcmp(in,KSPCR)) return KSPSolver::CR;
    else if(!strcmp(in,KSPPIPECR)) return KSPSolver::PIPECR;
    else if(!strcmp(in,KSPLSQR)) return KSPSolver::LSQR;
    else if(!strcmp(in,KSPPREONLY)) return KSPSolver::PREONLY;
    else if(!strcmp(in,KSPQCG)) return KSPSolver::QCG;
    else if(!strcmp(in,KSPBICG)) return KSPSolver::BICG;
    else if(!strcmp(in,KSPMINRES)) return KSPSolver::MINRES;
    else if(!strcmp(in,KSPSYMMLQ)) return KSPSolver::SYMMLQ;
    else if(!strcmp(in,KSPLCD)) return KSPSolver::LCD;
    else if(!strcmp(in,KSPPYTHON)) return KSPSolver::PYTHON;
    else if(!strcmp(in,KSPGCR)) return KSPSolver::GCR;
    else if(!strcmp(in,KSPPIPEGCR)) return KSPSolver::PIPEGCR;
    else if(!strcmp(in,KSPTSIRM)) return KSPSolver::TSIRM;
    else if(!strcmp(in,KSPCGLS)) return KSPSolver::CGLS;
    else if(!strcmp(in,KSPFETIDP)) return KSPSolver::FETIDP;
    else if(!strcmp(in,KSPHPDDM)) return KSPSolver::HPDDM;
    else{
      PZError<<__PRETTY_FUNCTION__
             <<"\nUnsupported type.\nAborting...\n";
      DebugStop();
    }
    return KSPSolver::PREONLY;
  }

  ::KSPType ConvertKSP(KSPSolver in)
  {
    switch(in){
    case KSPSolver::RICHARDSON: return KSPRICHARDSON;
    case KSPSolver::CHEBYSHEV: return KSPCHEBYSHEV;
    case KSPSolver::CG: return KSPCG;
    case KSPSolver::GROPPCG: return KSPGROPPCG;
    case KSPSolver::PIPECG: return KSPPIPECG;
    case KSPSolver::PIPECGRR: return KSPPIPECGRR;
    case KSPSolver::PIPELCG: return KSPPIPELCG;
    case KSPSolver::PIPEPRCG: return KSPPIPEPRCG;
    case KSPSolver::PIPECG2: return KSPPIPECG2;
    case KSPSolver::CGNE: return KSPCGNE;
    case KSPSolver::NASH: return KSPNASH;
    case KSPSolver::STCG: return KSPSTCG;
    case KSPSolver::GLTR: return KSPGLTR;
    case KSPSolver::FCG: return KSPFCG;
    case KSPSolver::PIPEFCG: return KSPPIPEFCG;
    case KSPSolver::GMRES: return KSPGMRES;
    case KSPSolver::PIPEFGMRES: return KSPPIPEFGMRES;
    case KSPSolver::FGMRES: return KSPFGMRES;
    case KSPSolver::LGMRES: return KSPLGMRES;
    case KSPSolver::DGMRES: return KSPDGMRES;
    case KSPSolver::PGMRES: return KSPPGMRES;
    case KSPSolver::TCQMR: return KSPTCQMR;
    case KSPSolver::BCGS: return KSPBCGS;
    case KSPSolver::IBCGS: return KSPIBCGS;
    case KSPSolver::FBCGS: return KSPFBCGS;
    case KSPSolver::FBCGSR: return KSPFBCGSR;
    case KSPSolver::BCGSL: return KSPBCGSL;
    case KSPSolver::PIPEBCGS: return KSPPIPEBCGS;
    case KSPSolver::CGS: return KSPCGS;
    case KSPSolver::TFQMR: return KSPTFQMR;
    case KSPSolver::CR: return KSPCR;
    case KSPSolver::PIPECR: return KSPPIPECR;
    case KSPSolver::LSQR: return KSPLSQR;
    case KSPSolver::PREONLY: return KSPPREONLY;
    case KSPSolver::QCG: return KSPQCG;
    case KSPSolver::BICG: return KSPBICG;
    case KSPSolver::MINRES: return KSPMINRES;
    case KSPSolver::SYMMLQ: return KSPSYMMLQ;
    case KSPSolver::LCD: return KSPLCD;
    case KSPSolver::PYTHON: return KSPPYTHON;
    case KSPSolver::GCR: return KSPGCR;
    case KSPSolver::PIPEGCR: return KSPPIPEGCR;
    case KSPSolver::TSIRM: return KSPTSIRM;
    case KSPSolver::CGLS: return KSPCGLS;
    case KSPSolver::FETIDP: return KSPFETIDP;
    case KSPSolver::HPDDM: return KSPHPDDM;
    }
  }


  Precond ConvertPrecond(PCType in){
    
    if(!strcmp(in,PCNONE)) return Precond::NONE;
    else if(!strcmp(in,PCJACOBI)) return Precond::JACOBI;
    else if(!strcmp(in,PCSOR)) return Precond::SOR;
    else if(!strcmp(in,PCLU)) return Precond::LU;
    else if(!strcmp(in,PCSHELL)) return Precond::SHELL;
    else if(!strcmp(in,PCBJACOBI)) return Precond::BJACOBI;
    else if(!strcmp(in,PCMG)) return Precond::MG;
    else if(!strcmp(in,PCEISENSTAT)) return Precond::EISENSTAT;
    else if(!strcmp(in,PCILU)) return Precond::ILU;
    else if(!strcmp(in,PCICC)) return Precond::ICC;
    else if(!strcmp(in,PCASM)) return Precond::ASM;
    else if(!strcmp(in,PCGASM)) return Precond::GASM;
    else if(!strcmp(in,PCKSP)) return Precond::KSP;
    else if(!strcmp(in,PCREDUNDANT)) return Precond::REDUNDANT;
    else if(!strcmp(in,PCCHOLESKY)) return Precond::CHOLESKY;
    else if(!strcmp(in,PCPBJACOBI)) return Precond::PBJACOBI;
    else if(!strcmp(in,PCVPBJACOBI)) return Precond::VPBJACOBI;
    else if(!strcmp(in,PCSVD)) return Precond::SVD;
    else if(!strcmp(in,PCBDDC)) return Precond::BDDC;
    else if(!strcmp(in,PCKACZMARZ)) return Precond::KACZMARZ;
    else if(!strcmp(in,PCTELESCOPE)) return Precond::TELESCOPE;
    else if(!strcmp(in,PCPATCH)) return Precond::PATCH;
    else if(!strcmp(in,PCLMVM)) return Precond::LMVM;
    else if(!strcmp(in,PCHMG)) return Precond::HMG;
    else if(!strcmp(in,PCDEFLATION)) return Precond::DEFLATION;
    else{
      PZError<<__PRETTY_FUNCTION__
             <<"\nUnsupported type.\nAborting...\n";
      DebugStop();
    }
    return Precond::NONE;
  }
  
  ::PCType ConvertPrecond(Precond in){
    switch(in){
    case Precond::NONE: return PCNONE;
    case Precond::JACOBI: return PCJACOBI;
    case Precond::SOR: return PCSOR;
    case Precond::LU: return PCLU;
    case Precond::SHELL: return PCSHELL;
    case Precond::BJACOBI: return PCBJACOBI;
    case Precond::MG: return PCMG;
    case Precond::EISENSTAT: return PCEISENSTAT;
    case Precond::ILU: return PCILU;
    case Precond::ICC: return PCICC;
    case Precond::ASM: return PCASM;
    case Precond::GASM: return PCGASM;
    case Precond::KSP: return PCKSP;
    case Precond::REDUNDANT: return PCREDUNDANT;
    case Precond::CHOLESKY: return PCCHOLESKY;
    case Precond::PBJACOBI: return PCPBJACOBI;
    case Precond::VPBJACOBI: return PCVPBJACOBI;
    case Precond::SVD: return PCSVD;
    case Precond::BDDC: return PCBDDC;
    case Precond::KACZMARZ: return PCKACZMARZ;
    case Precond::TELESCOPE: return PCTELESCOPE;
    case Precond::PATCH: return PCPATCH;
    case Precond::LMVM: return PCLMVM;
    case Precond::HMG: return PCHMG;
    case Precond::DEFLATION: return PCDEFLATION;
    }
  }

#ifdef WGMA_PETSC_CPLX
  template class EPSHandler<CSTATE>;
#else
  template class EPSHandler<STATE>;
#endif

#else
  template class EPSHandler<CSTATE>;
  template class EPSHandler<STATE>;
#endif
}
