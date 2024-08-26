//
// Created by Francisco Teixeira Orlandini on 12/13/17.
//
#include "slepcepshandler.hpp"

#include "TPZYSMPMatrix.h"
#include "TPZSimpleTimer.h"

#ifdef WGMA_USING_SLEPC
#include <slepceps.h>
#include <petsctime.h>
#include <petscksp.h>
//PETSC_SUCCESS was introduced in PETSC 3.19
#ifndef PETSC_SUCCESS
#define PETSC_SUCCESS 0
#endif
#endif


using namespace std::complex_literals;

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

  ::EPSWhich ConvertWhich(TPZEigenSort);

  KSPSolver ConvertKSP(::KSPType in);

  ::KSPType ConvertKSP(KSPSolver in);

  PC ConvertPrecond(PCType in);
  
  ::PCType ConvertPrecond(PC in);
  
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
      PetscErrorCode ierr;
#ifdef WGMA_PETSC_64BIT
      pzmat.GetData(ia,ja,aa);
#else
      TPZVec<int64_t> I, J;
      TPZVec<TVar> A;
      pzmat.GetData(I,J,A);
      
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
#endif
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
    /*
      this std::function is used to pass the TPZEigenSolver<TVar>::UserSortingFunc()
      to the EPSSetEigenvalueComparison function.

      Due to scope, it needs to be here
     */
    std::function<bool(CTVar,CTVar)> stdfunc{nullptr};
    
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
      if(!strcmp(pc_type,"lu")||!strcmp(pc_type,"cholesky")){
#ifdef PETSC_HAVE_MKL_PARDISO
        ierr = PCFactorSetMatSolverType(pc,MATSOLVERMKL_PARDISO);
#elif defined(PETSC_HAVE_MUMPS)
        ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
#endif
      }
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
      const ::EPSWhich eps_which = ConvertWhich(this->EigenSorting());
      ierr = EPSSetWhichEigenpairs(eps, eps_which);

      if(eps_which == ::EPSWhich::EPS_WHICH_USER){
        stdfunc = this->UserSortingFunc();
        if(stdfunc==nullptr){
          DebugStop();
        }
        ierr = EPSSetEigenvalueComparison(eps,
                                          [](PetscScalar ar, PetscScalar ai,
                                             PetscScalar br, PetscScalar bi,
                                             PetscInt *res, void *ctx){
                                            auto locfunc
                                              = *(std::function<bool(CTVar,CTVar)>*)(ctx);
                                            const bool bool_res = locfunc(ar,br);
                                            *res = bool_res ? -1 : 1;
                                            PetscFunctionReturn(PETSC_SUCCESS);
                                          },&stdfunc);
      }
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
      ierr = EPSSetUp(eps);
      CHKERRQ(ierr);
    }

    Vec *initVecArray{nullptr};
    PetscScalar **initVecMem{nullptr};
    PetscInt nInitVec = this->fInitVec.Rows();
    {
      TPZSimpleTimer initvecs("Computing initial vectors");
      if(nInitVec>0){
        PetscBool is_computed;
        STGetTransform(st,&is_computed);
        if(!is_computed){
          DebugStop();
        }
        std::cout<<"Creating PETSc initial vector...";
        const int neq = pzA.Rows();
        if(neq!=this->fInitVec.Rows()){
          DebugStop();
        }
        initVecMem = new PetscScalar *[nInitVec];
        constexpr PetscInt blocksize{1};
        nInitVec = this->fInitVec.Cols();
        initVecArray = new Vec[nInitVec];
        for(int i = 0; i < nInitVec; i++){
          PetscMalloc1(neq,&initVecMem[i]);
          Vec x;
          VecCreateSeqWithArray(MPI_COMM_WORLD,blocksize, neq,
                                fInitVec.Elem()+neq*i, &x);
          initVecArray[i] = Vec{};
          VecCreateSeqWithArray(MPI_COMM_WORLD,blocksize,neq,
                                initVecMem[i], &initVecArray[i]);
          STMatSolve(st,x,initVecArray[i]);
          VecDestroy(&x);
        }
        std::cout<<"Created!"<<std::endl;
        ierr = EPSSetInitialSpace(eps,nInitVec,initVecArray);
      }
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
      ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %ld\n",its);CHKERRQ(ierr);
      ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
      ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
      ierr = KSPGetTotalIterations(ksp, &lits);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %ld\n",lits);CHKERRQ(ierr);
      ierr = EPSGetType(eps, &type);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
      ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of requested eigenvalues: %ld\n", nev);CHKERRQ(ierr);
      ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%ld\n",(double)tol,maxit);CHKERRQ(ierr);
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
    if(nconv < this->fNEigenpairs){
      this->fNEigenpairs = nconv;
    }
    w.Resize(this->fNEigenpairs);
    std::set<int> skipped;
    //let us get the eigenvalues
    {
      PetscScalar eigr{0}, eigi{0};
      int count{0};
      for(int i = 0; i < nconv; i++){
        if(count>=this->fNEigenpairs){break;}
        EPSGetEigenvalue(eps,i, &eigr, &eigi);
        CSTATE val{0};
        if constexpr(std::is_same_v<PetscScalar,STATE>){
          val = eigr + 1i * eigi;
        }else{
          val = eigr;
        }
        if(std::abs(val)>1e-8){
          w[count] = val;
          count++;
        }else{
          skipped.insert(i);
        }
      }
      w.Resize(count);
      this->fNEigenpairs = count;
    }

    if(!calcVectors) return 0;
    
    eigenVectors.Resize(pzA.Rows(),this->fNEigenpairs);
    int count{0};
    for (int i = 0; i < nconv; ++i) {
      //found everyone already
      if(count >= this->fNEigenpairs){break;}
      //small spurious ev
      if(skipped.find(i)!=skipped.end()){continue;}
      if constexpr(std::is_same_v<PetscScalar,STATE>){
        Vec eigVecRe, eigVecIm;
        PetscScalar *eigVecReArray, *eigVecImArray;
        
        MatCreateVecs(petscA,&eigVecRe,nullptr);
        MatCreateVecs(petscA,&eigVecIm,nullptr);
        EPSGetEigenvector(eps,i,eigVecRe,eigVecIm);
        VecGetArray(eigVecRe,&eigVecReArray);
        VecGetArray(eigVecIm,&eigVecImArray);
        for (int j = 0; j < pzA.Rows(); ++j) {
          eigenVectors(j,count) = eigVecReArray[j] + 1i * eigVecImArray[j];
        }
        VecRestoreArray(eigVecRe,&eigVecReArray);
        VecRestoreArray(eigVecIm,&eigVecImArray);
      }else{
        Vec eigVec;
        PetscScalar  *eigVecArray;
        MatCreateVecs(petscA,&eigVec,nullptr);
        EPSGetEigenvector(eps,i,eigVec,nullptr);
        VecGetArray(eigVec,&eigVecArray);
        for (int j = 0; j < pzA.Rows(); ++j) {
          eigenVectors(j,count) = eigVecArray[j];
        }
        VecRestoreArray(eigVec,&eigVecArray);
      }
      count++;
    }

    EPSDestroy(&eps);
    
    MatDestroy(&petscA);
    MatDestroy(&petscB);

#ifndef WGMA_PETSC_64BIT
    //otherwise we havent copied the matrices
    if(iaP) PetscFree(iaP);
    if(jaP) PetscFree(jaP);
    if(aaP) PetscFree(aaP);
    if(ibP) PetscFree(ibP);
    if(jbP) PetscFree(jbP);
    if(abP) PetscFree(abP);
#endif
    if(nInitVec){
      for(int i = 0; i < nInitVec; i++){
        VecDestroy(&initVecArray[i]);
        PetscFree(initVecMem[i]);
        initVecMem[i] = nullptr;
      }
      delete [] initVecMem;
      delete [] initVecArray;
    }
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
    TPZLinearEigenSolver<TVar>::SetAsGeneralised(isGeneralised);
    if(fProbType == EPSProblemType::EPS_NOTSET){
      std::cout<<__PRETTY_FUNCTION__
               <<"\nWARNING: call SetProblemType instead"
               <<std::endl;
      if(isGeneralised) fProbType = EPSProblemType::EPS_GNHEP;
      else fProbType = EPSProblemType::EPS_NHEP;
    }
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
  void EPSHandler<TVar>::SetPrecond(const PC pre, RTVar zero){
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
    unreachable();
  };
  
#ifdef WGMA_USING_SLEPC
  ::EPSProblemType ConvertProblemType(EPSProblemType in)
  {
    switch(in){
    case EPSProblemType::EPS_NOTSET:
      DebugStop();
    case EPSProblemType::EPS_HEP: return EPS_HEP;
    case EPSProblemType::EPS_GHEP: return EPS_GHEP;
    case EPSProblemType::EPS_NHEP: return EPS_NHEP;
    case EPSProblemType::EPS_GNHEP: return EPS_GNHEP;
    case EPSProblemType::EPS_PGNHEP: return EPS_PGNHEP;
    case EPSProblemType::EPS_GHIEP: return EPS_GHIEP;
    }
    unreachable();
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
    default:
      DebugStop();
    }
    unreachable();
  }

  ::EPSConv ConvertConv(EPSConv in)
  {
    switch(in){
    case EPSConv::EPS_CONV_ABS: return EPS_CONV_ABS;
    case EPSConv::EPS_CONV_REL: return EPS_CONV_REL;
    case EPSConv::EPS_CONV_NORM: return EPS_CONV_NORM;
    }
    unreachable();
  }
  EPSConv ConvertConv(::EPSConv in)
  {
    switch(in){
    case EPS_CONV_ABS: return EPSConv::EPS_CONV_ABS;
    case EPS_CONV_REL: return EPSConv::EPS_CONV_REL;
    case EPS_CONV_NORM: return EPSConv::EPS_CONV_NORM;
    default:
      DebugStop();
    }
    unreachable();
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
    }
    unreachable();
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

  ::EPSWhich ConvertWhich(TPZEigenSort in)
  {
    switch(in){
    case TPZEigenSort::AbsDescending: return EPS_LARGEST_MAGNITUDE;
    case TPZEigenSort::AbsAscending: return EPS_SMALLEST_MAGNITUDE;
    case TPZEigenSort::RealDescending: return EPS_LARGEST_REAL;
    case TPZEigenSort::RealAscending: return EPS_SMALLEST_REAL;
    case TPZEigenSort::ImagDescending: return EPS_LARGEST_IMAGINARY;
    case TPZEigenSort::ImagAscending: return EPS_SMALLEST_IMAGINARY;
    case TPZEigenSort::TargetMagnitude: return EPS_TARGET_MAGNITUDE;
    case TPZEigenSort::TargetRealPart: return EPS_TARGET_REAL;
    case TPZEigenSort::TargetImagPart: return EPS_TARGET_IMAGINARY;
    case TPZEigenSort::UserDefined: return EPS_WHICH_USER;
    case TPZEigenSort::Invalid:
      DebugStop();
    }
    unreachable();
  }
  
  TPZEigenSort ConvertWhich(::EPSWhich in)
  {
    switch(in){
    case EPS_LARGEST_MAGNITUDE: return TPZEigenSort::AbsDescending;
    case EPS_SMALLEST_MAGNITUDE: return TPZEigenSort::AbsAscending;
    case EPS_LARGEST_REAL: return TPZEigenSort::RealDescending;
    case EPS_SMALLEST_REAL: return TPZEigenSort::RealAscending;
    case EPS_LARGEST_IMAGINARY: return TPZEigenSort::ImagDescending;
    case EPS_SMALLEST_IMAGINARY: return TPZEigenSort::ImagAscending;
    case EPS_TARGET_MAGNITUDE: return TPZEigenSort::TargetMagnitude;
    case EPS_TARGET_REAL: return TPZEigenSort::TargetRealPart;
    case EPS_TARGET_IMAGINARY: return TPZEigenSort::TargetImagPart;
    case EPS_WHICH_USER: return TPZEigenSort::UserDefined;
    default:
      DebugStop();
    }
    unreachable();
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
    unreachable();
  }


  PC ConvertPrecond(PCType in){
    
    if(!strcmp(in,PCNONE)) return PC::NONE;
    else if(!strcmp(in,PCJACOBI)) return PC::JACOBI;
    else if(!strcmp(in,PCSOR)) return PC::SOR;
    else if(!strcmp(in,PCLU)) return PC::LU;
    else if(!strcmp(in,PCSHELL)) return PC::SHELL;
    else if(!strcmp(in,PCBJACOBI)) return PC::BJACOBI;
    else if(!strcmp(in,PCMG)) return PC::MG;
    else if(!strcmp(in,PCEISENSTAT)) return PC::EISENSTAT;
    else if(!strcmp(in,PCILU)) return PC::ILU;
    else if(!strcmp(in,PCICC)) return PC::ICC;
    else if(!strcmp(in,PCASM)) return PC::ASM;
    else if(!strcmp(in,PCGASM)) return PC::GASM;
    else if(!strcmp(in,PCKSP)) return PC::KSP;
    else if(!strcmp(in,PCREDUNDANT)) return PC::REDUNDANT;
    else if(!strcmp(in,PCCHOLESKY)) return PC::CHOLESKY;
    else if(!strcmp(in,PCPBJACOBI)) return PC::PBJACOBI;
    else if(!strcmp(in,PCVPBJACOBI)) return PC::VPBJACOBI;
    else if(!strcmp(in,PCSVD)) return PC::SVD;
    else if(!strcmp(in,PCBDDC)) return PC::BDDC;
    else if(!strcmp(in,PCKACZMARZ)) return PC::KACZMARZ;
    else if(!strcmp(in,PCTELESCOPE)) return PC::TELESCOPE;
    else if(!strcmp(in,PCPATCH)) return PC::PATCH;
    else if(!strcmp(in,PCLMVM)) return PC::LMVM;
    else if(!strcmp(in,PCHMG)) return PC::HMG;
    else if(!strcmp(in,PCDEFLATION)) return PC::DEFLATION;
    else{
      PZError<<__PRETTY_FUNCTION__
             <<"\nUnsupported type.\nAborting...\n";
      DebugStop();
    }
    return PC::NONE;
  }
  
  ::PCType ConvertPrecond(PC in){
    switch(in){
    case PC::NONE: return PCNONE;
    case PC::JACOBI: return PCJACOBI;
    case PC::SOR: return PCSOR;
    case PC::LU: return PCLU;
    case PC::SHELL: return PCSHELL;
    case PC::BJACOBI: return PCBJACOBI;
    case PC::MG: return PCMG;
    case PC::EISENSTAT: return PCEISENSTAT;
    case PC::ILU: return PCILU;
    case PC::ICC: return PCICC;
    case PC::ASM: return PCASM;
    case PC::GASM: return PCGASM;
    case PC::KSP: return PCKSP;
    case PC::REDUNDANT: return PCREDUNDANT;
    case PC::CHOLESKY: return PCCHOLESKY;
    case PC::PBJACOBI: return PCPBJACOBI;
    case PC::VPBJACOBI: return PCVPBJACOBI;
    case PC::SVD: return PCSVD;
    case PC::BDDC: return PCBDDC;
    case PC::KACZMARZ: return PCKACZMARZ;
    case PC::TELESCOPE: return PCTELESCOPE;
    case PC::PATCH: return PCPATCH;
    case PC::LMVM: return PCLMVM;
    case PC::HMG: return PCHMG;
    case PC::DEFLATION: return PCDEFLATION;
    }
    unreachable();
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
