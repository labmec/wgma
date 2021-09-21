//
// Created by Francisco Teixeira Orlandini on 12/13/17.
//
#include "slepcepshandler.hpp"

#include "pzysmp.h"
#include "TPZSimpleTimer.h"

#include <slepceps.h>
#include <petsctime.h>
#include <petscksp.h>



/*******************
 *    GENERAL       *
 *******************/
namespace wgma::slepc{

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
  
  class EPSWrapper{
  public:
    EPS eps;
  };
  
  
  
  EPSHandler::EPSHandler() : fVerbose(true){
    static bool amItheFirst{true};
    if(amItheFirst){
      SlepcInitialize((int *)0, (char ***)0, (const char*)0,(const char*)0 );
      amItheFirst = false;
    }else{
      PZError<<"ERROR:\nOnly one instance of\n"
             <<"wgma::slepc::EPSHandler\n"
             <<"is allowed. Aborting...\n";
      DebugStop();
    }
    fEps = std::unique_ptr<EPSWrapper>(new EPSWrapper);
    EPSCreate( PETSC_COMM_WORLD, &(fEps->eps) );
  }
  
  EPSHandler::~EPSHandler() {
    EPSDestroy(&(fEps->eps));
    SlepcFinalize();
  }

  
  EPSHandler* EPSHandler::Clone() const {
    return (EPSHandler*)this;
  }

  
  void EPSHandler::SetNEigenpairs(int n) {
    std::cout<<__PRETTY_FUNCTION__
             <<"\nNote: setting ncv = 0 and mpd = 0.\n"
             <<"Call instead EPSHandler<T>::SetEPSDimensions";
    SetEPSDimensions(n, 0, 0);
  }
  
  
  bool EPSHandler::CheckMatrixTypes(){
    auto sparseA =
      TPZAutoPointerDynamicCast<TPZFYsmpMatrix<CSTATE>>(this->MatrixA());
    auto sparseB =
      TPZAutoPointerDynamicCast<TPZFYsmpMatrix<CSTATE>>(this->MatrixB());

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
  
  
  int EPSHandler::SolveEigenProblem(TPZVec<CSTATE> &w, TPZFMatrix<CSTATE> &ev){
    return SolveImpl(w, ev, true);
  }
  
  int EPSHandler::SolveEigenProblem(TPZVec<CSTATE> &w){
    TPZFMatrix<CSTATE> ev;
    return SolveImpl(w, ev, false);
  }

  
  int EPSHandler::SolveGeneralisedEigenProblem(TPZVec<CSTATE> &w,
                                               TPZFMatrix<CSTATE> &ev){
    return SolveImpl(w, ev, false);
  }

  
  int EPSHandler::SolveGeneralisedEigenProblem(TPZVec<CSTATE> &w){
    TPZFMatrix<CSTATE> ev;
    return SolveImpl(w, ev, false);
  }
  
  
  int EPSHandler::SolveImpl(TPZVec <CSTATE> &w, TPZFMatrix <CSTATE>& eigenVectors,
                                  bool calcVectors){
    CheckMatrixTypes();
    auto &pzA = dynamic_cast<TPZFYsmpMatrix<CSTATE>&>(this->MatrixA().operator*());
    auto &pzB = dynamic_cast<TPZFYsmpMatrix<CSTATE>&>(this->MatrixB().operator*());
    
    auto CreatePetscMat = [](TPZFYsmpMatrix<CSTATE> &pzmat,
                             Mat &mat,
                             PetscInt *&ia,
                             PetscInt *&ja,
                             PetscScalar *&aa) ->int{
      const int nRows = pzmat.Rows();
      const int nCols = pzmat.Cols();
      TPZVec<int64_t> I, J;
      TPZVec<CSTATE> A;
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
    
    PetscErrorCode ierr;
    EPSSetOperators(fEps->eps, petscA, petscB);

    {
      TPZSimpleTimer setup("EPSSetUp");
      EPSSetUp(fEps->eps);
    }

    if(fVerbose){
      EPSView(fEps->eps,PETSC_VIEWER_STDOUT_WORLD);
      PetscViewerAndFormat *vf;
      PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf);
      EPSMonitorSet(fEps->eps,
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
      ierr = EPSSolve(fEps->eps);CHKERRQ(ierr);
      ierr = PetscTime(&t2);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Elapsed Time in EPSSolve: %f\n",t2-t1);CHKERRQ(ierr);
      
    }
    /*
      Optional: Get some information from the solver and display it
    */
    PetscInt its,lits,maxit,nev;
    PetscReal tol;
    ::EPSType type;
    ::KSP ksp;
    ::ST st;
    ierr = EPSGetIterationNumber(fEps->eps, &its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
    ierr = EPSGetST(fEps->eps,&st);CHKERRQ(ierr);
    ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
    ierr = KSPGetTotalIterations(ksp, &lits);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %D\n",lits);CHKERRQ(ierr);
    ierr = EPSGetType(fEps->eps, &type);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
    ierr = EPSGetDimensions(fEps->eps,&nev,NULL,NULL);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n", nev);CHKERRQ(ierr);
    ierr = EPSGetTolerances(fEps->eps,&tol,&maxit);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Display solution and clean up
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /*
      Show detailed info unless -terse option is given by user
    */
    ::EPSConv eps_conv_test;
    EPSErrorType eps_error_type = EPS_ERROR_RELATIVE;
    EPSGetConvergenceTest(fEps->eps, &eps_conv_test);
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
      ierr = EPSConvergedReasonView(fEps->eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      ierr = EPSErrorView(fEps->eps,eps_error_type,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    } else {
      ierr = EPSErrorView(fEps->eps,eps_error_type,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
    PetscInt nconv;
    EPSGetConverged(fEps->eps, &nconv);
    const PetscInt nEigen = nev > nconv ? nconv : nev;
    w.Resize(nEigen);

    //let us sort the eigenvalues
    TPZManVector<int,20> indices;
    {
      PetscScalar eigr{0}, eigi{0};
      for(int i = 0; i < nEigen; i++){
        EPSGetEigenvalue(fEps->eps,i, &eigr, &eigi);
        //TODO: fix complex
        w[i] = eigr;
      }
      this->SortEigenvalues(w,indices);
    }

    
    eigenVectors.Resize(pzA.Rows(),nEigen);
    for (int i = 0; i < nEigen; ++i) {
      auto il = indices[i];
      //TODO: fix complex
      Vec eigVec;
      PetscScalar  *eigVecArray;
      MatCreateVecs(petscA,&eigVec,nullptr);
      EPSGetEigenvector(fEps->eps,il,eigVec,nullptr);
      VecGetArray(eigVec,&eigVecArray);
      for (int j = 0; j < pzA.Rows(); ++j) {
        eigenVectors(j,i) = eigVecArray[j];
      }
      VecRestoreArray(eigVec,&eigVecArray);
    }
    if(iaP) PetscFree(iaP);
    if(jaP) PetscFree(jaP);
    if(aaP) PetscFree(aaP);
    if(ibP) PetscFree(ibP);
    if(jbP) PetscFree(jbP);
    if(abP) PetscFree(abP);
    return 1;
  }

  
  void EPSHandler::SetProblemType(const EPSProblemType eps_problem)
  {
    ::EPSProblemType prob = ConvertProblemType(eps_problem);
    const PetscErrorCode ierr = EPSSetProblemType (fEps->eps, prob);
    if(ierr!=0){
      PZError<<"Invalid problem type.\n";
      PZError<<"Valid problems are:\n";
      PZError<<"EPS_HEP, EPS_GHEP, EPS_NHEP, EPS_GNHEP, EPS_PGNHEP or EPS_GHIEP\n";
      DebugStop();
    }
    this->SetAsGeneralised(!(prob == EPS_HEP || prob == EPS_NHEP));
  }

  
  EPSProblemType EPSHandler::GetProblemType() const
  {
    ::EPSProblemType eps_problem;
    const PetscErrorCode ierr = EPSGetProblemType (fEps->eps, &eps_problem);
    if(ierr!=0) DebugStop();
    EPSProblemType prob = ConvertProblemType(eps_problem);
    return prob;
  }

  
  void EPSHandler::SetTarget(const CSTATE target){
    this->fTarget = target;
    PetscErrorCode ierr = 1;
    ierr = EPSSetTarget (fEps->eps, target );
    ST st;
    ierr = EPSGetST(fEps->eps,&st);
    ierr = STSetType(st, STSINVERT);
    ierr = STSetShift(st, target);
    ierr = EPSSetST (fEps->eps, st);
    if(ierr!=0) DebugStop();
  }

  
  PetscScalar EPSHandler::Target() const {
    PetscScalar target;
    const PetscErrorCode ierr = EPSGetTarget (fEps->eps, &target );
    if(ierr!=0) DebugStop();
    return target;
  }

  
  void EPSHandler::SetWhichEigenpairs (const EPSWhich which){
    ::EPSWhich w = ConvertWhich(which);
    // set which portion of the eigenspectrum to solve for
    const PetscErrorCode ierr = EPSSetWhichEigenpairs (fEps->eps, w);
    if(ierr!=0) DebugStop();
  }

  
  EPSWhich EPSHandler::GetWhichEigenpairs () const {
    // set which portion of the eigenspectrum to solve for
    ::EPSWhich which;
    const PetscErrorCode ierr = EPSGetWhichEigenpairs (fEps->eps, &which);
    if(ierr!=0) DebugStop();
    EPSWhich w = ConvertWhich(which);
    return w;
  }

  
  void EPSHandler::SetKrylovOptions (const bool pLocking, const STATE restart){
    ::EPSType currentType;
    EPSGetType(fEps->eps, &currentType);
    if(strcmp(currentType,EPSKRYLOVSCHUR)){
      PZError<<"EPSType is not EPSKRYLOVSCHUR: Krylov settings will be ignored\n";
      return;
    }

    PetscBool locking = pLocking ? PETSC_TRUE : PETSC_FALSE;

    PetscErrorCode ierr = EPSKrylovSchurSetLocking (fEps->eps, locking);
    if(ierr!=0) DebugStop();
    ierr = EPSKrylovSchurSetRestart(fEps->eps, restart);
    if(ierr!=0) DebugStop();
  }

  
  void EPSHandler::GetKrylovOptions (bool &pLocking, STATE &restart) const{
    ::EPSType currentType;
    EPSGetType(fEps->eps, &currentType);
    if(strcmp(currentType,EPSKRYLOVSCHUR)){
      PZError<<"EPSType is not EPSKRYLOVSCHUR\n";
      DebugStop();
    }

    PetscBool locking;


    PetscErrorCode ierr = EPSKrylovSchurGetLocking(fEps->eps, &locking);
    if(ierr!=0) DebugStop();
    pLocking = locking;
    ierr = EPSKrylovSchurGetRestart(fEps->eps, &restart);
    if(ierr!=0) DebugStop();
  }

  
  void EPSHandler::SetTolerances(const STATE t, const int m_its) {
    auto tol = t;
    auto max_its = m_its;
    if(tol < 0) tol = PETSC_DEFAULT;
    if(max_its < 0) max_its = PETSC_DEFAULT;
    const PetscErrorCode ierr = EPSSetTolerances(fEps->eps,tol,max_its);
    if(ierr != 0) DebugStop();
  }

  
  void EPSHandler::GetTolerances(STATE &tol, int &max_its) const {
    const PetscErrorCode ierr = EPSGetTolerances(fEps->eps,&tol,&max_its);
    if(ierr != 0) DebugStop();
  }

  
  void EPSHandler::SetConvergenceTest(const EPSConv test) {
    ::EPSConv t = ConvertConv(test);
    const PetscErrorCode ierr = EPSSetConvergenceTest(fEps->eps,t);
    if(ierr != 0) DebugStop();
  }

  
  EPSConv EPSHandler::GetConvergenceTest() const {
    ::EPSConv test;
    const PetscErrorCode ierr = EPSGetConvergenceTest(fEps->eps,&test);
    if(ierr != 0) DebugStop();
    EPSConv t = ConvertConv(test);
    return t;
  }

  
  void EPSHandler::SetTrueResidual(const bool pOpt) {
    ST st;
    EPSGetST(fEps->eps, &st);
    STType type;
    STGetType(st, &type);

    if(strcmp(type,STSINVERT) && pOpt){
      PZError<<__PRETTY_FUNCTION__<<"is only available if STTYpe is STSINVERT\n";
      DebugStop();
    }
    const PetscBool opt= pOpt ? PETSC_TRUE : PETSC_FALSE;
    const PetscErrorCode ierr = EPSSetTrueResidual(fEps->eps,opt);
    if(ierr != 0) DebugStop();
  }

  
  void EPSHandler::SetType(const EPSType type) {
    ::EPSType t = ConvertType(type);
    const PetscErrorCode ierr = EPSSetType(fEps->eps,t);
    if(ierr != 0) DebugStop();
  }

  
  EPSType EPSHandler::GetType() const {
    ::EPSType type;
    const PetscErrorCode ierr = EPSGetType(fEps->eps,&type);
    if(ierr != 0) DebugStop();
    EPSType t = ConvertType(type);
    return t;
  }

  
  void EPSHandler::SetEPSDimensions(const int nev,
                                          const int nc,
                                          const int mp) {    
    this->fNEigenpairs = nev;

    const auto ncv = nc > 0 ? nc : PETSC_DEFAULT;
    const auto mpd = mp > 0 ? mp : PETSC_DEFAULT;
    const PetscErrorCode ierr = EPSSetDimensions(fEps->eps,nev,ncv,mpd);
    if(ierr != 0) DebugStop();
  }

  
  void EPSHandler::GetEPSDimensions(int &nev,
                                          int &ncv,
                                          int &mpd) const {
    const PetscErrorCode ierr = EPSGetDimensions(fEps->eps,&nev,&ncv,&mpd);
    if(ierr != 0) DebugStop();
  }

  
  void EPSHandler::SetVerbose(bool fVerbose) {
    EPSHandler::fVerbose = fVerbose;
  }

  
  void EPSHandler::SetLinearSolver(const KSPSolver solver){
    ::KSPType type = ConvertKSP(solver);
    
    PetscErrorCode ierr = 1;
    ::ST st;
    EPSGetST(fEps->eps, &st);
    ::KSP ksp;
    ierr = STGetKSP(st,&ksp);
    ierr = KSPSetType(ksp, type);
    if(ierr != 0){
      PZError<<__PRETTY_FUNCTION__
             <<"\nERROR.Aborting...\n";
      DebugStop();
    }
    STSetKSP(st,ksp);
    EPSSetST(fEps->eps,st);
  }

  
  void EPSHandler::SetLinearSolverTol(const STATE r, const STATE a,
                                      const STATE d, const int m){
    const auto rtol = r > 0 ? r : PETSC_DEFAULT;
    const auto atol = a > 0 ? a : PETSC_DEFAULT;
    const auto dtol = d > 0 ? d : PETSC_DEFAULT;
    const auto max_its = m > 0 ? m : PETSC_DEFAULT;
    
    PetscErrorCode ierr = 1;
    ::ST st;
    EPSGetST(fEps->eps, &st);
    ::KSP ksp;
    ierr = STGetKSP(st,&ksp);
    KSPSetTolerances(ksp , rtol , atol , dtol , max_its);
    if(ierr != 0){
      PZError<<__PRETTY_FUNCTION__
             <<"\nERROR.Aborting...\n";
      DebugStop();
    }
    STSetKSP(st,ksp);
    EPSSetST(fEps->eps,st);
  }

  
  void EPSHandler::GetLinearSolverTol(STATE &rtol, STATE &atol,
                                      STATE &dtol, int &max_its){
    PetscErrorCode ierr = 1;
    ::ST st;
    EPSGetST(fEps->eps, &st);
    ::KSP ksp;
    ierr = STGetKSP(st,&ksp);
    KSPGetTolerances(ksp , &rtol , &atol , &dtol , &max_its);
    if(ierr != 0){
      PZError<<__PRETTY_FUNCTION__
             <<"\nERROR.Aborting...\n";
      DebugStop();
    }
    STSetKSP(st,ksp);
    EPSSetST(fEps->eps,st);
  }
  
  
  void EPSHandler::SetPrecond(const Precond pre, STATE zero){
    PCType precond = ConvertPrecond(pre);
    PetscErrorCode ierr = 1;
    ::ST st;
    EPSGetST(fEps->eps, &st);
    ::KSP ksp;
    ierr = STGetKSP(st,&ksp);
    ::PC pc;
    ierr = KSPGetPC(ksp, &pc);
    ierr = PCSetType(pc, precond);
    ierr = PCFactorSetZeroPivot(pc,zero);
#ifdef PETSC_HAVE_MUMPS
    if(!strcmp(precond,"lu")||!strcmp(precond,"cholesky")){
      ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
    }
#endif
    ierr = KSPSetPC(ksp,pc);
    ierr = STSetKSP(st,ksp);
    ierr = EPSSetST(fEps->eps,st);
    if(ierr != 0) DebugStop();
  }

  ::EPSProblemType ConvertProblemType(EPSProblemType in)
  {
    switch(in){
    case EPSProblemType::EPS_HEP: return EPS_HEP;
    case EPSProblemType::EPS_GHEP: return EPS_GHEP;
    case EPSProblemType::EPS_NHEP: return EPS_NHEP;
    case EPSProblemType::EPS_GNHEP: return EPS_GNHEP;
    case EPSProblemType::EPS_PGNHEP: return EPS_PGNHEP;
    case EPSProblemType::EPS_GHIEP: return EPS_GHIEP;
    }
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
  
}
