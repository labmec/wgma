/**
slab_disc.cpp

This target performs the modal analysis of a dielectric slab
and then the subsequent scattering analysis at a waveguide discontinuity.
***/

// wgma includes
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include <wganalysis.hpp>
#include <scattering.hpp>
#include <util.hpp>
#include <slepcepshandler.hpp>
#include <post/orthosol.hpp>
// pz includes
#include <MMeshType.h>      //for MMeshType
#include <TPZSimpleTimer.h> //for TPZSimpleTimer
#include <pzlog.h>          //for TPZLogger
#include <TPZVTKGenerator.h>
#include <TPZVTKGeoMesh.h>
#include <pzbuildmultiphysicsmesh.h>

#include <regex>//for string search



#include <pzinterpolationspace.h>

constexpr bool optimizeBandwidth{false};
constexpr bool filterBoundaryEqs{true};
constexpr int nThreads{8};
constexpr int nThreadsDebug{0};
constexpr int vtkRes{0};

void CheckNodes(TPZAutoPointer<TPZGeoMesh> gmesh,
                std::map<int64_t,int64_t> &periodic_els);

void CompareDofs(TPZAutoPointer<TPZCompMesh> cmesh,
                 std::map<int64_t,int64_t> &periodic_els);

TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const CSTATE target, const int neigenpairs,
            TPZEigenSort sorting, bool usingSLEPC);

void ComputeModes(wgma::wganalysis::Wgma &an,
                  TPZAutoPointer<TPZCompMesh> mesh,
                  const bool orthogonalise,
                  const int nThreads);

void PostProcessModes(wgma::wganalysis::Wgma2D &an,
                      TPZAutoPointer<TPZCompMesh> cmesh,
                      std::string filename,
                      const int vtkres,
                      std::set<int> sols = {});



void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          wgma::wganalysis::Wgma& src_an,
                          wgma::wganalysis::Wgma& match_an,
                          const std::vector<int> &sources,
                          const std::set<int> &source_mats,
                          const std::vector<int> &nmodes,
                          const TPZVec<std::map<std::string, int>> &gmshmats,
                          const std::string &prefix);

int main(int argc, char *argv[]) {
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif
  /***********************
   * setting the problem *
   ***********************/
  // operational wavelength

  // the meshes were designed in micrometers, so lambda has to follow
  constexpr STATE lambda{1.5};

  
  constexpr STATE ncore{1.5};
  constexpr STATE nclad{1.000};

  constexpr STATE alphaPMLx {1.0};
  constexpr STATE alphaPMLy {1.0};
  constexpr STATE alphaPMLz {1.5};
  /*
    Given the small dimensions of the domain, scaling it can help in
    achieving good precision. Using 1./k0 as a scale factor results in
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  constexpr REAL scale{lambda / (2 * M_PI)};

  /******************
   *  fem options   *
   ******************/
  // polynomial order to be used in the approximation
  constexpr int pOrder{2};

  /******************
   * solver options *
   ******************/

  // how to sort eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  constexpr bool usingSLEPC {false};

  /*********************
   * exporting options *
   *********************/

  // whether to print the geometric mesh in .txt and .vtk formats
  constexpr bool printGMesh{true};
  // whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  // path for output files
  const std::string path {"res_slab_3d/"};
  // common prefix for both meshes and output files
  const std::string basisName{"slab_3d"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(prefix));


  constexpr int nEigenpairs_left{10};
  constexpr int nEigenpairs_right{50};

  constexpr CSTATE target{-ncore*ncore};

  /*********
   * begin *
   *********/




  /*************
   * geometry  *
   ************/
  // scoped-timer
  TPZSimpleTimer total("Total");

  // creates gmesh
  // file containing the .msh mesh
  const std::string meshfile{"meshes/"+basisName+".msh"};

  TPZVec<std::map<std::string, int>> gmshmats;
  constexpr bool verbosity_lvl{true};
  std::map<int64_t,int64_t> periodic_els;
  auto gmesh = wgma::gmeshtools::ReadPeriodicGmshMesh(meshfile, scale,gmshmats,
                                                      periodic_els,verbosity_lvl);

  CheckNodes(gmesh, periodic_els);
  // auto gmesh = wgma::gmeshtools::ReadGmshMesh(meshfile, scale,gmshmats,verbosity_lvl);

  // print wgma_gmesh to .txt and .vtk format
  if (printGMesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }


  constexpr bool true_periodic{true};
  /********************************
   * cmesh(modal analysis: left)  *
   ********************************/

  auto modal_l_cmesh = [gmesh,&gmshmats,&periodic_els, scale](){
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData modal_data;

    std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
    modal_mats["source_clad_left"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad, 1.);
    modal_mats["source_core_left"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    std::map<std::string, wgma::bc::type> modal_bcs;
    modal_bcs["source_left_bnd_pec"] = wgma::bc::type::PEC;
    modal_bcs["source_left_bnd_per_dep"] = wgma::bc::type::PERIODIC;
    modal_bcs["source_left_bnd_per_indep"] = wgma::bc::type::PERIODIC;
    //dimension of the modal analysis 
    constexpr int modal_dim{2};
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs, alphaPMLx,
                                            alphaPMLy, modal_data, modal_dim);
    //we must now filter the 2D PMLs
    std::vector<wgma::pml::data>  pmlvec;
    for(const auto &pml : modal_data.pmlvec){
      const std::string pattern{"source_[a-z]+_left"};
      const auto rx = std::regex{pattern, std::regex_constants::icase };
    
      const bool found_pattern = std::regex_search(*(pml.names.begin()), rx);
      if(found_pattern){pmlvec.push_back(pml);}
    }
    modal_data.pmlvec = pmlvec;
    constexpr bool verbose{true};
    if(true_periodic){
      return wgma::wganalysis::CMeshWgma2DPeriodic(gmesh, pOrder, modal_data, periodic_els, lambda, scale, verbose);
    }else{
      return wgma::wganalysis::CMeshWgma2D(gmesh, pOrder, modal_data, lambda, scale, verbose);
    }
  }();


  {
    std::string modal_left_file{prefix+"_modal_left_cmsh.vtk"};
    std::ofstream file(modal_left_file);
    TPZVTKGeoMesh::PrintCMeshVTK(modal_l_cmesh[0].operator->(), file, true);
  }
  /******************************
   * solve(modal analysis left) *
   ******************************/
  auto solver_left = SetupSolver(target, nEigenpairs_left, sortingRule, usingSLEPC);

  
  wgma::wganalysis::Wgma2D
    modal_l_an(modal_l_cmesh, nThreads, optimizeBandwidth, filterBoundaryEqs);
  modal_l_an.SetSolver(*solver_left);


  std::string modal_left_file{prefix+"_modal_left"};
  //no need to orthogonalise modes
  constexpr bool ortho{false};
  ComputeModes(modal_l_an, modal_l_cmesh[0], ortho, nThreadsDebug);

  PostProcessModes(modal_l_an, modal_l_cmesh[0], modal_left_file+"_ref", vtkRes);

  if(true_periodic){
    wgma::cmeshtools::RemovePeriodicity(modal_l_cmesh[1]);
    wgma::cmeshtools::RemovePeriodicity(modal_l_cmesh[2]);
    return 0;
  }
  
  const int isol = 0;
  //load first eigenvector
  modal_l_an.LoadSolution(isol);

  const auto nrows = modal_l_an.GetEigenvectors().Rows();
  
  TPZFMatrix<CSTATE> sol(nrows,1);
  modal_l_an.GetEigenvectors().GetSub(0, isol, nrows, 1, sol);
  //beta squared
  auto beta2 = modal_l_an.GetEigenvalues()[isol];


  auto CreateVTK = [prefix](TPZAutoPointer<TPZCompMesh> cmesh, std::string filename){
    TPZVec<std::string> fvars = {
      "Ez_real",
      "Ez_abs",
      "Et_real",
      "Et_abs"};
    return TPZVTKGenerator(cmesh, fvars, filename, vtkRes);
  };
  const std::string file_left{prefix+"_modal_left"};
  auto vtk_left = CreateVTK(modal_l_cmesh[0], file_left);
  vtk_left.Do();

  TPZManVector<TPZAutoPointer<TPZCompMesh>,2> meshvec = {modal_l_cmesh[1],
    modal_l_cmesh[2]};
  if(!true_periodic){
    std::cout<<"hcurl mesh"<<std::endl;
    CompareDofs(modal_l_cmesh[1], periodic_els);
    std::cout<<"h1 mesh"<<std::endl;
    CompareDofs(modal_l_cmesh[2],periodic_els);
  }
  
  auto CalcResidual = [](wgma::wganalysis::Wgma2D &ma,
                         TPZVec<TPZAutoPointer<TPZCompMesh>>meshvec,
                         CSTATE lambda, TPZVTKGenerator &vtk){
    ma.GetSolver().ResetMatrix();
    ma.Assemble(TPZEigenAnalysis::Mat::A);
    ma.Assemble(TPZEigenAnalysis::Mat::B);
    auto matA = ma.GetSolver().MatrixA();
    auto matB = ma.GetSolver().MatrixB();
    TPZFMatrix<CSTATE>  xfull = ma.GetMesh()->Solution();
    //if there are dependent connects, we need to eliminate these equatios
    xfull.Resize(ma.GetMesh()->NEquations(),1);
    std::cout<<" A ("<<matA->Rows()<<","<<matA->Cols()<<")\n";
    std::cout<<" B ("<<matB->Rows()<<","<<matB->Cols()<<")\n";
    std::cout<<" xfull ("<<xfull.Rows()<<","<<xfull.Cols()<<")\n";    
    TPZFMatrix<CSTATE> x(matA->Rows(),1);
    ma.StructMatrix()->EquationFilter().Gather(xfull, x);
    
    std::cout<<" x ("<<x.Rows()<<","<<x.Cols()<<")\n";
    auto res = (*matA) * x - (*matB) * x * lambda;

    const auto nrows = res.Rows();
    STATE resnorm = 0;
    for(auto ir  =0; ir < nrows; ir++){
      resnorm += std::real(
        res.GetVal(ir,0)*std::conj(res.GetVal(ir,0)
                                   ));
        }

    TPZFMatrix<CSTATE> resfull(xfull.Rows(),1);
    ma.StructMatrix()->EquationFilter().Scatter(res, resfull);
    ma.GetMesh()->LoadSolution(resfull);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec,ma.GetMesh());
    vtk.Do();
    ma.GetMesh()->LoadSolution(xfull);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec,ma.GetMesh());
    return resnorm;
    
  };

  const auto resbefore = CalcResidual(modal_l_an, meshvec, beta2, vtk_left);
  
  std::cout<<"Residual (before): "<<resbefore<<std::endl;
  //now we need to apply the periodicity to each atomic mesh
  auto SetPeriodic = [&periodic_els](TPZAutoPointer<TPZCompMesh> cmesh){
      auto gmesh = cmesh->Reference();
      gmesh->ResetReference();
      cmesh->LoadReferences();
      //let us copy the connects
      for(auto [dep, indep] : periodic_els){
        //geometric elements
        auto *dep_gel = gmesh->Element(dep);
        const auto *indep_gel = gmesh->Element(indep);
        //computational element
        auto *indep_cel = indep_gel->Reference();
        auto *dep_cel = dep_gel->Reference();
        //number of connects
        const auto n_dep_con = dep_cel->NConnects();
        const auto n_indep_con = indep_cel->NConnects();

        //now we create dependencies between connects
        for(auto ic = 0; ic < n_indep_con; ic++){
          const auto indep_ci = indep_cel->ConnectIndex(ic);
          const auto dep_ci = dep_cel->ConnectIndex(ic);

          auto &dep_con = dep_cel->Connect(ic);
          const auto ndof = dep_con.NDof(cmesh);
          if(ndof==0) {continue;}
          constexpr int64_t ipos{0};
          constexpr int64_t jpos{0};
      
          TPZFMatrix<REAL> mat(ndof,ndof);
          mat.Identity();
          dep_con.AddDependency(dep_ci, indep_ci, mat, ipos,jpos,ndof,ndof);
        } 
      }
      cmesh->CleanUpUnconnectedNodes();
    };

  SetPeriodic(modal_l_cmesh[0]);
  SetPeriodic(modal_l_cmesh[1]);
  SetPeriodic(modal_l_cmesh[2]);


  wgma::wganalysis::Wgma2D
    periodic_an(modal_l_cmesh, nThreads, optimizeBandwidth, filterBoundaryEqs);
  periodic_an.SetSolver(*solver_left);


  const std::string file_periodic{prefix+"_modal_left_periodic"};
  auto vtk_periodic = CreateVTK(modal_l_cmesh[0], file_periodic);
  TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec,modal_l_cmesh[0]);
  vtk_periodic.Do();

  const auto resafter = CalcResidual(periodic_an,meshvec, beta2, vtk_periodic);
  std::cout<<"Residual (after): "<<resafter<<std::endl;
 


  wgma::cmeshtools::RemovePeriodicity(modal_l_cmesh[1]);
  wgma::cmeshtools::RemovePeriodicity(modal_l_cmesh[2]);
  return 0;
}

void CheckNodes(TPZAutoPointer<TPZGeoMesh> gmesh,
                std::map<int64_t,int64_t> &periodic_els)
{
  for(auto [dep, indep] : periodic_els){
    auto depgel = gmesh->Element(dep);
    auto indepgel = gmesh->Element(indep);
    constexpr int nnodes = 2;
    if(depgel->NNodes() != 2 || indepgel->NNodes() != 2){
      PZError<<"WRONG NUMBER OF NODES"<<std::endl;
      DebugStop();
    }
    const bool dep_orient = depgel->Node(0).Id() > depgel->Node(1).Id();
    const bool indep_orient = indepgel->Node(0).Id() > indepgel->Node(1).Id();
    if(dep_orient != indep_orient){
      std::cout<<"dep gel "<<depgel->Index()
               <<" indep gel "<<indepgel->Index()
               <<std::endl;
      std::cout<<"dep gel nodes:";
      for(auto in = 0; in < nnodes; in++){
        std::cout<<'\t'<<depgel->Node(in).Id();
      }
      std::cout<<std::endl;
      std::cout<<"indep gel nodes:";
      for(auto in = 0; in < nnodes; in++){
        std::cout<<'\t'<<indepgel->Node(in).Id();
      }
      DebugStop(); 
    }
    
  }
}

void CompareDofs(TPZAutoPointer<TPZCompMesh> cmesh,
                 std::map<int64_t,int64_t> &periodic_els)
{
  auto gmesh = cmesh->Reference();
  ///set reference
  gmesh->ResetReference();
  cmesh->LoadReferences();
  const TPZFMatrix<CSTATE> &meshsol = cmesh->Solution();
  TPZBlock &block = cmesh->Block();
  for(auto [dep, indep] : periodic_els){
    auto depgel = gmesh->Element(dep);
    auto indepgel = gmesh->Element(indep);
    auto depcel = dynamic_cast<TPZInterpolationSpace*>(depgel->Reference());
    auto indepcel = dynamic_cast<TPZInterpolationSpace*>(indepgel->Reference());

    

    auto PrintDofs = [&meshsol, &block] (TPZInterpolationSpace *cel1,
                                         TPZInterpolationSpace *cel2){
      const auto ncon = cel1->NConnects();
      
      for (auto icon = 0; icon < ncon; icon++){
        const auto &con1 = cel1->Connect(icon);
        const auto seqnum1 = con1.SequenceNumber();
        const auto blpos1 = block.Position(seqnum1);

        const auto &con2 = cel2->Connect(icon);
        const auto seqnum2 = con2.SequenceNumber();
        const auto blpos2 = block.Position(seqnum2);
        
        const auto nvar = block.Size(seqnum1);
        for(auto i = 0; i < nvar; i++){
          auto sol1 = meshsol.GetVal(blpos1 + i, 0);
          auto sol2 = meshsol.GetVal(blpos2 + i, 0);
          auto diff = std::abs(sol1-sol2);
          std::cout<<'\t'
                   <<sol1<<'\t'
                   <<sol2<<'\t'
                   <<diff<<std::endl;
        }
        std::cout<<std::endl;
      }

      TPZMaterialDataT<CSTATE> data1, data2;
      constexpr int npts{5};
      if(cel1->Reference()->Node(0).Id() == 299){
        std::cout<<"NOWWWWWWWWWW"<<std::endl;
      }
      cel1->InitMaterialData(data1);
      cel2->InitMaterialData(data2);
      for(int i = 0;i < npts; i++){
        REAL xi = -1 + ((REAL)i)/(npts-1)*2;
        TPZManVector<REAL,1> xipt = {xi};
        TPZManVector<REAL,3> x1(3,0.), x2(3,0.);
        cel1->Reference()->X(xipt, x1);
        cel2->Reference()->X(xipt, x2);
        cel1->ComputeShape(xipt, data1);
        cel2->ComputeShape(xipt, data2);
        std::cout<<"\txi: "<<xi<<" x1: "<<x1[1]<<" x2: "<<x2[1]<<std::endl;
        const int nshape = data1.phi.Rows();
        for(int is = 0; is < nshape; is++){
          std::cout<<"\t\tshape "<<is<<": "
                   <<data1.phi.GetVal(is,0)<<" "<<data2.phi.GetVal(is,0)<<'\n';
        }
      }
      
    };

    std::cout<<"dep gel "<<depgel->Index()<<" indep gel "<<indepgel->Index()<<std::endl;
    const auto nnodes = depgel->NNodes();
    std::cout<<"dep gel nodes:\n";
    TPZManVector<REAL,3> x(3,0.);
    for(auto in = 0; in < nnodes; in++){
      std::cout<<'\t';
      depgel->Node(in).Print(std::cout);
    }
    std::cout<<std::endl;
    std::cout<<"indep gel nodes:\n";
    for(auto in = 0; in < nnodes; in++){
      std::cout<<'\t';
      indepgel->Node(in).Print(std::cout);
    }
    std::cout<<std::endl;
    PrintDofs(depcel, indepcel);
  }
}

#include <slepcepshandler.hpp>
#include <TPZKrylovEigenSolver.h>
//utility functions
TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const CSTATE target,const int neigenpairs,
            TPZEigenSort sorting, bool usingSLEPC)
{

#ifndef WGMA_USING_SLEPC
  if(usingSLEPC){
    std::cout<<"wgma was not configured with slepc. defaulting to: "
             <<"TPZKrylovSolver"<<std::endl;
    usingSLEPC = false;
  }
#endif

  TPZAutoPointer<TPZEigenSolver<CSTATE>> solver{nullptr};
  const int krylovDim = std::max(5,2*neigenpairs);
  if (usingSLEPC){
    using namespace wgma::slepc;
    /*
      The following are suggested SLEPc settings.
      NOTE: -1 stands for PETSC_DECIDE
    */
    
    constexpr STATE eps_tol = -1;//PETSC_DECIDE
    constexpr int eps_max_its = -1;//PETSC_DECIDE
    constexpr EPSConv eps_conv_test = EPSConv::EPS_CONV_REL;
    constexpr EPSWhich eps_which = EPSWhich::EPS_TARGET_REAL;
    
    constexpr Precond pc = Precond::LU;
    constexpr KSPSolver linsolver = KSPSolver::PREONLY;
    constexpr STATE ksp_rtol = -1;//PETSC_DECIDE
    constexpr STATE ksp_atol = -1;//PETSC_DECIDE
    constexpr STATE ksp_dtol = -1;//PETSC_DECIDE
    constexpr STATE ksp_max_its = -1;//PETSC_DECIDE
    constexpr bool eps_true_residual = false;
    constexpr EPSProblemType eps_prob_type = EPSProblemType::EPS_GNHEP;//do NOT change
    constexpr EPSType eps_solver_type = EPSType::KRYLOVSCHUR;
    constexpr bool eps_krylov_locking = true;
    constexpr STATE eps_krylov_restart = 0.7;
    constexpr STATE eps_mpd = -1;//PETSC_DECIDE
    constexpr bool eps_verbosity = false;
    
    
    auto eps_solver = new EPSHandler<CSTATE>;
    eps_solver->SetType(eps_solver_type);
    eps_solver->SetProblemType(eps_prob_type);
    eps_solver->SetEPSDimensions(neigenpairs, krylovDim, eps_mpd);
    eps_solver->SetWhichEigenpairs(eps_which);
    eps_solver->SetTarget(target);
    eps_solver->SetTolerances(eps_tol,eps_max_its);
    eps_solver->SetConvergenceTest(eps_conv_test);
    eps_solver->SetKrylovOptions(eps_krylov_locking,eps_krylov_restart);
    eps_solver->SetVerbose(eps_verbosity);
    eps_solver->SetTrueResidual(eps_true_residual);
    
    eps_solver->SetLinearSolver(linsolver);
    eps_solver->SetLinearSolverTol(ksp_rtol,ksp_atol,ksp_dtol,ksp_max_its);
    eps_solver->SetPrecond(pc, 1e-16);

    solver = eps_solver;
  }else{
    auto krylov_solver = new TPZKrylovEigenSolver<CSTATE>;
    TPZSTShiftAndInvert<CSTATE> st;
    krylov_solver->SetSpectralTransform(st);
    krylov_solver->SetTarget(target);
    krylov_solver->SetKrylovDim(krylovDim);
    krylov_solver->SetNEigenpairs(neigenpairs);
    krylov_solver->SetAsGeneralised(true);
    
    solver = krylov_solver;
  }
  
  solver->SetEigenSorting(sorting);
  return solver;
}

void ComputeModes(wgma::wganalysis::Wgma &an,
                  TPZAutoPointer<TPZCompMesh> cmesh,
                  const bool orthogonalise,
                  const int nThreads)
{
  
  TPZSimpleTimer analysis("Modal analysis");
  static constexpr bool computeVectors{true};
  an.Run(computeVectors);
  //let us discard all solutions with negative real part of beta

  auto &ev = an.GetEigenvalues();
  auto &evec = an.GetEigenvectors();

  //load all obtained modes into the mesh
  an.LoadAllSolutions();
  if(orthogonalise){
    //leave empty for all valid matids
    std::set<int> matids {};
    auto ortho = wgma::post::OrthoSol(cmesh, matids, nThreads);
    //orthogonalise the modes
    auto normsol = ortho.Orthogonalise();
    //let us set the orthogonalised modes
    an.SetEigenvectors(normsol);
  }

}

void PostProcessModes(wgma::wganalysis::Wgma2D &an,
                      TPZAutoPointer<TPZCompMesh> cmesh,
                      std::string filename,
                      const int vtkres,
                      std::set<int> sols){
  TPZSimpleTimer postProc("Post processing");
    
  const std::string file = filename;
  TPZVec<std::string> fvars = {
      "Ez_real",
      "Ez_abs",
      "Et_real",
      "Et_abs"};
  auto vtk = TPZVTKGenerator(cmesh, fvars, file, vtkres);

  if(sols.size() == 0){
    const auto nsol = an.GetEigenvectors().Cols();
    for(auto is = 0; is < nsol; is++){
      sols.insert(is);
    }
  }
  
  std::cout<<"Exporting "<<sols.size()<<" solutions"<<std::endl;
  for(auto isol : sols){
    an.LoadSolution(isol);
    vtk.Do();
  }
}

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          wgma::wganalysis::Wgma& src_an,
                          wgma::wganalysis::Wgma& match_an,
                          const std::vector<int> &sources,
                          const std::set<int> &source_mats,
                          const std::vector<int> &nmodes,
                          const TPZVec<std::map<std::string, int>> &gmshmats,
                          const std::string &prefix)
{
  auto ev = match_an.GetEigenvectors();
  const auto evr = ev.Rows();
  const auto evc = ev.Cols();

  auto gmesh = scatt_mesh->Reference();
  
  auto match_mesh = match_an.GetMesh();

  //set up post processing
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_imag",
    "Field_abs"};
  const std::string scatt_file = prefix+"_scatt";
  auto vtk = TPZVTKGenerator(scatt_mesh, fvars, scatt_file, vtkRes);
  
  for(auto nm : nmodes){

    if(nm){
      if(nm>evc){
        std::cout<<"nm "<<nm<<" bigger than n of computed modes "
                 <<evc<<" skipping.."<<std::endl;
        continue;
      }
      const TPZFMatrix<CSTATE> new_ev(evr,nm,ev.Elem(),evr*evc);
      match_an.SetEigenvectors(new_ev);
      //this "dummy" connect will impose the restriction over the 1D line
      auto indep_con_id = scatt_mesh->AllocateNewConnect(nm,1,1);
      auto &indep_con = scatt_mesh->ConnectVec()[indep_con_id];
      indep_con.IncrementElConnected();
      
      /*we will iterate over the elements of the modal analysis mesh
        and find their counterpart on the scattering mesh*/
      const auto nshape = indep_con.NShape();
      //we load all the modes in the comp mesh
      match_an.LoadAllSolutions();
      //the geometric mesh will point to the scattering mesh
      gmesh->ResetReference();
      scatt_mesh->LoadReferences();

      //so we can access the dofs
      TPZFMatrix<CSTATE> &modal_sol = match_mesh->Solution();
      auto &modal_block = match_mesh->Block();
      for(auto modal_el : match_mesh->ElementVec()){
        if(modal_el->Dimension()< match_mesh->Dimension()){continue;}
        auto scatt_el = modal_el->Reference()->Reference();
        if(!scatt_el){continue;}
        const auto ncon = scatt_el->NConnects();
        if(ncon != modal_el->NConnects()){
          DebugStop();
        }
        for(auto icon = 0; icon < ncon; icon++){
          auto &scatt_con = scatt_el->Connect(icon);
          auto &modal_con = modal_el->Connect(icon);
        
          const int64_t modal_dfseq = modal_con.SequenceNumber();
          const int64_t modal_pos = modal_block.Position(modal_dfseq);
          const int modal_dfvar = modal_block.Size(modal_dfseq);
          scatt_con.AddDependency(scatt_el->ConnectIndex(icon), indep_con_id, modal_sol, modal_pos, 0, modal_dfvar, nshape);
        }
      }
      scatt_mesh->ComputeNodElCon();
      scatt_mesh->CleanUpUnconnectedNodes();
    }

    auto scatt_an = wgma::scattering::Analysis(scatt_mesh, nThreads,
                                               optimizeBandwidth, filterBoundaryEqs);

    //in the first run we assemble the whole algebraic system. afterwards we can
    //just compute the rhs
    bool firstrun=true;

    /*********************
     * cmesh(scattering) *
     *********************/

    auto src_mesh = src_an.GetMesh();
    //get id of source materials
    wgma::scattering::SourceWgma src;
    src.id = source_mats;
    src.modal_cmesh = src_mesh;
    const int nsols = sources.size();
  
    for(int isol = 0; isol < nsols; isol++){
      std::cout<<"running source "<<isol+1<<" out of "<<nsols<<std::endl;
      src_an.LoadSolution(sources[isol]);
      auto beta = std::sqrt(src_an.GetEigenvalues()[isol]);
      wgma::scattering::LoadSource(scatt_mesh, src);
      wgma::scattering::SetPropagationConstant(scatt_mesh, beta);
      
      if(firstrun){
        firstrun=false;
        scatt_an.Assemble();
      }else{
        scatt_an.AssembleRhs(source_mats);
      }
      scatt_an.Solve();
      vtk.Do();
    }

    wgma::cmeshtools::RemovePeriodicity(scatt_mesh);
    scatt_mesh->ComputeNodElCon();
    scatt_mesh->CleanUpUnconnectedNodes();
  }

  match_an.SetEigenvectors(ev);
}