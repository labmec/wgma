/**
slab_disc.cpp

This target performs the modal analysis of a dielectric slab
and then the subsequent scattering analysis at a waveguide discontinuity.
***/

// wgma includes
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include <TPZEigenAnalysis.h>
#include <wganalysis.hpp>
#include <scattering.hpp>
#include <util.hpp>
#include <post/solutionnorm.hpp>
#include <post/waveguideportbc.hpp>
#include <post/waveguidecoupling.hpp>
// pz includes
#include <MMeshType.h>      //for MMeshType
#include <TPZSimpleTimer.h> //for TPZSimpleTimer
#include <pzlog.h>          //for TPZLogger
#include <TPZVTKGenerator.h>
#include <pzinterpolationspace.h>

#include <regex>//for string search
#include <thread>


//!minimum shared sim data
struct SimData{
  //!.msh mesh
  std::string meshfile;
  //!.csv file containing cylinder info
  std::string cylfile;
  //!wavelength
  STATE lambda{1.5};
  //!geometric scaling (floating point precision)
  REAL scale{1};
  //!reffractive index of core
  STATE ncore{1};
  //!reffractive index of clad
  STATE nclad{1};
  //!polynomial order
  int porder{-1};
  //!pml attenuation constant in radial direction
  CSTATE alphaPMLr{0};
  //!pml attenuation constant in z-direction
  CSTATE alphaPMLz{0};
  //!whether to filter dirichlet eqs
  bool filterBoundEqs{true};
  //!renumber equations
  bool optimizeBandwidth{true};
  //!output geometric mesh in .txt and .vtk files
  bool printGmesh{false};
  //!post process fields
  bool exportVtk{false};
  //!whether to compute coupling mat
  bool couplingmat{false};
  //!vtk resolution
  int vtkRes{0};
  //!number of threads
  int nThreads{(int)std::thread::hardware_concurrency()};
  //!prefix for both meshes and output files
  std::string prefix{""};
};

//!needed data from modal analysis to create waveguide port bc
struct WgbcData{
  TPZAutoPointer<TPZCompMesh> cmesh;
  TPZFMatrix<CSTATE> wgbc_k;
  TPZVec<CSTATE> wgbc_f;
};


// Sets geometric info regarding all the cylinders in the mesh
TPZVec<wgma::gmeshtools::CylinderData> SetUpCylData(std::string_view filename,
                                                    const REAL scale);
TPZAutoPointer<wgma::wganalysis::Wgma2D>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const SimData& simdata,
  const CSTATE target,
  const int nEigenpairs,
  const TPZEigenSort sortingRule,
  bool usingSLEPC,
  const std::string &name);

void SolveScattering(TPZAutoPointer<TPZGeoMesh>gmesh,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> modal_l_an,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> modal_r_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const SimData& simdata);

SimData GetSimData()
{

  // path for output files
  const std::string path {"res_sf3d_cyl/"};
  // common prefix for both meshes and output files
  const std::string basisName{"sf3d_cyl"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  
  SimData data;
  data.meshfile="meshes/sf3d_cyl.msh";
  data.cylfile="meshes/sf3d_cyl_cyldata.csv";
  data.lambda = 4.00;
  /*
    Given the small dimensions of the domain, scaling it can help in
    achieving good precision. Using 1./k0 as a scale factor results in
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  data.scale = data.lambda/(2*M_PI);
  data.ncore = 1.4457;
  data.nclad = 1.4378;
  data.alphaPMLr = {sqrt(0.4*0.4+0.4*0.4), 0.0};
  data.alphaPMLz = {0.4, 0.0};
  data.porder = 1;
  data.filterBoundEqs = true;
  data.printGmesh=true;
  data.exportVtk = true;
  data.couplingmat = true;
  data.vtkRes=0;
  data.prefix = prefix;
  return std::move(data);
}

int main(int argc, char *argv[]) {

#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif
  const SimData simdata=GetSimData();

  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(simdata.prefix));
  
  /******************
   * eigensolver options *
   ******************/

  // how to sort eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  bool usingSLEPC {true};
  constexpr int nEigenpairs_left{50};
  constexpr int nEigenpairs_right{50};
  //if we put the exact target PETSc will complain about zero pivot in LU factor
  const CSTATE target{-1.0000001*simdata.ncore*simdata.ncore};

  /*********
   * begin *
   *********/




  /*************
   * geometry  *
   ************/
  // scoped-timer
  TPZSimpleTimer total("Total",true);

  TPZVec<std::map<std::string, int>> gmshmats;
  constexpr bool verbosity_lvl{false};
  auto gmesh = wgma::gmeshtools::ReadGmshMesh(simdata.meshfile, simdata.scale,
                                              gmshmats,
                                              verbosity_lvl);

  /**
     in order to exactly represent all the circles in the mesh,
     the python script that generated the mesh also generated this .csv file
    **/
    auto cyldata = SetUpCylData(simdata.cylfile, simdata.scale);
    wgma::gmeshtools::SetExactCylinderRepresentation(gmesh, cyldata);
    
  // print wgma_gmesh to .txt and .vtk format
  if (simdata.printGmesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = simdata.prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }

  TPZAutoPointer<wgma::wganalysis::Wgma2D> modal_l_an, modal_r_an;
  /********************************
   * cmesh(modal analysis):left   *
   ********************************/

  {
    TPZSimpleTimer timer("Modal analysis",true);
    modal_l_an = ComputeModalAnalysis(gmesh, gmshmats, simdata,
                                      target, nEigenpairs_left,
                                      sortingRule, usingSLEPC, "left");
  
    /********************************
     * cmesh(modal analysis):right   *
     ********************************/

    modal_r_an = ComputeModalAnalysis(gmesh, gmshmats, simdata,
                                      target, nEigenpairs_left,
                                      sortingRule, usingSLEPC, "right");
  }
  SolveScattering(gmesh, modal_l_an,  modal_r_an, gmshmats, simdata);
  return 0;
}

// Sets geometric info regarding all the cylindes in the mesh
TPZVec<wgma::gmeshtools::CylinderData> SetUpCylData(std::string_view filename,
                                               const REAL scale) {
  std::ifstream read(filename.data());
  if (!read) {
    std::cout << "Couldn't open the file " << filename << std::endl;
    DebugStop();
  }

  auto getNextLineAndSplitIntoTokens =
      [](std::istream &str) -> std::vector<std::string> {
    std::vector<std::string> result;
    std::string line;
    std::getline(str, line);

    std::stringstream lineStream(line);
    std::string cell;

    while (std::getline(lineStream, cell, ',')) {
      result.push_back(cell);
    }
    // This checks for a trailing comma with no data after it.
    if (!lineStream && cell.empty()) {
      // If there was a trailing comma then add an empty element.
      result.push_back("");
    }
    return result;
  };

  auto line = getNextLineAndSplitIntoTokens(read); // header
  line = getNextLineAndSplitIntoTokens(read);
  // we expect xc, yc, zc, xaxis, yaxis, zaxis, r, and matid
  TPZVec<wgma::gmeshtools::CylinderData> cyls;
  const auto factor = 1./scale;
  while (line.size() == 8) {
    wgma::gmeshtools::CylinderData cyl;
    cyl.m_xc = std::stod(line[0]) * factor;
    cyl.m_yc = std::stod(line[1]) * factor;
    cyl.m_zc = std::stod(line[2]) * factor;
    cyl.m_xaxis = std::stod(line[3]) * factor;
    cyl.m_yaxis = std::stod(line[4]) * factor;
    cyl.m_zaxis = std::stod(line[5]) * factor;
    cyl.m_radius = std::stod(line[6]) * factor;
    cyl.m_matid = std::stoi(line[7]);
    const int ncyls = cyls.size();
    cyls.Resize(ncyls + 1);
    cyls[ncyls] = cyl;

    line = getNextLineAndSplitIntoTokens(read);
  }
  return cyls;
}


void ComputeModes(wgma::wganalysis::Wgma2D &an,
                  const int nThreads);

void PostProcessModes(wgma::wganalysis::Wgma2D &an,
                      std::string filename,
                      const int vtkres,
                      std::set<int> sols = {});


TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const CSTATE target, const int nEigen,
            TPZEigenSort sorting, bool &usingSLEPC);


TPZAutoPointer<wgma::wganalysis::Wgma2D>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const SimData& simdata,
  const CSTATE target,
  const int nEigenpairs,
  const TPZEigenSort sortingRule,
  bool usingSLEPC,
  const std::string &name)
{
  auto modal_cmesh = [gmesh,&gmshmats,&simdata,name](){
    // setting up cmesh data
    const auto &nclad = simdata.nclad;
    const auto &ncore = simdata.ncore;
    const auto &alphaPMLr = simdata.alphaPMLr;
    const auto &pOrder = simdata.porder;
    const auto &lambda = simdata.lambda;
    const auto &scale = simdata.scale;
    wgma::cmeshtools::PhysicalData modal_data;

    std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
    modal_mats["src_"+name+"_clad"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad, 1.);
    modal_mats["src_"+name+ "_core"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    std::map<std::string, wgma::bc::type> modal_bcs;
    modal_bcs["src_"+name+"_bnd"] = wgma::bc::type::PEC;
    //dimension of the modal analysis 
    constexpr int modal_dim{2};
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                            {alphaPMLr,0},
                                            modal_data, modal_dim);
    //we must now filter the 2D PMLs
    std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
    for(const auto &pml : modal_data.pmlvec){
      const std::string pattern{"src_"+name+"_clad"};
      const auto rx = std::regex{pattern, std::regex_constants::icase };
    
      const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
      if(found_pattern){pmlvec.push_back(pml);}
    }
    modal_data.pmlvec = pmlvec;
    return wgma::wganalysis::CMeshWgma2D(gmesh,pOrder, modal_data,lambda, scale);
  }();

  TPZAutoPointer<wgma::wganalysis::Wgma2D>
    modal_an = new wgma::wganalysis::Wgma2D(modal_cmesh, simdata.nThreads,
                                              simdata.optimizeBandwidth,
                                              simdata.filterBoundEqs);
  {
    auto solver =
      SetupSolver(target, nEigenpairs, sortingRule, usingSLEPC);
    modal_an->SetSolver(*solver);
  }
  std::string modal_file{simdata.prefix+"_modal_"+name};
  if(usingSLEPC){
    modal_file += "_slepc";
  }else{
    modal_file += "_krylov";
  }
  ComputeModes(*modal_an, simdata.nThreads);
  if(simdata.exportVtk){
    PostProcessModes(*modal_an, modal_file, simdata.vtkRes);
  }
  return modal_an;
}


#include <TPZKrylovEigenSolver.h>
#include <TPZBndCond.h>
#include <pzstepsolver.h>
#include <slepcepshandler.hpp>
#include <precond.hpp>

void SolveWithPML(TPZAutoPointer<TPZCompMesh> scatt_cmesh,
                  wgma::wganalysis::Wgma2D& src_an,
                  const TPZVec<CSTATE> &src_coeffs,
                  const SimData &simdata);
void SetupPrecond(wgma::scattering::Analysis &scatt_an);

//utility functions
TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const CSTATE target,const int neigenpairs,
            TPZEigenSort sorting, bool &usingSLEPC)
{

#ifndef WGMA_USING_SLEPC
  if(usingSLEPC){
    std::cout<<"wgma was not configured with slepc. defaulting to: "
             <<"TPZKrylovSolver"<<std::endl;
    usingSLEPC = false;
  }
#endif

  TPZAutoPointer<TPZEigenSolver<CSTATE>> solver{nullptr};
  const int krylovDim{5*neigenpairs};
  if (usingSLEPC){
    using namespace wgma::slepc;
    /*
      The following are suggested SLEPc settings.
      NOTE: -1 stands for PETSC_DECIDE
    */
    
    constexpr STATE eps_tol = 1e-18;//PETSC_DECIDE
    constexpr int eps_max_its = -1;//PETSC_DECIDE
    constexpr EPSConv eps_conv_test = EPSConv::EPS_CONV_REL;
    constexpr EPSWhich eps_which = EPSWhich::EPS_TARGET_REAL;
    
    constexpr PC pc = PC::LU;
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
    constexpr bool eps_verbosity = true;
    
    
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


void ComputeModes(wgma::wganalysis::Wgma2D &an,
                  const int nThreads)
{
  
  TPZSimpleTimer analysis("Modal analysis");
  static constexpr bool computeVectors{true};
  
  an.Run(computeVectors);

  an.LoadAllSolutions();
  auto cmesh = an.GetMesh();
  constexpr bool conj{false};
  //leave empty for all valid matids
  std::set<int> matids {};
  wgma::post::SolutionNorm<wgma::post::MultiphysicsIntegrator>(cmesh,matids,conj,nThreads).Normalise();

  
  TPZFMatrix<CSTATE> &mesh_sol=cmesh->Solution();
  //we update analysis object
  an.SetEigenvectors(mesh_sol);
  an.LoadAllSolutions();
}
void PostProcessModes(wgma::wganalysis::Wgma2D &an,
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
  auto cmesh = an.GetMesh();
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

void SolveScattering(TPZAutoPointer<TPZGeoMesh>gmesh,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> src_an,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> match_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const SimData& simdata)
{
    
  auto CreateScattMesh = [gmesh,&gmshmats,&simdata]
    (bool extendDomains){

    const CSTATE &ncore = simdata.ncore;
    const CSTATE &nclad = simdata.nclad;
  
    const CSTATE alphaPMLr = simdata.alphaPMLr;
    const CSTATE alphaPMLz = simdata.alphaPMLz;

    const auto &pOrder = simdata.porder;
    const auto &lambda = simdata.lambda;
    const auto &scale = simdata.scale;
    
    const std::string &prefix = simdata.prefix;

    
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData scatt_data;
    std::map<std::string, std::pair<CSTATE, CSTATE>> scatt_mats;
    scatt_mats["core_left"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    scatt_mats["core_right"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    scatt_mats["cladding"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad,1.);
    std::map<std::string, wgma::bc::type> scatt_bcs;
    scatt_bcs["scatt_bnd_mid"] = wgma::bc::type::PEC;
    if(extendDomains){
      scatt_bcs["scatt_bnd_right_left"] = wgma::bc::type::PEC;
    }

    
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, scatt_mats, scatt_bcs,
                                            {alphaPMLr,0,alphaPMLz}, scatt_data);


    //materials that will represent our source
    std::set<int> src_ids;
    /*
      if the domain is extended, we will create a current source at
      source_core_left and source_clad_left.
      otherwise, the source is modelled as a waveguide port boundary
      condition and there is no need for source mats
     */
    if (extendDomains){
      const std::string srcMat[] = {"src_left_core", "src_left_clad"};
      const auto matdim = 2;
      for(const auto &mat : srcMat){
        src_ids.insert(gmshmats[matdim].at(mat));
      }
      //now we add the 1d pml mats to source mats
      for(const auto &pml : scatt_data.pmlvec){
        const std::string pattern{"src_left_clad"};
        const auto rx = std::regex{pattern, std::regex_constants::icase };
    
        const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
        if(found_pattern){
          const auto id = gmshmats[matdim].at(*pml->names.begin());
          src_ids.insert(id);
        }
      }
    }
    


    /*
      probe mats are regions of the domain in which we want to be able
      to evaluate our solution
      they are also used to ensure that the computational elements are created
      so every region that will be used in a waveguide port bc must be
      also inserted as a probe mat
    */
    std::vector<std::string> probeMats;
    {
      probeMats.push_back("src_right_clad");
      probeMats.push_back("src_right_core");
      probeMats.push_back("src_left_clad");
      probeMats.push_back("src_left_core");
      probeMats.push_back("probe_clad");
      probeMats.push_back("probe_core");
      //now for 1d pml mats
      for(const auto &pml : scatt_data.pmlvec){
        const std::string pattern_left{"src_left_clad"};
        const auto rx_left =
          std::regex{pattern_left, std::regex_constants::icase };
        const std::string right_pattern{"src_right_clad"};
        const auto rx_right =
          std::regex{right_pattern, std::regex_constants::icase };
        const std::string pattern_probe{"probe_clad"};
        const auto rx_probe =
          std::regex{pattern_probe, std::regex_constants::icase };
        
        const bool found_pattern =
          std::regex_search(*(pml->names.begin()), rx_left) ||
          std::regex_search(*(pml->names.begin()), rx_right) ||
          std::regex_search(*(pml->names.begin()), rx_probe);

        if(found_pattern){
          probeMats.push_back(*pml->names.begin());
        }
      }
      for(auto &mat : probeMats){
        const auto matdim = 2;
        const auto id = gmshmats[matdim].at(mat);
        if(src_ids.count(id)==0){
          //just to avoid inserting materials twice when extending domains
          scatt_data.probevec.push_back({id,matdim});
        }
      }
    }
    
    

    
    
    /*when using waveguide port bc we need to filter out some PML regions*/
    if (!extendDomains){
      std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
      for(const auto &pml : scatt_data.pmlvec){
        const std::string pattern_left{"zm"};
        const auto rx_left = std::regex{pattern_left, std::regex_constants::icase };
        const std::string pattern_right{"zp"};
        const auto rx_right = std::regex{pattern_right, std::regex_constants::icase };
    
        const bool found_pattern =
          std::regex_search(*(pml->names.begin()), rx_left) ||
          std::regex_search(*(pml->names.begin()), rx_right);
        
        if(found_pattern){//we discard these regions
          continue;
        }
        pmlvec.push_back(pml);
      }
      scatt_data.pmlvec = pmlvec;
    }

    return wgma::scattering::CMeshScattering3D(gmesh, pOrder, scatt_data,src_ids,
                                               lambda,scale);
  };
  /*********************
   * solve(scattering) *  
   *********************/  
  TPZSimpleTimer tscatt("Scattering");


  const auto n_eigenpairs_left = src_an->GetEigenvalues().size();
  /*
    the source is written as a linear combination of the modes
    this vector contains the coefficients of such combination
   */
  TPZVec<CSTATE> src_coeffs(n_eigenpairs_left,0);
  src_coeffs[0] = 1;

  //set up post processing
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_imag",
    "Field_abs"};

  
  auto scatt_mesh_pml = CreateScattMesh(true);
  //solve using PML as a reference solution
  {
    const std::string suffix = "pml";
    const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
    auto vtk = TPZVTKGenerator(scatt_mesh_pml, fvars, scatt_file, simdata.vtkRes);
    SolveWithPML(scatt_mesh_pml,src_an,src_coeffs,simdata);
    vtk.Do();
  }

  //let us just test the creation of the other mesh

  auto scatt_mesh_wgbc = CreateScattMesh(false);
  //solve using PML as a reference solution
  {
    const std::string suffix = "wgbc";
    const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
    auto vtk = TPZVTKGenerator(scatt_mesh_wgbc, fvars, scatt_file, simdata.vtkRes);
    vtk.Do();
  }
  
  // //now we solve varying the number of modes used in the wgbc
  
  // //index of the number of modes to be used to restrict the dofs on waveguide bcs
  // TPZVec<int> nmodes = {0,1,2,5,10,15,20,100,300};
  // src_an.LoadAllSolutions();
  // match_an.LoadAllSolutions();

  // //as an initial test, one could just simply project the solution and check
  // //the results
  // constexpr bool test_proj{false};
  // if(test_proj){
  //   ProjectSolIntoRestrictedMesh(src_an, scatt_mesh_pml, simdata, nmodes);
  // }
  
  // auto scatt_mesh_wgbc = CreateScattMesh(false);
  // const std::string suffix = "wgbc";
  // const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
  // auto vtk = TPZVTKGenerator(scatt_mesh_wgbc, fvars, scatt_file, simdata.vtkRes);

  // //compute wgbc coefficients
  // WgbcData src_data;
  // src_data.cmesh = src_an.GetMesh();
  // ComputeWgbcCoeffs(src_an, src_adjoint_sol,
  //                   src_data.wgbc_k, src_data.wgbc_f, false, src_coeffs);

  // WgbcData match_data;
  // match_data.cmesh = match_an.GetMesh();
  // ComputeWgbcCoeffs(match_an, match_adjoint_sol,
  //                   match_data.wgbc_k, match_data.wgbc_f,true, {});
    
    
  // //here we will store the error between pml approx and wgbc approx
  // /*
  // TPZAutoPointer<TPZCompMesh> error_mesh = [gmesh,&gmshmats, &simdata](){
  //   // setting up cmesh data
  //   const auto &nclad = simdata.nclad;
  //   const auto &ncore = simdata.ncore;
  //   const auto &alphaPMLx = simdata.alphaPMLx;
  //   const auto &alphaPMLy = simdata.alphaPMLy;
  //   const auto &mode = simdata.mode;
  //   const auto &pOrder = simdata.porder;
  //   const auto &lambda = simdata.lambda;
  //   const auto &scale = simdata.scale;
    
  //   wgma::cmeshtools::PhysicalData modal_data;
  //   std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
  //   modal_mats["eval_clad_left"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad, 1.);
  //   modal_mats["eval_core_left"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
  //   std::map<std::string, wgma::bc::type> modal_bcs;
  //   modal_bcs["eval_left_bnd"] = wgma::bc::type::PEC;
  //   //dimension of the modal analysis 
  //   constexpr int modal_dim{1};
  //   wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
  //                                           {alphaPMLx,alphaPMLy}, modal_data, modal_dim);
  //   //we must now filter the 1D PMLs
  //   std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
  //   for(const auto &pml : modal_data.pmlvec){
  //     const std::string pattern{"eval_clad_left"};
  //     const auto rx = std::regex{pattern, std::regex_constants::icase };
    
  //     const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
  //     if(found_pattern){pmlvec.push_back(pml);}
  //   }
  //   modal_data.pmlvec = pmlvec;
  //   return wgma::wganalysis::CMeshWgma1D(gmesh,mode,pOrder,modal_data,
  //                                        lambda, scale);
  // }();
  // */
  // TPZAutoPointer<TPZCompMesh> error_mesh = src_an.GetMesh()->Clone();
  
  // TPZFMatrix<CSTATE> sol_pml = error_mesh->Solution();
  // {
  //   Extract1DSolFrom2DMesh(error_mesh, scatt_mesh_pml, sol_pml);
  // }
  // //just to get the same size, we will zero it later
  // auto sol_wgbc = sol_pml;
  // const std::string error_file = simdata.prefix+"_error_";
  // auto vtk_error = TPZVTKGenerator(error_mesh, fvars, error_file, simdata.vtkRes);
  // std::map<int,STATE> error_res;
  // for(int im = 0; im < nmodes.size(); im++){
  //   const int nm = nmodes[im];
  //   if(!nm){continue;}
  //   RestrictDofsAndSolve(scatt_mesh_wgbc, src_data, match_data,
  //                        src_coeffs, nm,simdata);
  //   //plot
  //   vtk.Do();
  //   //now we compute the error
  //   Extract1DSolFrom2DMesh(error_mesh, scatt_mesh_wgbc, sol_wgbc);
  //   sol_wgbc -= sol_pml;
  //   error_mesh->LoadSolution(sol_wgbc);
  //   vtk_error.Do();
  //   auto normsol =
  //     wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
  //   normsol.SetNThreads(std::thread::hardware_concurrency());
  //   const auto norm = std::real(normsol.ComputeNorm()[0]);
  //   std::cout<<"nmodes "<<nm<<" error "<<norm<<std::endl;
  //   error_res.insert({nm,norm});
      
  //   //removing restrictions
  //   wgma::cmeshtools::RemovePeriodicity(scatt_mesh_wgbc);
  //   scatt_mesh_wgbc->ComputeNodElCon();
  //   scatt_mesh_wgbc->CleanUpUnconnectedNodes();
  // }
  // std::cout<<"++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  // for(auto [nm,error] : error_res){
  //   std::cout<<"nmodes "<<nm<<" norm error: "<<error<<std::endl;
  // }
  // std::cout<<"++++++++++++++++++++++++++++++++++++++++"<<std::endl;
}

void SetupPrecond(wgma::scattering::Analysis &scatt_an) {
  TPZSimpleTimer solve("SetupPrecond", true);
  constexpr REAL tol = 5e-8;
      
  auto &solver = dynamic_cast<TPZStepSolver<CSTATE>&>(scatt_an.GetSolver());

  TPZAutoPointer<TPZMatrixSolver<CSTATE>> precond;
  {
    TPZSimpleTimer pre("precond",true);
    const auto eqfilt = scatt_an.StructMatrix()->EquationFilter();
    auto scatt_cmesh = scatt_an.GetMesh();
    //we must filter out the dirichlet BCs
    std::set<int> bnd_ids;
    {
      for(auto [id,mat] : scatt_cmesh->MaterialVec()){
        auto *bnd = dynamic_cast<TPZBndCond *>(mat);
        if(bnd && bnd->Type()==0){
          bnd_ids.insert(id);
        }
      }
    }
    
    TPZVec<int64_t> eqgraph, eqgraphindex;
    wgma::precond::CreateZaglBlocks(scatt_cmesh,bnd_ids, eqfilt, eqgraph,eqgraphindex);

    TPZVec<int> colors(eqgraphindex.size()-1,0);
    auto mat = scatt_an.GetSolver().Matrix();
    const int numc =
      wgma::precond::ColorEqGraph(eqgraph,eqgraphindex,
                                  *mat,eqfilt.NActiveEquations(),colors);

    std::cout<<"created "<<eqgraphindex.size()-1
             <<" blocks split into "
             <<numc<<" colors"<<std::endl;
    TPZVec<int64_t> sparse_blocks = {0};
    precond = new wgma::precond::BlockPrecond(mat,
                                              std::move(eqgraph),
                                              std::move(eqgraphindex),
                                              colors, numc,
                                              sparse_blocks);
  }
  const int64_t n_iter = {300};
  const int n_vecs = {30};
  constexpr int64_t from_current{0};
  solver.SetGMRES(n_iter, n_vecs, *precond, tol, from_current);
}

void SolveWithPML(TPZAutoPointer<TPZCompMesh> scatt_cmesh,
                  wgma::wganalysis::Wgma2D& src_an,
                  const TPZVec<CSTATE> &src_coeffs,
                  const SimData &simdata)
{
  constexpr bool sym{false};
  auto scatt_an = wgma::scattering::Analysis(scatt_cmesh, simdata.nThreads,
                                             simdata.optimizeBandwidth,
                                             simdata.filterBoundEqs,
                                             sym);
    
  auto src_mesh = src_an.GetMesh();
  //get id of source materials
  wgma::scattering::SourceWgma src;
  for(auto [id,mat] : src_mesh->MaterialVec()){
    src.id.insert(id);
  }
  src.modal_cmesh = src_mesh;

  TPZFMatrix<CSTATE> sol(scatt_cmesh->NEquations(),1);
  //we will add the components of the solution one by one
  const int nsol = src_coeffs.size();
  bool first_assemble{true};
  for(int isol = 0; isol < nsol; isol++){
    if(IsZero(src_coeffs[isol])){continue;}
    src_an.LoadSolution(isol);
    auto beta = std::sqrt(-src_an.GetEigenvalues()[isol]);
    wgma::scattering::LoadSource2D(scatt_cmesh, src, src_coeffs[isol]);
    wgma::scattering::SetPropagationConstant(scatt_cmesh, beta);
    if(first_assemble){
      first_assemble = false;
      scatt_an.Assemble();
      SetupPrecond(scatt_an);
    }else{
      scatt_an.AssembleRhs(src.id);
    }
    scatt_an.Solve();
    scatt_an.LoadSolution();
    TPZFMatrix<CSTATE> &curr_sol = scatt_an.Solution();
    sol+=curr_sol;
  }
}