/**
sf3d_cyl.cpp

This target is used to validate the waveguide port boundary condition
using a slab waveguide as an example.
It consists of a "shrinked" version on slab_disc.cpp.
***/


#include "post/orthowgsol.hpp"     // for OrthoWgSol
#include "post/wgnorm.hpp"         // for WgNorm
#include <post/solutionnorm.hpp>   // for SolutionNorm
#include <post/waveguideportbc.hpp>// for WaveguidePortBC
#include <post/waveguidecoupling.hpp>//for WaveguideCoupling
#include <materials/solutionprojection.hpp>
#include "wganalysis.hpp"          // for Wgma2D, CMeshWgma2D2D
#include "scattering.hpp"          // for CMeshScattering3D
#include "bctypes.hpp"             // for type, type::PEC
#include "cmeshtools.hpp"          // for PhysicalData, SetupGmshMaterialData
#include "gmeshtools.hpp"          // for PrintGeoMesh, ReadPeriodicGmshMesh
#include "pmltypes.hpp"            // for data
#include "post/integrator.hpp"     // for SingleSpaceIntegrator
#include "slepcepshandler.hpp"     // for EPSHandler, EPSConv, EPSProblemType
#include "util.hpp"                // for CreatePath, ExtractPath
#include <json_util.hpp>

#include <TPZKrylovEigenSolver.h>  // for TPZKrylovEigenSolver
#include <TPZMatrixSolver.h>       // for TPZPardisoSolver
#include <TPZSimpleTimer.h>        // for TPZSimpleTimer
#include <TPZVTKGenerator.h>       // for TPZVTKGenerator
#include <TPZParallelUtils.h>
#include <pzbuildmultiphysicsmesh.h>
#include <Electromagnetics/TPZWgma.h>
#include <precond.hpp>
#include <pzstepsolver.h>
#include <pzcmesh.h>               // for TPZCompMesh
#include <pzerror.h>               // for PZError
#include <pzfmatrix.h>             // for TPZFMatrix
#include <pzgmesh.h>               // for TPZGeoMesh
#include <pzlog.h>                 // for TPZLogger
#include <TPZPardisoSolver.h>
#include <TPZNullMaterial.h>
#include <pzintel.h>
#include <pzelementgroup.h>
#include <pzvec_extras.h>
#include <regex>                   // for regex_search, match_results<>::_Un...


using namespace std::complex_literals;
//!minimum shared sim data
struct SimData{
  //!.msh mesh
  std::string meshfile;
  //!.csv file containing cylinder info
  std::string cylfile;
  //!wavelength
  STATE lambda{4.0};
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
  //!number of eigenvalues computed on left port
  int n_eigenpairs_left;
  //!number of eigenvalues computed on right port
  int n_eigenpairs_right;
  //!perform validation test on waveguide cross section (no discontinuity)
  bool check_mode_propagation;
  //!if set to true, wpbc is compared against PML in the mode propagation test
  bool compare_pml;
  //! whether to use direct solver
  bool direct_solver;
  //!pairs of mode index/coefficient to be used as source
  std::vector<std::pair<int,double>> source_coeffs;
  //!index of the number of modes to be used to restrict the dofs on the wpbc(in)
  TPZVec<int> n_modes_left;
  //!index of the number of modes to be used to restrict the dofs on the wpbc(out)
  TPZVec<int> n_modes_right;
  //!whether to compute coupling mat
  bool couplingmat{false};
  //!whether to filter dirichlet eqs
  bool filter_bnd_eqs{true};
  //!renumber equations (for modal analysis, scattering is always false)
  bool optimize_bandwidth{true};
  //!output geometric mesh in .txt and .vtk files
  bool print_gmesh{false};
  //!post process modal fields
  bool export_vtk_modes{false};
  //!export .csv with computed eigenvalues
  bool export_csv_modes{false};
  //!export csv file of wpbc comparison
  bool export_csv_error{false};
  //!post process scatt fields
  bool export_vtk_scatt{false};
  //!post process error of scatt field at waveguide port
  bool export_vtk_error{false};
  
  //!vtk resolution
  int vtk_res{0};
  //!number of threads
  int n_threads{(int)std::thread::hardware_concurrency()};
  //!prefix for both meshes and output files
  std::string prefix{""};
};

//! Reads sim data from file
SimData ReadSimData(const std::string &dataname);

//!needed data from modal analysis to create waveguide port bc
struct WpbcData{
  TPZAutoPointer<TPZCompMesh> cmesh;
  TPZFMatrix<CSTATE> wgbc_k;
  TPZVec<CSTATE> wgbc_f;
};


TPZVec<wgma::gmeshtools::CylinderData> SetUpCylData(std::string_view filename,
                                                    const REAL scale);

wgma::cmeshtools::PhysicalData
FillDataForModalAnalysis(const TPZVec<std::map<std::string, int>> &gmshmats,
                         const SimData& simdata,
                         const std::string &name);

TPZAutoPointer<wgma::wganalysis::Wgma2D>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const SimData& simdata,
  const CSTATE target,
  int &nEigenpairs,
  const TPZEigenSort sortingRule,
  bool usingSLEPC,
  const std::string &name);

std::map<int,int>
SplitMaterialsNearWpbc(const TPZAutoPointer<TPZCompMesh> &modal_mesh,
                       std::set<int> &all_matids);
void
SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
                TPZAutoPointer<wgma::wganalysis::Wgma2D> &src_an,
                TPZAutoPointer<wgma::wganalysis::Wgma2D> &match_an,
                const TPZVec<std::map<std::string, int>> &gmshmats,
                const std::map<int,int> &split_mats,
                const std::map<int64_t,int64_t> &periodic_els,
                const SimData &simdata);
void
SolveModePropagation(TPZAutoPointer<TPZGeoMesh> gmesh,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> &src_an,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> &match_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const std::map<int,int> &split_mats,
                     const std::map<int64_t,int64_t> &periodic_els,
                     const SimData &simdata);

SimData GetSimData()
{

  // path for output files
  const std::string path {"res_sf3d_cyl/"};
  // common prefix for both meshes and output files
  const std::string basisName{"sf3d_cyl"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  
  SimData data;
  data.meshfile="meshes/sf3d_validation.msh";
  data.cylfile="meshes/sf3d_validation_cyldata.csv";
  data.lambda = 4.0;
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
  constexpr REAL pmlval=0.4;
  data.alphaPMLr = {sqrt(pmlval*pmlval+pmlval*pmlval), 0.0};
  data.alphaPMLz = {1.0, 0.0};
  data.check_mode_propagation=true;
  data.compare_pml=false;
  data.source_coeffs = {{0,1}};
  data.porder = 2;
  data.n_eigenpairs_left = 5;
  data.n_eigenpairs_right = 5;
  data.n_modes_left = {5};
  data.n_modes_right = {5};
  data.filter_bnd_eqs = true;
  data.print_gmesh=false;
  data.direct_solver = false;
  data.export_vtk_modes = false;
  data.export_csv_modes = true;
  data.export_csv_error = true;
  data.export_vtk_scatt = true;
  data.export_vtk_error = true;
  data.vtk_res=0;
  data.prefix = prefix;
  data.n_threads=std::thread::hardware_concurrency();
  return std::move(data);
}

int main(int argc, char *argv[]) {

  wgma::wganalysis::using_tbb_mat=true;
  wgma::scattering::using_tbb_mat=true;
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif
  if(argc>2){
    PZError<<"Unexpected number of parameters. USAGE: ./sf3d_cyl param_file\n"
           <<"or ./sf3d_Cyl to run with default"<<std::endl;
    return -1;
  }
  
  SimData simdata= [argc, &argv](){
    if(argc==1){
      auto sd = GetSimData();
      return sd;
    }
    const std::string dataname = argv[1];
    return ReadSimData(dataname);
  }();
  
  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(simdata.prefix));
  
  /******************
   * eigensolver options *
   ******************/

  // how to sort eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::UserDefined};
  constexpr bool usingSLEPC {true};

  const CSTATE target_left{-1.00001*simdata.ncore*simdata.ncore};
  const CSTATE target_right{-1.00001*simdata.ncore*simdata.ncore};
  /*********
   * begin *
   *********/




  /*************
   * geometry  *
   ************/
  // scoped-timer
  TPZSimpleTimer total("Total");

  TPZVec<std::map<std::string, int>> gmshmats;
  constexpr bool verbosity_lvl{false};
  std::map<int64_t,int64_t> periodic_els;
  TPZAutoPointer<TPZGeoMesh> gmesh{nullptr};
  {
    TPZSimpleTimer timer("ReadMesh",true);
    gmesh = wgma::gmeshtools::ReadPeriodicGmshMesh(simdata.meshfile, simdata.scale,
                                                   gmshmats, periodic_els,
                                                   verbosity_lvl);
  }

  constexpr bool cyl3d{true};
  if(cyl3d){
  /**
     in order to exactly represent all the circles in the mesh,
     the python script that generated the mesh also generated this .csv file
    **/
    auto cyldata = SetUpCylData(simdata.cylfile, simdata.scale);
    wgma::gmeshtools::SetExactCylinderRepresentation(gmesh, cyldata);
  }
  
  // print wgma_gmesh to .txt and .vtk format
  if (simdata.print_gmesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = simdata.prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }


  const CSTATE epsilon_clad{simdata.nclad*simdata.nclad};
  const CSTATE epsilon_core{simdata.ncore*simdata.ncore};
  
  //left port modal analysis
  TPZAutoPointer<wgma::wganalysis::Wgma2D>
    modal_l_an{nullptr};
  {
    modal_l_an =  ComputeModalAnalysis(gmesh,gmshmats,simdata,target_left,
                                       simdata.n_eigenpairs_left,sortingRule,
                                       usingSLEPC,"src_left");
  }

  //right port modal analysis
  TPZAutoPointer<wgma::wganalysis::Wgma2D>
    modal_r_an{nullptr};
  {
    modal_r_an =  ComputeModalAnalysis(gmesh,gmshmats,simdata,target_right,
                                       simdata.n_eigenpairs_right,sortingRule,
                                       usingSLEPC,"src_right");
  }

  std::set<int> all_matids;
  for(auto &mats : gmshmats){//dim
    for(auto &[name,id] : mats){//name,matid
      all_matids.insert(id);
    }
  }
  
  auto modal_l_map =
    SplitMaterialsNearWpbc(modal_l_an->GetMesh(),all_matids);
  auto modal_r_map =
    SplitMaterialsNearWpbc(modal_r_an->GetMesh(),all_matids);

  //now we combine the maps but inverting key->value, so we have new_mat->old_mat
  std::map<int,int> split_mats;
  for(auto [old_mat,new_mat] : modal_l_map){
    split_mats[new_mat] = old_mat;
  }
  for(auto [old_mat,new_mat] : modal_r_map){
    split_mats[new_mat] = old_mat;
  }

  //these methods will delete the wgma2d analysis instances!!!
  if(!simdata.check_mode_propagation){
    SolveScattering(gmesh, modal_l_an,  modal_r_an, gmshmats,
                    split_mats,periodic_els, simdata);
  }else{
    SolveModePropagation(gmesh, modal_l_an,  modal_r_an, gmshmats,
                         split_mats,periodic_els, simdata);
  }
  return 0;
}

void
ComputeCouplingMat(wgma::wganalysis::Wgma2D &an,
                   std::string filename,
                   const int nthreads,
                   const bool conj);

void
PostProcessModes(wgma::wganalysis::Wgma2D &an,
                 std::string filename,
                 const int vtkres,
                 const int nthreads);

void
TransformModes(wgma::wganalysis::Wgma2D& an);

TPZAutoPointer<TPZCompMesh>
CreateScattMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                const TPZVec<std::map<std::string, int>> &gmshmats,
                const std::map<int,int> &split_mats,
                std::set<int> &new_mats,
                const SimData &simdata, bool usingPML);

[[nodiscard]] std::set<int>
UpdatePhysicalDataSplittedMats(TPZAutoPointer<TPZGeoMesh> &gmesh,
                               wgma::cmeshtools::PhysicalData& data,
                               const std::map<int,int> &matid_map,
                               const std::set<int> &volids,
                               const int dim);

void
CreateElementGroups(TPZCompMesh *cmesh,const std::set<int>& mat_ids);

void ComputeWpbcCoeffs(wgma::wganalysis::Wgma2D& an,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       const bool positive_z, const TPZVec<CSTATE> &coeff,
                       const int nthreads);

void
TransferSolutionBetweenPeriodicMeshes(TPZAutoPointer<TPZCompMesh> dest_mesh,
                                      TPZAutoPointer<TPZCompMesh> src_mesh,
                                      const std::map<int64_t,int64_t>& periodic_els);

void SolveWithPML(TPZAutoPointer<TPZCompMesh> scatt_cmesh,
                  wgma::wganalysis::Wgma2D& src_an,
                  const TPZVec<CSTATE> &src_coeffs,
                  const SimData &simdata);

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WpbcData& src_data,
                          WpbcData& match_data,
                          const TPZVec<CSTATE> &source_coeffs,
                          const int nmodes_src,
                          const int nmodes_match,
                          const std::set<int> &mats_near_wpbc,
                          const SimData &simdata);

void
AddWaveguidePortContribution(wgma::scattering::Analysis &scatt_an, 
                             const int64_t indep_con_id,
                             const int nm,
                             const TPZFMatrix<CSTATE> &wgbc_k,
                             const TPZVec<CSTATE> &wgbc_f);


STATE ComputeProjAndError(TPZAutoPointer<TPZCompMesh> proj_mesh,
                          TPZAutoPointer<TPZCompMesh> error_mesh,
                          TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          const TPZFMatrix<CSTATE> &projected_modes,
                          const int nm,
                          const std::set<int64_t> &bound_connects,
                          TPZAutoPointer<TPZVTKGenerator> vtk_proj,
                          TPZAutoPointer<TPZVTKGenerator> vtk_error,
                          const std::string &name,
                          const SimData &simdata
                          );


void SetupPrecond(wgma::scattering::Analysis &scatt_an,
                  const std::set<int64_t> &indep_cons);

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

wgma::cmeshtools::PhysicalData
FillDataForModalAnalysis(const TPZVec<std::map<std::string, int>> &gmshmats,
                         const SimData& simdata,
                         const std::string &name)
{
  // setting up cmesh data
  const auto &nclad = simdata.nclad;
  const auto &ncore = simdata.ncore;
  const auto &alphaPMLr = simdata.alphaPMLr;
    
  wgma::cmeshtools::PhysicalData modal_data;
  std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
  modal_mats[name+"_clad"] = std::pair<CSTATE, CSTATE>(nclad*nclad, 1.);
  modal_mats[name+"_core"] = std::pair<CSTATE, CSTATE>(ncore*ncore, 1.);
  std::map<std::string, wgma::bc::type> modal_bcs;
  modal_bcs[name+"_bnd"] = wgma::bc::type::PEC;
  //dimension of the modal analysis 
  constexpr int modal_dim{2};
  wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                          {alphaPMLr,0}, modal_data, modal_dim);
  //we must now filter the 2D PMLs
  std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
  for(const auto &pml : modal_data.pmlvec){
    const std::string pattern{name+"_clad"};
    const auto rx = std::regex{pattern, std::regex_constants::icase };
    
    const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
    if(found_pattern){pmlvec.push_back(pml);}
  }
  modal_data.pmlvec = pmlvec;
  return modal_data;
}

TPZAutoPointer<wgma::wganalysis::Wgma2D>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const SimData& simdata,
  const CSTATE target,
  int &nEigenpairs,
  const TPZEigenSort sortingRule,
  bool usingSLEPC,
  const std::string &name)
{
  auto modal_data =
    FillDataForModalAnalysis(gmshmats,simdata,name);

  const auto &pOrder = simdata.porder;
  const auto &lambda = simdata.lambda;
  const auto &scale = simdata.scale;
  
  auto modal_cmesh = wgma::wganalysis::CMeshWgma2D(gmesh,pOrder,modal_data,
                                       lambda, scale);
  /******************************
   * solve(modal analysis left) *
   ******************************/
  constexpr bool verbose{true};
  const int krylovDim =
    nEigenpairs < 20 ? nEigenpairs*5 : (int)std::ceil(1.25*nEigenpairs);
  auto solver = wgma::wganalysis::SetupSolver(target, nEigenpairs, sortingRule, usingSLEPC,krylovDim,verbose);
  if(sortingRule==TPZEigenSort::UserDefined){
    solver->SetUserSortingFunc([](CSTATE a, CSTATE b)->bool{
      const auto sqrt_a_im = std::fabs(std::imag(std::sqrt(-a)));
      const auto sqrt_b_im = std::fabs(std::imag(std::sqrt(-b)));
      return sqrt_a_im < sqrt_b_im;
    });
  }

  TPZAutoPointer<wgma::wganalysis::Wgma2D> an =
    new wgma::wganalysis::Wgma2D(modal_cmesh, simdata.n_threads,
                                     false,
                                     simdata.filter_bnd_eqs);
  an->SetSolver(*solver);

  std::string modalfile{simdata.prefix+"_modal_"+name};

  {
    TPZSimpleTimer analysis("Modal analysis");
    an->Assemble();
    static constexpr bool computeVectors{true};
  

    {
      TPZSimpleTimer timer("Solve",true);
      an->Solve(computeVectors,verbose);
    }
  
    //load all obtained modes into the mesh
    an->LoadAllSolutions();

    TPZVec<CSTATE> betavec = an->GetEigenvalues();
    nEigenpairs = betavec.size();
    std::cout<<nEigenpairs<<" eigenpairs have converged"<<std::endl;
    for(auto &b : betavec){b = sqrt(-b);}
    if(simdata.export_csv_modes){
      std::ostringstream eigeninfo;
      typedef std::numeric_limits< double > dbl;
      eigeninfo.precision(dbl::max_digits10);
      for(auto &b : betavec){
        const auto pos_sign = std::imag(b) > 0 ? "+" : "-";
        eigeninfo<<std::fixed<<std::real(b)<<pos_sign<<std::abs(std::imag(b))<<"j\n";
      }
      std::ofstream eigenfile(simdata.prefix+"_evalues_"+name+".csv",std::ios::trunc);
      eigenfile<<eigeninfo.str();
      eigenfile.close();
    }
    
    {
      TPZSimpleTimer timer("Ortho",true);
      constexpr STATE tol{1e-14};
      constexpr bool conj{false};
      const int n_ortho = wgma::post::OrthoWgSol(an,tol,conj);
      std::cout<<"orthogonalised  "<<n_ortho<<" degenerate eigenpairs"<<std::endl;
    }

    if(simdata.export_vtk_modes){
      PostProcessModes(*an, modalfile, simdata.vtk_res,simdata.n_threads);
      an->LoadAllSolutions();
    }
    
    if(simdata.couplingmat){

      
      std::string couplingfile{simdata.prefix+"_coupling_"+name+".csv"};
      ComputeCouplingMat(*an,couplingfile,simdata.n_threads,false);
      couplingfile = simdata.prefix+"_coupling_"+name+"_conj.csv";
      ComputeCouplingMat(*an,couplingfile,simdata.n_threads,true);
      an->LoadAllSolutions();
    }

    /*
      in the modal analysis we perform a change of variables
      now we transform back the solutions
    */
    TransformModes(*an);
    
    //now we normalise them
    auto cmesh = an->GetMesh();
    // //leave empty for all valid matids
    // std::set<int> matids {};
    // constexpr bool conj{true};
    // auto norm =
    //   wgma::post::WgNorm<wgma::post::MultiphysicsIntegrator>(cmesh,matids,
    //                                                conj,simdata.n_threads);
    // norm.SetNThreads(simdata.n_threads);    
    // norm.SetBeta(betavec);
    // norm.SetWavelength(simdata.lambda/simdata.scale);
    // norm.Normalise();
    // TPZSimpleTimer timer2("LoadAllSolutions",true);
    // TPZFMatrix<CSTATE> &mesh_sol=cmesh->Solution();
    // //we update analysis object
    // an->SetEigenvectors(mesh_sol);
    an->LoadAllSolutions();
  }
  //we dont need them anymore, let us free up memory
  an->GetSolver().SetMatrixA(nullptr);
  an->GetSolver().SetMatrixB(nullptr);
  return an;
}

//! Reads sim data from file
SimData ReadSimData(const std::string &dataname){
  using json = nlohmann::json;
  std::ifstream f(dataname);
  json data = json::parse(f);
  SimData sd;
  sd.meshfile = data["meshfile"];
  sd.cylfile = data["cylfile"];
  sd.prefix =  data["prefix"];
  sd.lambda =  data["wavelength"];
  sd.scale = data["scale"];
  sd.ncore = data["ncore"];
  sd.nclad = data["nclad"];
  sd.couplingmat = data["export_coupling_mat"];
  sd.check_mode_propagation = data["check_mode_propagation"];
  sd.compare_pml = data["compare_pml"];
  sd.direct_solver = data.value("direct_solver",false);
  std::vector<double> tmpvec_double;
  tmpvec_double = data["alpha_pml_r"].get<std::vector<double>>();
  if(tmpvec_double.size()!=2){
    DebugStop();
  }
  sd.alphaPMLr = {tmpvec_double[0], tmpvec_double[1]};
  tmpvec_double = data["alpha_pml_z"].get<std::vector<double>>();
  if(tmpvec_double.size()!=2){
    DebugStop();
  }
  sd.alphaPMLz = {tmpvec_double[0], tmpvec_double[1]};
  sd.porder = data["porder"];
  sd.source_coeffs = data["source_coeffs"].get<std::vector<std::pair<int,double>>>();
  sd.n_eigenpairs_left = data["n_eigenpairs_left"];
  sd.n_eigenpairs_right = data["n_eigenpairs_right"];
  std::vector<int> tmpvec_int;
  tmpvec_int = data["n_modes_left"].get<std::vector<int>>();
  for(auto nm : tmpvec_int){sd.n_modes_left.push_back(nm);}
  tmpvec_int = data["n_modes_right"].get<std::vector<int>>();
  for(auto nm : tmpvec_int){sd.n_modes_right.push_back(nm);}
  sd.filter_bnd_eqs = data.value("filter_bnd_eqs",true);
  sd.print_gmesh=data.value("print_gmesh",true);
  sd.optimize_bandwidth = data.value("optimize_bandwidth",true);
  sd.export_csv_modes = data.value("export_csv_modes",true);
  sd.export_csv_error = data.value("export_csv_error",true);
  sd.export_vtk_modes = data.value("export_vtk_modes",false);
  sd.export_vtk_scatt = data.value("export_vtk_scatt",true);
  sd.export_vtk_error = data.value("export_vtk_error",true);
  sd.vtk_res = data.value("vtk_res",(int)0);
  sd.n_threads = data.value("n_threads",(int)std::thread::hardware_concurrency());
  return sd;
}


void SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> &src_an,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> &match_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const std::map<int,int> &split_mats,
                     const std::map<int64_t,int64_t> &periodic_els,
                     const SimData &simdata)
{

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
  for(auto [i, alpha] : simdata.source_coeffs){
    if(i >= src_coeffs.size()){
      std::cout<<"ERROR: src coefficient bigger than computed number of modes\n"
               <<"i: "<<i<<" alpha "<<alpha<<std::endl;
      DebugStop();
    }
    src_coeffs[i] = alpha;
  }


  //index of the number of modes to be used to restrict the dofs on waveguide bcs
  auto nmodes_left = simdata.n_modes_left;
  auto nmodes_right = simdata.n_modes_right;
  {
    //we need to check the highest order mode used as a source
    int last_mode = -1;
    for(int i = 0; i < src_coeffs.size();i++){
      if(src_coeffs[i] != (CSTATE)0){
        last_mode = i;
      }
    }

    TPZVec<int> valid_nmodes_left;
    for (auto im : nmodes_left){
      if(im > last_mode){
        valid_nmodes_left.push_back(im);
      }
    }
    nmodes_left = valid_nmodes_left;
  }

  
  src_an->LoadAllSolutions();
  match_an->LoadAllSolutions();

  //set up post processing
  TPZVec<std::string> fvars_3d = {
    "Field_real",
    "Field_imag",
    "Field_abs"};
  
  if(simdata.compare_pml){
    std::set<int> mats_near_wpbc;
    TPZAutoPointer<TPZCompMesh> scatt_mesh_pml =
      CreateScattMesh(gmesh,gmshmats,split_mats,mats_near_wpbc,simdata,true);
    const std::string suffix = "pml";
    const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
    SolveWithPML(scatt_mesh_pml,src_an,src_coeffs,simdata);
    //now we load back all the solutions into the src mesh
    src_an->LoadAllSolutions();
    if(simdata.export_vtk_scatt){
      auto vtk = TPZVTKGenerator(scatt_mesh_pml, fvars_3d, scatt_file, simdata.vtk_res);
      vtk.SetNThreads(simdata.n_threads);
      vtk.Do();
    }
  }



  std::set<int> mats_near_wpbc;
  TPZAutoPointer<TPZCompMesh> scatt_mesh_wpbc =
    CreateScattMesh(gmesh,gmshmats,split_mats,mats_near_wpbc,simdata,false);

  const std::string suffix = "wpbc";
  const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
  auto vtk = TPZVTKGenerator(scatt_mesh_wpbc, fvars_3d, scatt_file, simdata.vtk_res);
  vtk.SetNThreads(simdata.n_threads);
  

  
  //compute wgbc coefficients
  WpbcData src_data;
  WpbcData match_data;

  {
    TPZSimpleTimer timer("wpbc coeffs",true);
    src_data.cmesh = src_an->GetHCurlMesh();
    ComputeWpbcCoeffs(src_an,  src_data.wgbc_k,
                      src_data.wgbc_f, false, src_coeffs,
                      simdata.n_threads);

    //only hcurl mesh is needed from now on, we can delete analysis object
    src_an=nullptr;
    match_data.cmesh = match_an->GetHCurlMesh();
    ComputeWpbcCoeffs(match_an, match_data.wgbc_k,
                      match_data.wgbc_f,true, {},
                      simdata.n_threads);
    match_an=nullptr;
  }

  //set up post processing vars
  TPZVec<std::string> fvars_2d = {
    "Solution_abs",
    "Solution_real",
    "Solution_imag"};
  

  auto CreateErrorMesh = [&simdata, &gmesh, &gmshmats, &fvars_2d,
                          &periodic_els](TPZAutoPointer<TPZCompMesh> &error_mesh,
                                         const TPZAutoPointer<TPZCompMesh> &hcurl_mesh,
                                         TPZFMatrix<CSTATE> &projected_modes,
                                         TPZAutoPointer<TPZVTKGenerator> &vtk_error,
                                         const std::string &name)
  {
    auto modal_data = FillDataForModalAnalysis(gmshmats, simdata, name);
    std::set<int> volmats, pmlmats;
    wgma::wganalysis::SetupModalAnalysisMaterials(gmesh,modal_data,volmats,pmlmats);

    //we want hcurl mesh
    constexpr bool is_h1{false};
    error_mesh =
      wgma::wganalysis::CreateAtomicWgma2D(gmesh,is_h1,simdata.porder,volmats,pmlmats,
                                           modal_data.bcvec, modal_data.probevec);
    //now we transfer the modal solution from the WPBC to the error mesh and store it
    TransferSolutionBetweenPeriodicMeshes(error_mesh, hcurl_mesh, periodic_els);
    projected_modes = error_mesh->Solution();
    const int neq = error_mesh->Solution().Rows();
    error_mesh->Solution().Resize(neq, 1);
    if(simdata.export_vtk_error){
      const std::string error_file = simdata.prefix+"_error_"+name;
      vtk_error = new TPZVTKGenerator(error_mesh, fvars_2d, error_file, simdata.vtk_res);
      vtk_error->SetNThreads(simdata.n_threads);
    }
  };


  /*
    we want to analyse the solution close to the ports and see if the
    number of modes is sufficient to represent it.
   */
  TPZAutoPointer<TPZCompMesh> error_mesh_left{nullptr}, error_mesh_right{nullptr};
  TPZFMatrix<CSTATE> projected_modes_left, projected_modes_right;
  TPZAutoPointer<TPZVTKGenerator> vtk_error_left{nullptr}, vtk_error_right{nullptr};
  CreateErrorMesh(error_mesh_left,src_data.cmesh,projected_modes_left,vtk_error_left,"probe_left");
  CreateErrorMesh(error_mesh_right,match_data.cmesh,projected_modes_right,vtk_error_right,"probe_right");
  
  

  auto SetupProjMesh = [&simdata, &fvars_2d](TPZAutoPointer<TPZCompMesh> &proj_mesh,
                                          TPZAutoPointer<TPZCompMesh> &error_mesh,
                                          std::set<int64_t> &bound_connects,
                                          TPZAutoPointer<TPZVTKGenerator> &vtk_proj,
                                          const std::string &name){
    // this will be the restricted mesh close to the wg port
    proj_mesh = error_mesh->Clone();
    wgma::cmeshtools::FindDirichletConnects(proj_mesh, bound_connects);
    if(simdata.export_vtk_error){
      const std::string proj_file = simdata.prefix+"_proj_"+name;
      vtk_proj = new TPZVTKGenerator(proj_mesh,fvars_2d, proj_file, simdata.vtk_res);
      vtk_proj->SetNThreads(simdata.n_threads);
    }
  };
  /*
    dirichlet boundary connects should not be restricted, otherwise
    this will result in all the equations on the same dependency
    being removed as well
  */
  std::set<int64_t> bound_connects_left, bound_connects_right;
  TPZAutoPointer<TPZCompMesh> proj_mesh_left{nullptr}, proj_mesh_right{nullptr};
  TPZAutoPointer<TPZVTKGenerator> vtk_proj_left{nullptr}, vtk_proj_right{nullptr};
  
  SetupProjMesh(proj_mesh_left,error_mesh_left,bound_connects_left,vtk_proj_left,"left");
  SetupProjMesh(proj_mesh_right,error_mesh_right,bound_connects_right,vtk_proj_right,"right");

  std::map<std::pair<int,int>,std::pair<STATE,STATE>> wpbc_error_res;


  std::ostringstream errorinfo;
  errorinfo.precision(std::numeric_limits<STATE>::max_digits10);
  
  for(int im_left = 0; im_left < nmodes_left.size(); im_left++){
    const int nm_left = nmodes_left[im_left];
    for(int im_right = 0; im_right < nmodes_right.size(); im_right++){
      const int nm_right = nmodes_right[im_right];
      RestrictDofsAndSolve(scatt_mesh_wpbc, src_data, match_data,
                           src_coeffs, nm_left,nm_right,
                           mats_near_wpbc,simdata);
      //plot
      if(simdata.export_vtk_scatt){vtk.Do();}
      const auto error_left =
        ComputeProjAndError(proj_mesh_left,error_mesh_left,scatt_mesh_wpbc,
                            projected_modes_left,nm_left,bound_connects_left,
                            vtk_proj_left,vtk_error_left,"left",
                            simdata);
      const auto error_right = 
        ComputeProjAndError(proj_mesh_right,error_mesh_right,scatt_mesh_wpbc,
                            projected_modes_right,nm_right,bound_connects_right,
                            vtk_proj_right,vtk_error_right,"right",
                            simdata);
      wpbc_error_res.insert({{nm_left,nm_right},{error_left,error_right}});

      if(simdata.export_csv_error){
        errorinfo<<nm_left<<','<<nm_right<<','
                 <<std::fixed<<error_left<<','
                 <<std::fixed<<error_right<<'\n';
      }
      //removing restrictions
      wgma::cmeshtools::RemovePeriodicity(scatt_mesh_wpbc);
      scatt_mesh_wpbc->ComputeNodElCon();
      scatt_mesh_wpbc->CleanUpUnconnectedNodes();
    }
  }
  if(simdata.export_csv_error){
    std::ofstream errorfile(simdata.prefix+"_error.csv",std::ios::trunc);
    errorfile<<"nm(left),nm(right),error(left),error(right)"<<std::endl;
    errorfile<<errorinfo.str();
    errorinfo.clear();
  }
  std::cout<<"+++++++++++++WPBC COMPARISON+++++++++++++"<<std::endl;
  for(auto [nm,error] : wpbc_error_res){
    std::cout<<"nmodes (left) "<<nm.first<<" nmodes (right) "<<nm.second
             <<" norm error (left): "<<error.first<<" norm error (right): "<<error.second <<std::endl;
  }
  std::cout<<"+++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
  
}

void SolveModePropagation(TPZAutoPointer<TPZGeoMesh> gmesh,
                          TPZAutoPointer<wgma::wganalysis::Wgma2D> &src_an,
                          TPZAutoPointer<wgma::wganalysis::Wgma2D> &match_an,
                          const TPZVec<std::map<std::string, int>> &gmshmats,
                          const std::map<int,int> &split_mats,
                          const std::map<int64_t,int64_t> &periodic_els,
                          const SimData &simdata)
{
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
  for(auto [i, alpha] : simdata.source_coeffs){
    if(i >= src_coeffs.size()){
      std::cout<<"ERROR: src coefficient bigger than computed number of modes\n"
               <<"i: "<<i<<" alpha "<<alpha<<std::endl;
      DebugStop();
    }
    src_coeffs[i] = alpha;
  }


  //index of the number of modes to be used to restrict the dofs on waveguide bcs
  auto nmodes_left = simdata.n_modes_left;
  auto nmodes_right = simdata.n_modes_right;
  {
    //we need to check the highest order mode used as a source
    int last_mode = -1;
    for(int i = 0; i < src_coeffs.size();i++){
      if(src_coeffs[i] != (CSTATE)0){
        last_mode = i;
      }
    }

    std::set<int> valid_nmodes_left;
    for (auto im : nmodes_left){
      if(im > last_mode){
        valid_nmodes_left.insert(im);
      }
    }
    const int nmodes_left_size = valid_nmodes_left.size();
    nmodes_left.resize(0);
    for(auto im : valid_nmodes_left){
      nmodes_left.push_back(im);
    }
  }

  /*
    first thing is to compute the propagation of the waveguide mode in the
    error mesh
   */
  TPZFMatrix<CSTATE> ref_sol;
  STATE norm_sol_ref{1};
  
  TPZAutoPointer<TPZCompMesh> error_mesh{nullptr};
  TPZAutoPointer<TPZVTKGenerator> vtk_proj_eval{nullptr};

  //set up post processing vars
  TPZVec<std::string> fvars_2d = {
    "Solution_abs",
    "Solution_real",
    "Solution_imag"};

  {

    {
      auto modal_data = FillDataForModalAnalysis(gmshmats, simdata, "probe_left");
      std::set<int> volmats, pmlmats;
      wgma::wganalysis::SetupModalAnalysisMaterials(gmesh,modal_data,volmats,pmlmats);

      //we want hcurl mesh
      constexpr bool is_h1{false};
      error_mesh =
        wgma::wganalysis::CreateAtomicWgma2D(gmesh,is_h1,simdata.porder,volmats,pmlmats,
                                             modal_data.bcvec, modal_data.probevec);
    }
    //just to set size
    ref_sol = error_mesh->Solution();

    //now error mesh will contain all the modes
    TransferSolutionBetweenPeriodicMeshes(error_mesh, src_an->GetHCurlMesh(), periodic_els);
    //let us get the distance between src and error mesh
    const STATE dist =
      error_mesh->Element(0)->Reference()->Node(0).Coord(2)-
      src_an->GetHCurlMesh()->Element(0)->Reference()->Node(0).Coord(2);
    //now we compute the expected (propagated) solution
    ref_sol.Redim(ref_sol.Rows(),1);

    
    for(int i = 0; i < src_coeffs.size();i++){
      CSTATE *ref_sol_ptr = ref_sol.Elem();
      if(src_coeffs[i] != (CSTATE)0){
        TPZFMatrix<CSTATE> &mode = error_mesh->Solution();
        const auto neq = mode.Rows();
        const auto offset = neq*i;
        const auto beta = std::sqrt(-src_an->GetEigenvalues()[i]);
        const auto coeff = src_coeffs[i]*std::exp(-1i*beta*dist);
        std::cout<<"beta "<<beta<<" dist "<<dist<<std::endl;
        CSTATE *mode_ptr = mode.Elem() + offset;
        for(int ieq = 0; ieq < neq; ieq++){
          *ref_sol_ptr++ += coeff*(*mode_ptr++);
        }
      }
    }

    error_mesh->LoadSolution(ref_sol);
    auto normsol =
      wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
    normsol.SetNThreads(simdata.n_threads);
    norm_sol_ref = std::real(normsol.ComputeNorm()[0]);
    if(simdata.export_vtk_error){
      const std::string proj_eval_file = simdata.prefix+"_proj_eval";
      TPZVTKGenerator vtk_proj_eval(error_mesh, fvars_2d, proj_eval_file, simdata.vtk_res);
      vtk_proj_eval.SetNThreads(simdata.n_threads);
      vtk_proj_eval.Do();
    }
    
  }

  
  //set up post processing vars
  TPZVec<std::string> fvars_3d = {
    "Field_real",
    "Field_imag",
    "Field_abs"};

  std::set<int> mats_near_wpbc;
  TPZAutoPointer<TPZCompMesh> scatt_mesh_wpbc =
    CreateScattMesh(gmesh,gmshmats,split_mats,mats_near_wpbc,simdata,false);


  const std::string suffix = "wpbc";
  const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
  auto vtk = TPZVTKGenerator(scatt_mesh_wpbc, fvars_3d, scatt_file, simdata.vtk_res);
  vtk.SetNThreads(simdata.n_threads);

  TPZFMatrix<CSTATE> sol_pml;
  STATE norm_sol_pml{1}, norm_error_pml{1};
  if(simdata.compare_pml){
    //ok first we solve with a PML
    std::set<int> mats_near_wpbc;
    TPZAutoPointer<TPZCompMesh> scatt_mesh_pml =
      CreateScattMesh(gmesh,gmshmats,split_mats,mats_near_wpbc,simdata,true);
    const std::string suffix = "pml";
    const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
    SolveWithPML(scatt_mesh_pml,src_an,src_coeffs,simdata);
    //now we load back all the solutions into the src mesh
    src_an->LoadAllSolutions();
    if(simdata.export_vtk_scatt){
      auto vtk = TPZVTKGenerator(scatt_mesh_pml, fvars_3d, scatt_file, simdata.vtk_res);
      vtk.SetNThreads(simdata.n_threads);
      vtk.Do();
    }
    //just to set size
    sol_pml = error_mesh->Solution();

    wgma::cmeshtools::ExtractSolFromMesh(error_mesh, scatt_mesh_pml, sol_pml);
    error_mesh->LoadSolution(sol_pml);
    
    auto normsol =
      wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
    normsol.SetNThreads(simdata.n_threads);
    norm_sol_pml = std::real(normsol.ComputeNorm()[0]);
    std::cout<<"norm sol pml: "<<norm_sol_pml<<std::endl;

    //export PML solution at error mesh
    if(simdata.export_vtk_error){
      const std::string pml_eval_file = simdata.prefix+"_pml_eval";
      TPZVTKGenerator vtk_pml_eval(error_mesh, fvars_2d, pml_eval_file, simdata.vtk_res);
      vtk_pml_eval.SetNThreads(simdata.n_threads);
      vtk_pml_eval.Do();
    }
    auto norm_error =
        wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
    norm_error.SetNThreads(simdata.n_threads);
    norm_error_pml = std::real(norm_error.ComputeNorm()[0])/norm_sol_ref;
    std::cout<<"norm error pml: "<<norm_error_pml<<std::endl;
    
    if(simdata.export_vtk_error){
      sol_pml-=ref_sol;
      const std::string error_file = simdata.prefix+"_pml_error";
      TPZVTKGenerator vtk_error(error_mesh, fvars_2d, error_file, simdata.vtk_res);
      vtk_error.SetNThreads(simdata.n_threads);
      vtk_error.Do();
    }
  }


  src_an->LoadAllSolutions();
  match_an->LoadAllSolutions();

  

  //compute wgbc coefficients
  WpbcData src_data;
  WpbcData match_data;

  {
    TPZSimpleTimer timer("wpbc coeffs",true);
    src_data.cmesh = src_an->GetHCurlMesh();
    ComputeWpbcCoeffs(src_an,  src_data.wgbc_k,
                      src_data.wgbc_f, false, src_coeffs,
                      simdata.n_threads);

    //only hcurl mesh is needed from now on, we can delete analysis object
    src_an=nullptr;
    match_data.cmesh = match_an->GetHCurlMesh();
    ComputeWpbcCoeffs(match_an, match_data.wgbc_k,
                      match_data.wgbc_f,true, {},
                      simdata.n_threads);
    match_an = nullptr;
  }
  
  TPZFMatrix<CSTATE> sol_wpbc = ref_sol;//just to have the same size
  
  std::map<std::pair<int,int>,STATE> wpbc_error_res;

  std::ostringstream errorinfo;
  errorinfo.precision(std::numeric_limits<STATE>::max_digits10);


  TPZAutoPointer<TPZVTKGenerator> vtk_wpbc_eval{nullptr}, vtk_wpbc_error{nullptr};
  if(simdata.export_vtk_error){
    const std::string wpbc_eval_file = simdata.prefix+"_wpbc_eval";
    vtk_wpbc_eval = new TPZVTKGenerator(error_mesh, fvars_2d, wpbc_eval_file, simdata.vtk_res);
    const std::string wpbc_error_file = simdata.prefix+"_wpbc_error";
    vtk_wpbc_error = new TPZVTKGenerator(error_mesh, fvars_2d, wpbc_error_file, simdata.vtk_res);
  }
  
  for(int im_left = 0; im_left < nmodes_left.size(); im_left++){
    const int nm_left = nmodes_left[im_left];
    for(int im_right = 0; im_right < nmodes_right.size(); im_right++){
      const int nm_right = nmodes_right[im_right];
      RestrictDofsAndSolve(scatt_mesh_wpbc, src_data, match_data,
                           src_coeffs, nm_left,nm_right,
                           mats_near_wpbc,simdata);
      //plot
      if(simdata.export_vtk_scatt){vtk.Do();}

      wgma::cmeshtools::ExtractSolFromMesh(error_mesh, scatt_mesh_wpbc, sol_wpbc);
      error_mesh->LoadSolution(sol_wpbc);
      if(simdata.export_vtk_error){vtk_wpbc_eval->Do();}
      sol_wpbc -= ref_sol;
      error_mesh->LoadSolution(sol_wpbc);
      if(simdata.export_vtk_error){vtk_wpbc_error->Do();}
      auto normsol =
        wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
      normsol.SetNThreads(simdata.n_threads);
      const auto norm = std::real(normsol.ComputeNorm()[0]);
      const auto error = norm/norm_sol_ref;
      std::cout<<"nmodes (left) "<<nm_left<< " nmodes(right) "<<nm_right
               <<" error "<<error<<std::endl;
      wpbc_error_res.insert({{nm_left,nm_right},error});

      if(simdata.export_csv_error){
        errorinfo<<nm_left<<','<<nm_right<<','
                 <<std::fixed<<error<<'\n';
      }
      //removing restrictions
      wgma::cmeshtools::RemovePeriodicity(scatt_mesh_wpbc);
      scatt_mesh_wpbc->ComputeNodElCon();
      scatt_mesh_wpbc->CleanUpUnconnectedNodes();
    }
  }
  if(simdata.export_csv_error){
    std::ofstream errorfile(simdata.prefix+"_error.csv",std::ios::trunc);
    errorfile<<"nm(left),nm(right),error"<<std::endl;
    errorfile<<errorinfo.str();
    errorinfo.clear();
  }

  if(simdata.compare_pml){
    std::cout<<"+++++++++++++PML COMPARISON+++++++++++++"<<std::endl;
    for(auto [nm,error] : wpbc_error_res){
      std::cout<<"nmodes (left) "<<nm.first<<" nmodes (right) "<<nm.second
               <<" norm rel error: "<<error <<" percentage PML error: "<<(error/norm_error_pml)*100<<std::endl;
    }
    std::cout<<"+++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  }
  
  
}


void ComputeCouplingMat(wgma::wganalysis::Wgma2D &an,
                        std::string filename,
                        const int nthreads,
                        const bool conj)
{
  
  using namespace wgma::post;

  std::set<int> matids;
  WaveguideCoupling<MultiphysicsIntegrator> integrator(an.GetMesh(),
                                                       matids,
                                                       conj,
                                                       nthreads
                                                       );
  integrator.SetNThreads(nthreads);
  integrator.ComputeCoupling();
  TPZFMatrix<CSTATE> couplingmat;
  integrator.GetCoupling(couplingmat);
  std::ofstream matfile(filename);
  couplingmat.Print("",matfile,ECSV);
}

void TransformModes(wgma::wganalysis::Wgma2D& an)
{
  using namespace std::complex_literals;
  TPZAutoPointer<TPZCompMesh> h1_mesh = an.GetH1Mesh();
  TPZAutoPointer<TPZCompMesh> hcurl_mesh = an.GetHCurlMesh();
  TPZAutoPointer<TPZCompMesh> mf_mesh = an.GetMesh();
  TPZFMatrix<CSTATE> &hcurl_sol = hcurl_mesh->Solution();
  TPZFMatrix<CSTATE> &h1_sol = h1_mesh->Solution();

  TPZVec<CSTATE> betavec = an.GetEigenvalues();
  for(auto &b : betavec){b = std::sqrt(-b);}
    
  const int nsol = hcurl_sol.Cols();
  {
    const int nrow = hcurl_sol.Rows();
    for(int isol = 0; isol < nsol; isol++){
      const auto beta = betavec[isol];
      CSTATE *sol_ptr = &hcurl_sol.g(0,isol);
      for(int irow = 0; irow < nrow; irow++){
        *sol_ptr++ /=beta;
      }
    }
  }
  {
    const int nrow = h1_sol.Rows();
    for(int isol = 0; isol < nsol; isol++){
      CSTATE *sol_ptr = &h1_sol.g(0,isol);
      for(int irow = 0; irow < nrow; irow++){
        *sol_ptr++ *= 1i;
      }
    }
  }
    
  TPZManVector<TPZAutoPointer<TPZCompMesh>,2> meshvec(2);
  meshvec[TPZWgma::H1Index()] = h1_mesh;
  meshvec[TPZWgma::HCurlIndex()] = hcurl_mesh;    
  TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec,mf_mesh);
  //we update analysis object
  TPZFMatrix<CSTATE> &mesh_sol=mf_mesh->Solution();
  an.SetEigenvectors(mesh_sol);
  an.LoadAllSolutions();
}


TPZAutoPointer<TPZCompMesh>
CreateScattMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                const TPZVec<std::map<std::string, int>> &gmshmats,
                const std::map<int,int> &split_mats,
                std::set<int> &mats_near_wpbc,
                const SimData &simdata, bool usingPML)
{

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
    if(usingPML){
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
    if (usingPML){
      const std::string srcMat[] = {"src_left_core", "src_left_clad"};
      constexpr auto matdim{2};
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
      probeMats.push_back("probe_left_clad");
      probeMats.push_back("probe_left_core");
      probeMats.push_back("probe_right_clad");
      probeMats.push_back("probe_right_core");
      //now for 1d pml mats
      for(const auto &pml : scatt_data.pmlvec){
        const std::string pattern_left{"src_left_clad"};
        const auto rx_left =
          std::regex{pattern_left, std::regex_constants::icase };
        const std::string pattern_right{"src_right_clad"};
        const auto rx_right =
          std::regex{pattern_right, std::regex_constants::icase };
        const std::string pattern_probe_left{"probe_left_clad"};
        const auto rx_probe_left =
          std::regex{pattern_probe_left, std::regex_constants::icase };
        const std::string pattern_probe_right{"probe_right_clad"};
        const auto rx_probe_right =
          std::regex{pattern_probe_right, std::regex_constants::icase };
        const bool found_pattern =
          std::regex_search(*(pml->names.begin()), rx_left) ||
          std::regex_search(*(pml->names.begin()), rx_right) ||
          std::regex_search(*(pml->names.begin()), rx_probe_left) ||
          std::regex_search(*(pml->names.begin()), rx_probe_right);

        if(found_pattern){
          probeMats.push_back(*pml->names.begin());
        }
      }
      for(auto &mat : probeMats){
        const auto matdim = 2;
        const auto id = gmshmats[matdim].at(mat);
        if(src_ids.count(id)==0){
          scatt_data.probevec.push_back({id,matdim});
        }
      }
    }
    
    

    
    
    //due to the waveguide port bc we need to filter out some PML regions
    if (!usingPML){
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

    constexpr int dim{3};
    std::set<int> volids;
    for(auto [id,dummy1,dummy2] : scatt_data.matinfovec){
      volids.insert(id);
    }
    
    mats_near_wpbc =
      UpdatePhysicalDataSplittedMats(gmesh, scatt_data, split_mats,
                                     volids, dim);
    auto cmesh =
      wgma::scattering::CMeshScattering3D(gmesh, pOrder, scatt_data,src_ids,
                                          lambda,scale,true);
    return cmesh;
 };

void PostProcessModes(wgma::wganalysis::Wgma2D &an,
                      std::string filename,
                      const int vtkres,
                      const int nthreads){
  TPZSimpleTimer postProc("Post processing");
    
  const std::string file = filename;
  TPZVec<std::string> fvars = {
      "Ez_real",
      "Ez_abs",
      "Et_real",
      "Et_abs"};
  auto cmesh = an.GetMesh();
  auto vtk = TPZVTKGenerator(cmesh, fvars, file, vtkres);
  vtk.SetNThreads(nthreads);

  std::set<int> sols;
  const auto nsol = std::min((int64_t)20,an.GetEigenvectors().Cols());
  for(auto is = 0; is < nsol; is++){
    sols.insert(is);
  }
  
  std::cout<<"Exporting "<<sols.size()<<" solutions"<<std::endl;
  for(auto isol : sols){
    an.LoadSolution(isol);
    vtk.Do();
  }
}

void AddWaveguidePortContribution(wgma::scattering::Analysis &scatt_an, 
                                  const int64_t indep_con_id,
                                  const int nm,
                                  const TPZFMatrix<CSTATE> &wgbc_k,
                                  const TPZVec<CSTATE> &wgbc_f)
{
  auto mat = scatt_an.GetSolver().Matrix();
  TPZFMatrix<CSTATE>& fvec = scatt_an.Rhs();
  auto scatt_mesh = scatt_an.GetMesh();
  const auto &indep_con = scatt_mesh->ConnectVec()[indep_con_id];
  const auto &block = scatt_mesh->Block();
  const auto seqnum = indep_con.SequenceNumber();
  const auto pos_orig = block.Position(seqnum);
  const auto sz = block.Size(seqnum);
  if(sz!=nm){DebugStop();}

  int64_t pos_filt{-1};
  {
    TPZManVector<int64_t,1> posvec_orig(1,pos_orig), posvec_filt(1,pos_orig);
    scatt_an.StructMatrix()->EquationFilter().Filter(posvec_orig, posvec_filt);
    pos_filt = posvec_filt[0];
    //they are sequential
  }


  TPZManVector<int64_t,600> src_index(nm,0), dest_index(nm,0);
  std::iota(src_index.begin(),src_index.end(),0);
  std::iota(dest_index.begin(),dest_index.end(),pos_orig);

  
  TPZFMatrix<CSTATE> dummy_rhs(nm,1,const_cast<CSTATE*>(wgbc_f.begin()),nm);
  fvec.AddFel(dummy_rhs,src_index,dest_index);
  
  std::iota(dest_index.begin(),dest_index.end(),pos_filt);
  mat->AddKel(const_cast<TPZFMatrix<CSTATE>&>(wgbc_k),src_index,dest_index);
}

void ComputeWpbcCoeffs(wgma::wganalysis::Wgma2D& an,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       const bool positive_z, const TPZVec<CSTATE> &coeff,
                       const int nthreads){
  auto mesh = an.GetMesh();

  TPZFMatrix<CSTATE>& sol_orig = mesh->Solution();
  const int neq = sol_orig.Rows();
  const int nsol = sol_orig.Cols();
  
  wgma::post::WaveguidePortBC<wgma::post::MultiphysicsIntegrator> wgbc(mesh);
  wgbc.SetNThreads(nthreads);
  TPZManVector<CSTATE,1000> betavec = an.GetEigenvalues();
  for(auto &b : betavec){b = std::sqrt(-b);}
  if(coeff.size()){
    wgbc.SetSrcCoeff(coeff);
  }
  wgbc.SetBeta(betavec);
  wgbc.ComputeContribution();
  wgbc.GetContribution(wgbc_k,wgbc_f);
}
#include <TPZSpStructMatrix.h>
#include <TPZStructMatrixOMPorTBB.h>
#include "TPZYSMPMatrix.h"
void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WpbcData& src_data,
                          WpbcData& match_data,
                          const TPZVec<CSTATE> &source_coeffs,
                          const int nmodes_src,
                          const int nmodes_match,
                          const std::set<int> &mats_near_wpbc,
                          const SimData &simdata)
{

  auto gmesh = scatt_mesh->Reference();
  

  auto match_mesh = match_data.cmesh;
  auto src_mesh = src_data.cmesh;

  /**
     no dirichlet connects must be restricted!!
   **/
  std::set<int64_t> boundConnects;
  wgma::cmeshtools::FindDirichletConnects(scatt_mesh, boundConnects);



  const int64_t indep_con_id_match =
    wgma::cmeshtools::RestrictDofs(scatt_mesh, match_mesh, nmodes_match, boundConnects);
  const int64_t indep_con_id_src =
    wgma::cmeshtools::RestrictDofs(scatt_mesh, src_mesh, nmodes_src, boundConnects);

  constexpr bool group{true};
  if(group){
    //already calls expandsolution
    CreateElementGroups(scatt_mesh.operator->(), mats_near_wpbc);
  }
  
  constexpr bool sym{false};
  //either we solve by iterative method or we send it to pardiso to order, so...
  constexpr bool optimize_bandwidth{false};
  auto scatt_an = wgma::scattering::Analysis(scatt_mesh, simdata.n_threads,
                                             optimize_bandwidth,
                                             simdata.filter_bnd_eqs,
                                             sym);


  std::cout<<"nmodes on source boundary: "<<nmodes_src<<std::endl;
  std::cout<<"nmodes on outgoing boundary: "<<nmodes_match<<std::endl;

  //we precomputed it already
  scatt_an.StructMatrix()->SetComputeRhs(false);

  auto strmtrx =
    TPZAutoPointerDynamicCast<TPZSpStructMatrix<CSTATE,
                                                TPZStructMatrixOMPorTBB<CSTATE>>>(
                                                  scatt_an.StructMatrix());
  if(!strmtrx){DebugStop();}
  //300 is the number of shape functions in a k4 hexahedron
  const int maxsz = std::max(nmodes_match,nmodes_src) + 540;
  const int bufsz = maxsz*maxsz;
  strmtrx->BufferSizeForUserMatrix(bufsz);
    
  std::cout<<"Assembling..."<<std::endl;
  TPZSimpleTimer timer("WPBC:Assemble+solve",true);
  {
    std::set<int> mat_ids;
    for(auto [id,mat]: scatt_mesh->MaterialVec()){
      auto nullmat = dynamic_cast<TPZNullMaterial<CSTATE>*>(mat);
      if(!nullmat){
        mat_ids.insert(id);
      }
    }
    TPZSimpleTimer tassemble("Assemble",true);
    scatt_an.Assemble();
  
    //for now we unwrap the groups as they seem to interfere with the solving stage
    if(group){
      const auto nel = scatt_mesh->ElementVec().NElements();
      for(auto index = 0; index < nel; index++){
        auto cel = scatt_mesh->ElementVec()[index];
        auto group = dynamic_cast<TPZElementGroup*>(cel);
        if(group){
          //this call will delete the element groups
          group->Unwrap();
          scatt_mesh->ElementVec()[index] = nullptr;
        }
      }
    }

    //now that we have removed the groups we can get a new sparse matrix
    //with reduced sparsity
    TPZSimpleTimer transf_mat("New sparse matrix");
    auto sparse_old =
      TPZAutoPointerDynamicCast<TPZFYsmpMatrix<CSTATE>>(scatt_an.GetSolver().Matrix());
    auto sparse_new = 
      dynamic_cast<TPZFYsmpMatrix<CSTATE>*>(strmtrx->Create());
    int64_t *ia_old,*ja_old,*ia_new,*ja_new;
    CSTATE *aa_old,*aa_new;
    sparse_old->GetData(ia_old, ja_old, aa_old);
    sparse_new->GetData(ia_new, ja_new, aa_new);
    const auto nr = sparse_old->Rows();
    for(auto ir = 0; ir < nr; ir++){
      const auto first_new = ia_new[ir];
      const auto last_new = ia_new[ir+1];
      int64_t ij_old = ia_old[ir];
      for(auto ij = first_new; ij < last_new; ij++){
        const auto col = ja_new[ij];
        while(ja_old[ij_old] < col){ij_old++;}
#ifdef PZDEBUG
        if(ja_old[ij_old] != col){
          DebugStop();
        }
#endif
        aa_new[ij] = aa_old[ij_old];
      }
    }
    scatt_an.GetSolver().SetMatrix(sparse_new);
  }

  {
    TPZSimpleTimer twpbc("AddWPBC",true);
    //now we must add the waveguide port terms
    AddWaveguidePortContribution(scatt_an, indep_con_id_match,
                                 nmodes_match, match_data.wgbc_k, match_data.wgbc_f);
    AddWaveguidePortContribution(scatt_an, indep_con_id_src,
                                 nmodes_src, src_data.wgbc_k, src_data.wgbc_f);
  }
  
  
  if(simdata.direct_solver){
    auto *pardiso = scatt_an.GetSolver().GetPardisoControl();
    pardiso->SetMatrixType(SymProp::NonSym, TPZPardisoSolver<CSTATE>::MProperty::EIndefinite);
    auto param = pardiso->GetParam();
    param[4] = 2;
    param[9] = 15;
    pardiso->SetParam(param);
  }
  else{
    SetupPrecond(scatt_an, {indep_con_id_src, indep_con_id_match});
  }
  TPZSimpleTimer tsolve("Solve",true);
  scatt_an.Solve();
}

void
TransferSolutionBetweenPeriodicMeshes(TPZAutoPointer<TPZCompMesh> dest_mesh,
                                      TPZAutoPointer<TPZCompMesh> src_mesh,
                                      const std::map<int64_t,int64_t>& periodic_els)
{
  TPZGeoMesh *gmesh = dest_mesh->Reference();
  const TPZFMatrix<CSTATE> &src_sol = src_mesh->Solution();
  const auto &src_block = src_mesh->Block();
  const int nrow = src_sol.Rows();
  const int ncol = src_sol.Cols();
  TPZFMatrix<CSTATE> &dest_sol = dest_mesh->Solution();
  const auto &dest_block = dest_mesh->Block();
  dest_sol.Redim(nrow,ncol);
  //gmesh will point to src_mesh
  src_mesh->LoadReferences();

  const auto &elvec = dest_mesh->ElementVec();
  const int nel = elvec.NElements();
  
  pzutils::ParallelFor(0,nel,[&gmesh, &dest_mesh, &src_mesh, &elvec,
                              &periodic_els, &src_block, &dest_block,
                              &src_sol, &dest_sol, ncol](int iel){
    auto dest_cel = elvec[iel];
    //we skip boundary els
    if(dest_cel->Dimension()!=dest_mesh->Dimension()){
      return;
    }
    const int64_t dest_gel_index = dest_cel->ReferenceIndex();
    auto periodic_it = periodic_els.find(dest_gel_index);
    if(periodic_it==periodic_els.end()){
      std::cout<<__PRETTY_FUNCTION__
               <<"\nCould not find periodic pair of element "
               <<dest_gel_index<<std::endl;
      DebugStop();
    }
    const int64_t src_gel_index = periodic_it->second;
    auto src_cel = gmesh->ElementVec()[src_gel_index]->Reference();
    if(!src_cel){
      DebugStop();
    }
    if(dest_cel->NConnects()!=src_cel->NConnects()){
      DebugStop();
    }
    const int ncon = dest_cel->NConnects();
    for(int icon = 0; icon < ncon; icon++){
      const auto src_seqnum = src_cel->Connect(icon).SequenceNumber();
      const auto src_pos = src_block.Position(src_seqnum);
      const auto src_sz = src_block.Size(src_seqnum);
      
      const auto dest_seqnum = dest_cel->Connect(icon).SequenceNumber();
      const auto dest_pos = dest_block.Position(dest_seqnum);
      const auto dest_sz = dest_block.Size(dest_seqnum);
      
      if(src_sz != dest_sz){
        PZError<<__PRETTY_FUNCTION__<<'\n'
               <<"src sz: "<<src_sz<<" dest sz: "<<dest_sz<<'\n'
               <<"src_cel index: "<<src_cel->Index()<<'\t'
               <<"dest_cel index: "<<dest_cel->Index()<<'\t'
               <<"icon: "<<icon<<std::endl;
        DebugStop();
      }
      const auto neq = src_cel->Connect(icon).NDof(src_mesh);
      if(neq != src_sz){
        PZError<<__PRETTY_FUNCTION__<<'\n'
               <<"src sz: "<<src_sz<<" src_cel->Connect(icon).NDof: "<<neq<<'\n'
               <<"src_cel index: "<<src_cel->Index()<<'\t'
               <<"dest_cel index: "<<dest_cel->Index()<<'\t'
               <<"icon: "<<icon<<std::endl;
        DebugStop();
      }
      for(int icol = 0; icol < ncol; icol++){
        for(int ieq = 0; ieq < neq; ieq++){
          const auto val = src_sol.Get(src_pos+ieq,icol);
          dest_sol.Put(dest_pos+ieq,icol,val);
        }
      }
    }
  });
}

STATE ComputeProjAndError(TPZAutoPointer<TPZCompMesh> proj_mesh,
                          TPZAutoPointer<TPZCompMesh> error_mesh,
                          TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          const TPZFMatrix<CSTATE> &projected_modes,
                          const int nm,
                          const std::set<int64_t> &bound_connects,
                          TPZAutoPointer<TPZVTKGenerator> vtk_proj,
                          TPZAutoPointer<TPZVTKGenerator> vtk_error,
                          const std::string &name,
                          const SimData &simdata
                          ){
  error_mesh->LoadSolution(projected_modes);
  //now we restrict the proj mesh
  int64_t indep_con =
    wgma::cmeshtools::RestrictDofs(proj_mesh, error_mesh, nm, bound_connects);
      
  { //we want just 1 column
    const int neq = error_mesh->Solution().Rows();
    error_mesh->Solution().Resize(neq, 1);
  }
  //now we compute projection
  constexpr bool sym{false};
  auto proj_an = wgma::scattering::Analysis(proj_mesh, simdata.n_threads,
                                            simdata.optimize_bandwidth,
                                            simdata.filter_bnd_eqs,
                                            sym);
  //now we load desired solution into proj_mesh so we can project it
  TPZFMatrix<CSTATE> &sol_proj = proj_mesh->Solution();
  wgma::cmeshtools::ExtractSolFromMesh(proj_mesh, scatt_mesh, sol_proj);
  //get reference solution (just copy, not reference)
  TPZFMatrix<CSTATE> sol_ref = error_mesh->Solution();
  wgma::cmeshtools::ExtractSolFromMesh(error_mesh, scatt_mesh, sol_ref);

  STATE norm_sol_wpbc{1};
  {
    error_mesh->LoadSolution(sol_ref);
    auto normsol =
      wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
    normsol.SetNThreads(simdata.n_threads);
    norm_sol_wpbc = std::real(normsol.ComputeNorm()[0]);
    constexpr STATE tol = 1000*std::numeric_limits<STATE>::epsilon();
    if(norm_sol_wpbc < tol){norm_sol_wpbc=1.;}
    std::cout<<"norm sol wpbc: "<<norm_sol_wpbc<<std::endl;
  }
  proj_an.Assemble();
  proj_an.Solve();
      
  //now we copy it to the error mesh
  TPZFMatrix<CSTATE> near_proj_error = sol_ref;
  wgma::cmeshtools::ExtractSolFromMesh(error_mesh, proj_mesh, near_proj_error);

      
  near_proj_error -= sol_ref;
  error_mesh->LoadSolution(near_proj_error);
  if(simdata.export_vtk_error){
    vtk_proj->Do();
    vtk_error->Do();
  }
  auto normsol =
    wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
    
  normsol.SetNThreads(simdata.n_threads);
  const auto error = std::real(normsol.ComputeNorm()[0]);
  const auto rel_error = error/norm_sol_wpbc;
  std::cout<<"nmodes ("<<name<<") "<<nm
           <<" error "<<error<<" rel error "<<rel_error<<std::endl;

  //now we remove restrictions
  wgma::cmeshtools::RemovePeriodicity(proj_mesh);
  proj_mesh->ComputeNodElCon();
  proj_mesh->CleanUpUnconnectedNodes();
  return rel_error;
}

void SolveWithPML(TPZAutoPointer<TPZCompMesh> scatt_cmesh,
                  wgma::wganalysis::Wgma2D& src_an,
                  const TPZVec<CSTATE> &src_coeffs,
                  const SimData &simdata)
{
  constexpr bool sym{false};
  auto scatt_an = wgma::scattering::Analysis(scatt_cmesh, simdata.n_threads,
                                             simdata.optimize_bandwidth,
                                             simdata.filter_bnd_eqs,
                                             sym);
    
  auto src_mesh = src_an.GetMesh();
  //get id of source materials
  wgma::scattering::SourceWgma src;
  for(auto [id,mat] : src_mesh->MaterialVec()){
    src.id.insert(id);
  }
  src.modal_cmesh = src_mesh;

  TPZFMatrix<CSTATE> sol = scatt_cmesh->Solution();
  {
    const int neq = sol.Rows();
    sol.Resize(neq,1);
    sol.Zero();
  }
  //we will add the components of the solution one by one
  const int nsol = src_coeffs.size();
  bool first_assemble{true};
  TPZSimpleTimer timer("PML:Assemble+solve",true);
  for(int isol = 0; isol < nsol; isol++){
    if(src_coeffs[isol]==(CSTATE)0){continue;}
    src_an.LoadSolution(isol);
    auto beta = std::sqrt(-src_an.GetEigenvalues()[isol]);
    wgma::scattering::LoadSource2D(scatt_cmesh, src, src_coeffs[isol]);
    wgma::scattering::SetPropagationConstant(scatt_cmesh, beta);
    if(first_assemble){
      {
        TPZSimpleTimer timer("Assemble",true);
        scatt_an.Assemble();
      }
      if(simdata.direct_solver){
        auto *pardiso = scatt_an.GetSolver().GetPardisoControl();
        const auto sp = sym ? SymProp::Sym : SymProp::NonSym;
        pardiso->SetMatrixType(sp, TPZPardisoSolver<CSTATE>::MProperty::EIndefinite);
        pardiso->SetMessageLevel(1);
        auto param = pardiso->GetParam();
        param[3] = 81;//tol 10^-8, CGS replaces LU
        param[4] = 2;
        param[9] = 15;
        pardiso->SetParam(param);
      }
      else{
        SetupPrecond(scatt_an, {});
      }
    }else{
      scatt_an.AssembleRhs(src.id);
    }
    scatt_an.Solve();
    if(first_assemble && simdata.direct_solver){
      auto *pardiso = scatt_an.GetSolver().GetPardisoControl();
      pardiso->SetMessageLevel(0);
      auto param = pardiso->GetParam();
      const auto peak_mem = param[14];
      std::cout<<"Pardiso used "<<peak_mem<<"kB during factorization stage"<<std::endl;
    }
    scatt_an.LoadSolution();
    TPZFMatrix<CSTATE> &curr_sol = scatt_cmesh->Solution();
    sol+=curr_sol;
    first_assemble = false;
  }
  scatt_an.LoadSolution(sol);
}

void SetupPrecond(wgma::scattering::Analysis &scatt_an,
                  const std::set<int64_t> &indep_cons) {
  TPZSimpleTimer solve("SetupPrecond", true);
  constexpr REAL tol = 5e-6;
      
  auto &solver = dynamic_cast<TPZStepSolver<CSTATE>&>(scatt_an.GetSolver());

  TPZAutoPointer<TPZMatrixSolver<CSTATE>> precond;
  {
    const auto &eqfilt = scatt_an.StructMatrix()->EquationFilter();
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
    wgma::precond::CreateZaglBlocks(scatt_cmesh,bnd_ids, eqfilt, eqgraph,eqgraphindex,
                                    indep_cons);

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

std::set<int> UpdatePhysicalDataSplittedMats(TPZAutoPointer<TPZGeoMesh> &gmesh,
                                             wgma::cmeshtools::PhysicalData& data,
                                             const std::map<int,int> &matid_map,
                                             const std::set<int> &orig_volids,
                                             const int dim){

  std::set<int> mats_not_found;
  std::set<int> mats_found;
  std::set<int> volids = orig_volids;
  //first we need to find all the volumetric materials
  for(auto [new_id,old_id] : matid_map){
    bool found{false};
    for(auto [vol_id,er,ur] :data.matinfovec){
      if(vol_id == old_id){
        found=true;
        mats_found.insert(new_id);
	volids.insert(new_id);
        data.matinfovec.push_back({new_id,er,ur});
        break;
      }
    }
  }

  //now we search for pmls
  for(auto [new_id,old_id] : matid_map){
    bool found{false};
    if(mats_found.find(new_id)!=mats_found.end()){continue;}
    for(auto &pmldata : data.pmlvec){
      if(pmldata->ids.find(old_id)!=pmldata->ids.end()){
	const auto pml_old_ids = pmldata->ids;
	pmldata->ids.insert(new_id);
	auto cart_pml = TPZAutoPointerDynamicCast<wgma::pml::cart::data>(pmldata);
	std::optional<int> neigh_mat_res;
	if(cart_pml){
	  REAL boundPosX{0}, boundPosY{0}, boundPosZ{0}, dX{0}, dY{0}, dZ{0};
	  wgma::gmeshtools::FindPMLWidth(gmesh, cart_pml->ids, cart_pml->t,
					 boundPosX, dX,
					 boundPosY, dY,
					 boundPosZ, dZ);
	  neigh_mat_res =
	    wgma::gmeshtools::FindCartPMLNeighbourMaterial(gmesh, dim, new_id, volids,
							   boundPosX,boundPosY,boundPosZ);
	}else{
	  auto cyl_pml = TPZAutoPointerDynamicCast<wgma::pml::cyl::data>(pmldata);
	  if(!cyl_pml){
	    DebugStop();
	  }
	  REAL rMin{0}, rMax{0}, boundPosZ{0}, dZ{0};
	  wgma::gmeshtools::FindPMLWidth(gmesh, cyl_pml->ids, cyl_pml->t,
					 rMin, rMax,
					 boundPosZ, dZ);
	  neigh_mat_res =
	    wgma::gmeshtools::FindCylPMLNeighbourMaterial(gmesh, dim, new_id, volids, rMin, boundPosZ);
            
	}
	//i hope we found it...
	if(neigh_mat_res.has_value()==false){
	  std::cout<<"Could not find neighbour of material "<<old_id<<" new id "<<new_id<<std::endl;
	  DebugStop();
	}
	const int neigh_id = neigh_mat_res.value();
	pmldata->neigh[new_id] = neigh_id;
	for(auto id : pml_old_ids){
	  pmldata->neigh[id] = neigh_id;
	}
	mats_found.insert(new_id);
	found=true;
	break;
      }
    }
    if(!found){ mats_not_found.insert(old_id);}
  }
  if(mats_not_found.size()){
    std::cout<<__PRETTY_FUNCTION__
             <<"\nCould not find the following mats:\n";
    for(auto m : mats_not_found){std::cout<<m<<' ';}
    std::cout<<"\nso they were skipped."<<std::endl;
  }
  return mats_found;
}
  
  
std::map<int,int> SplitMaterialsNearWpbc(const TPZAutoPointer<TPZCompMesh> &modal_mesh,
                                         std::set<int> &all_matids){
  //auxiliary function to create a unique matid
  auto FindFreeMatId = [](const std::set<int> &mat_ids) -> int{
    //sets are always sorted, so we know it is the minimum value
    const int minval = * mat_ids.begin();
    const int maxval = * (mat_ids.end()--);
    for(int i = minval; i < maxval; i++){
      if(mat_ids.count(i)==0){
        return i;
      }
    }
    return maxval+1;
  };
  
  //first we get all the 2d materials from the modal analysis
  std::set<int> modalmats;
  const int modaldim = modal_mesh->Dimension();
  for(auto [id,mat] : modal_mesh->MaterialVec()){
    if(mat->Dimension() == modaldim){
      modalmats.insert(id);
    }
  }
  
  auto gmesh = modal_mesh->Reference();


  //original mat id -> new mat id
  std::map<int,int> matid_map;
  std::set<int> new_ids;
  
  
  for(auto gel_2d : gmesh->ElementVec()){
    if(gel_2d->Dimension()!=modaldim){continue;}
    const int modalid = gel_2d->MaterialId();
    if(modalmats.find(modalid)==modalmats.end()){continue;}
    //now we know it is a modal element
    const int nsides = gel_2d->NSides();
    const int nnodes = gel_2d->NCornerNodes();
    for(int is = nnodes; is < nsides; is++){
      TPZGeoElSide gelside(gel_2d,is);
      TPZGeoElSide neigh = gelside.Neighbour();
      while(neigh!=gelside){
	auto gel_3d = neigh.Element();
	if(gel_3d && gel_3d->Dimension()==modaldim+1){
	  const int matid = gel_3d->MaterialId();
	  if(new_ids.find(matid) == new_ids.end()){
	    if(matid_map.find(matid) == matid_map.end()){
	      //we need to insert it into the map
	      const int new_id = FindFreeMatId(all_matids);
	      //it is no longer free
	      all_matids.insert(new_id);
	      new_ids.insert(new_id);
	      matid_map[matid]=new_id;
	    }
	    //already in the map
	    const int new_id = matid_map.at(matid);
	    gel_3d->SetMaterialId(new_id);
	  }
	}
	neigh=neigh.Neighbour();
      }
    }
  }
  return matid_map;
}

/*group elements with id in mat_ids, forming patches*/
void CreateElementGroups(TPZCompMesh *cmesh,const std::set<int> &mat_ids){

  cmesh->LoadReferences();
  auto gmesh = cmesh->Reference();
  std::set<int64_t> already_grouped;
  int ngroups{0};
  int biggest_group{-1};
  for(auto el : cmesh->ElementVec()){
    auto grp = dynamic_cast<TPZElementGroup*>(el);
    if(!el || !el->HasMaterial(mat_ids) || grp){continue;}
    const auto gel = el->Reference();
    const auto matid = gel->MaterialId();
    const auto first_edge = gel->NCornerNodes();
    const auto last_edge = first_edge+gel->NSides(1);
    //we dont use set to avoid dynamic mem alloc
    //so we must remove duplicates afterwards
    TPZManVector<TPZCompEl*,200> group_candidate;
    //every face neighbour is also an edge neighbour
    for(auto edge = first_edge; edge < last_edge; edge++){
      TPZGeoElSide gelside(gel,edge);
      TPZGeoElSide neighside = gelside.Neighbour();
      while(neighside!=gelside){
        auto gel_neigh = neighside.Element();          
        if(gel_neigh && gel_neigh->MaterialId() == matid){
          auto neigh = gel_neigh->Reference();
          if(neigh && neigh->Mesh() == cmesh){
            const auto n_index = neigh->Index();
            if(already_grouped.find(n_index)==already_grouped.end()){
              //maybe we are dealing with condensed elements?
              auto real_neigh = cmesh->ElementVec()[n_index];
              group_candidate.push_back(real_neigh);
            }
          }
        }
        neighside=neighside.Neighbour();
      }
    }
    group_candidate.push_back(el);
    //sort and remove duplicates
    RemoveDuplicates(group_candidate);
    const auto nel_in_group = group_candidate.size();
    if(nel_in_group>1){
      auto elgroup = new TPZElementGroup(*cmesh);
      for(auto elg : group_candidate){
        elgroup->AddElement(elg);
        already_grouped.insert(elg->Index());
      }
      ngroups++;
      if(biggest_group < nel_in_group ){
        biggest_group = nel_in_group;
      }
    }
  }
  
  if(ngroups){
    std::cout<<"Created "<<ngroups<<" groups, biggest one: "<<biggest_group<<std::endl;
  }
  cmesh->ExpandSolution();
}
