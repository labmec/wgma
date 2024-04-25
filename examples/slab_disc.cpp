/**
slab_disc.cpp

This target is used to validate the waveguide port boundary condition
using a slab waveguide as an example.
It consists of a "shrinked" version on slab_disc.cpp.
***/


#include "post/orthowgsol.hpp"     // for OrthoWgSol
#include "post/wgnorm.hpp"         // for WgNorm
#include <post/solutionnorm.hpp>   // for SolutionNorm
#include <post/waveguideportbc.hpp>// for WaveguidePortBC
#include <materials/solutionprojection.hpp>
#include "wganalysis.hpp"          // for WgmaPlanar, CMeshWgma1D
#include "scattering.hpp"          // for CMeshScattering2D
#include "bctypes.hpp"             // for type, type::PEC
#include "cmeshtools.hpp"          // for PhysicalData, SetupGmshMaterialData
#include "gmeshtools.hpp"          // for PrintGeoMesh, ReadPeriodicGmshMesh
#include "modetypes.hpp"           // for mode, mode::TE
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
#include <pzcmesh.h>               // for TPZCompMesh
#include <pzerror.h>               // for PZError
#include <pzfmatrix.h>             // for TPZFMatrix
#include <pzgmesh.h>               // for TPZGeoMesh
#include <pzlog.h>                 // for TPZLogger

#include <regex>                   // for regex_search, match_results<>::_Un...


using namespace std::complex_literals;
//!minimum shared sim data
struct SimData{
  //!.msh mesh
  std::string meshfile;
  //!wavelength
  STATE lambda{1.5};
  //!geometric scaling (floating point precision)
  REAL scale{1};
  //!reffractive index of core
  STATE ncore{1};
  //!reffractive index of clad
  STATE nclad{1};
  //!mode currently analysed
  wgma::planarwg::mode mode{wgma::planarwg::mode::TE};
  //!polynomial order
  int porder{-1};
  //!pml attenuation constant in x-direction
  CSTATE alphaPMLx{0};
  //!pml attenuation constant in y-direction
  CSTATE alphaPMLy{0};
  //!number of eigenvalues computed on left port
  int n_eigenpairs_left;
  //!number of eigenvalues computed on right port
  int n_eigenpairs_right;
  //!if set to true, wpbc is compared against PML for validation
  bool compare_pml;
  //!pairs of mode index/coefficient to be used as source
  std::vector<std::pair<int,double>> source_coeffs;
  //!index of the number of modes to be used to restrict the dofs on the wpbc(in)
  TPZVec<int> n_modes_left;
  //!index of the number of modes to be used to restrict the dofs on the wpbc(out)
  TPZVec<int> n_modes_right;
  //!whether to filter dirichlet eqs
  bool filter_bnd_eqs{true};
  //!renumber equations
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

TPZAutoPointer<TPZCompMesh>
CreateModalAnalysisMesh(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const SimData& simdata,
  const CSTATE epsilon_clad,
  const CSTATE epsilon_core,
  const std::string &name);

TPZAutoPointer<wgma::wganalysis::WgmaPlanar>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const CSTATE epsilon_clad,
  const CSTATE epsilon_core,
  const SimData& simdata,
  const CSTATE target,
  const int nEigenpairs,
  const TPZEigenSort sortingRule,
  bool usingSLEPC,
  const std::string &name);

void
SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
                wgma::wganalysis::WgmaPlanar &src_an,
                wgma::wganalysis::WgmaPlanar &match_an,
                const TPZVec<std::map<std::string, int>> &gmshmats,
                const std::map<int64_t,int64_t> &periodic_els,
                const SimData &simdata);
void
SolveScatteringWithPML(TPZAutoPointer<TPZGeoMesh> gmesh,
                       wgma::wganalysis::WgmaPlanar &src_an,
                       wgma::wganalysis::WgmaPlanar &match_an,
                       const TPZVec<std::map<std::string, int>> &gmshmats,
                       const std::map<int64_t,int64_t> &periodic_els,
                       const SimData &simdata);

SimData GetSimData()
{

  // path for output files
  const std::string path {"res_slab_disc/"};
  // common prefix for both meshes and output files
  const std::string basisName{"slab_disc"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  
  SimData data;
  data.meshfile="meshes/slab_disc.msh";
  data.lambda = 1.55;
  /*
    Given the small dimensions of the domain, scaling it can help in
    achieving good precision. Using 1./k0 as a scale factor results in
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  data.scale = data.lambda/(2*M_PI);
  data.mode = wgma::planarwg::mode::TE;
  data.ncore = 1.55;
  data.nclad = 1.00;
  data.alphaPMLx = {28, 0.0};
  data.alphaPMLy = {42, 0.0};
  data.compare_pml=false;
  data.source_coeffs = {{0,1}};
  data.porder = 4;
  data.n_eigenpairs_left = 400;
  data.n_eigenpairs_right = 400;
  data.n_modes_left = {350};
  data.n_modes_right = {350};
  data.filter_bnd_eqs = true;
  data.print_gmesh=true;
  data.export_vtk_modes = true;
  data.export_csv_modes = true;
  data.export_csv_error = true;
  data.export_vtk_scatt = true;
  data.export_vtk_error = true;
  data.vtk_res=0;
  data.prefix = prefix;
  return std::move(data);
}

int main(int argc, char *argv[]) {

#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif
  if(argc>2){
    PZError<<"Unexpected number of parameters. USAGE: ./meta_surf_1d param_file\n"
           <<"or ./meta_surf_1d to run with hardcoded params"<<std::endl;
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
  const int nEigenpairs_left = simdata.n_eigenpairs_left;
  const int nEigenpairs_right = simdata.n_eigenpairs_right;

  const CSTATE target_left{simdata.ncore*simdata.ncore};
  const CSTATE target_right{simdata.ncore*simdata.ncore};
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
  auto gmesh = wgma::gmeshtools::ReadPeriodicGmshMesh(simdata.meshfile, simdata.scale,
                                                      gmshmats, periodic_els,
                                                      verbosity_lvl);

  // print wgma_gmesh to .txt and .vtk format
  if (simdata.print_gmesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = simdata.prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }


  const CSTATE epsilon_clad{simdata.nclad*simdata.nclad};
  const CSTATE epsilon_core{simdata.ncore*simdata.ncore};
  
  //left port modal analysis
  TPZAutoPointer<wgma::wganalysis::WgmaPlanar>
    modal_l_an{nullptr};
  {
    modal_l_an =  ComputeModalAnalysis(gmesh,gmshmats,epsilon_clad,
                                       epsilon_core,simdata,target_left,
                                       nEigenpairs_left,sortingRule,
                                       usingSLEPC,"source_left");
  }

  //right port modal analysis
  TPZAutoPointer<wgma::wganalysis::WgmaPlanar>
    modal_r_an{nullptr};
  {
    modal_r_an =  ComputeModalAnalysis(gmesh,gmshmats,epsilon_clad,
                                       epsilon_core,simdata,target_right,
                                       nEigenpairs_left,sortingRule,
                                       usingSLEPC,"source_right");
  }

  if(!simdata.compare_pml){
    SolveScattering(gmesh, modal_l_an,  modal_r_an, gmshmats, periodic_els, simdata);
  }else{
    SolveScatteringWithPML(gmesh, modal_l_an,  modal_r_an, gmshmats, periodic_els, simdata);
  }
  return 0;
}

void
PostProcessModes(wgma::wganalysis::WgmaPlanar &an,
                 std::string filename,
                 const int vtkres);

TPZAutoPointer<TPZCompMesh>
CreateScattMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                const TPZVec<std::map<std::string, int>> &gmshmats,
                const SimData &simdata, bool usingPML);

void ComputeWpbcCoeffs(wgma::wganalysis::WgmaPlanar& an,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       const bool positive_z, const TPZVec<CSTATE> &coeff,
                       const wgma::planarwg::mode mode,
                       const int nthreads);

void
TransferSolutionBetweenPeriodicMeshes(TPZAutoPointer<TPZCompMesh> dest_mesh,
                                      TPZAutoPointer<TPZCompMesh> src_mesh,
                                      const std::map<int64_t,int64_t>& periodic_els);

void SolveWithPML(TPZAutoPointer<TPZCompMesh> scatt_cmesh,
                  wgma::wganalysis::WgmaPlanar& src_an,
                  const TPZVec<CSTATE> &src_coeffs,
                  const SimData &simdata);

void
ReplaceMaterialsForProjection(TPZAutoPointer<TPZCompMesh> proj_mesh);

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WpbcData& src_data,
                          WpbcData& match_data,
                          const TPZVec<CSTATE> &source_coeffs,
                          const int nmodes_src,
                          const int nmodes_match,
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

TPZAutoPointer<TPZCompMesh>
CreateModalAnalysisMesh(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const SimData& simdata,
  const CSTATE epsilon_clad,
  const CSTATE epsilon_core,
  const std::string &name)
{
  // setting up cmesh data
  const auto &nclad = simdata.nclad;
  const auto &ncore = simdata.ncore;
  const auto &alphaPMLx = simdata.alphaPMLx;
  const auto &alphaPMLy = simdata.alphaPMLy;
  const auto &mode = simdata.mode;
  const auto &pOrder = simdata.porder;
  const auto &lambda = simdata.lambda;
  const auto &scale = simdata.scale;
    
  wgma::cmeshtools::PhysicalData modal_data;
  std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
  modal_mats[name+"_clad"] = std::pair<CSTATE, CSTATE>(epsilon_clad, 1.);
  modal_mats[name+"_core"] = std::pair<CSTATE, CSTATE>(epsilon_core, 1.);
  std::map<std::string, wgma::bc::type> modal_bcs;
  modal_bcs[name+"_bnd"] = wgma::bc::type::PEC;
  //dimension of the modal analysis 
  constexpr int modal_dim{1};
  wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                          {alphaPMLx,alphaPMLy}, modal_data, modal_dim);
  //we must now filter the 1D PMLs
  std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
  for(const auto &pml : modal_data.pmlvec){
    const std::string pattern{name+"_clad"};
    const auto rx = std::regex{pattern, std::regex_constants::icase };
    
    const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
    if(found_pattern){pmlvec.push_back(pml);}
  }
  modal_data.pmlvec = pmlvec;
  return wgma::wganalysis::CMeshWgma1D(gmesh,mode,pOrder,modal_data,
                                       lambda, scale);
}

TPZAutoPointer<wgma::wganalysis::WgmaPlanar>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const CSTATE epsilon_clad,
  const CSTATE epsilon_core,
  const SimData& simdata,
  const CSTATE target,
  const int nEigenpairs,
  const TPZEigenSort sortingRule,
  bool usingSLEPC,
  const std::string &name)
{


  auto modal_cmesh =
    CreateModalAnalysisMesh(gmesh,gmshmats,simdata,
                            epsilon_clad,epsilon_core,
                            name);
  /******************************
   * solve(modal analysis left) *
   ******************************/
  constexpr bool verbose{false};
  auto solver = wgma::wganalysis::SetupSolver(target, nEigenpairs, sortingRule, usingSLEPC,nEigenpairs*10,verbose);
  if(sortingRule==TPZEigenSort::UserDefined){
    solver->SetUserSortingFunc([](CSTATE a, CSTATE b)->bool{
      const auto sqrt_a_im = std::fabs(std::imag(std::sqrt(a)));
      const auto sqrt_b_im = std::fabs(std::imag(std::sqrt(b)));
      return sqrt_a_im < sqrt_b_im;
    });
  }

  TPZAutoPointer<wgma::wganalysis::WgmaPlanar> an =
    new wgma::wganalysis::WgmaPlanar(modal_cmesh, simdata.n_threads,
                                     simdata.optimize_bandwidth,
                                     simdata.filter_bnd_eqs);
  an->SetSolver(*solver);

  std::string modalfile{simdata.prefix+"_modal_"+name};

  {
    TPZSimpleTimer analysis("Modal analysis");
    an->Assemble();
    static constexpr bool computeVectors{true};
  

    an->Solve(computeVectors,verbose);
  
    //load all obtained modes into the mesh
    an->LoadAllSolutions();

    {
      TPZSimpleTimer timer("Ortho",true);
      constexpr STATE tol{1e-14};
      constexpr bool conj{false};
      const int n_ortho = wgma::post::OrthoWgSol(an,tol,conj);
      std::cout<<"orthogonalised  "<<n_ortho<<" eigenvectors"<<std::endl;
    }


    //now we normalise them
    auto cmesh = an->GetMesh();
    //leave empty for all valid matids
    std::set<int> matids {};
    constexpr bool conj{true};
    auto norm =
      wgma::post::WgNorm<wgma::post::SingleSpaceIntegrator>(cmesh,matids,
                                                            conj,0);
    norm.SetNThreads(simdata.n_threads);
    TPZVec<CSTATE> betavec = an->GetEigenvalues();
    for(auto &b : betavec){b = sqrt(b);}
    norm.SetBeta(betavec);
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
    norm.SetWavelength(simdata.lambda/simdata.scale);
    norm.Normalise();
    TPZFMatrix<CSTATE> &mesh_sol=cmesh->Solution();
    //we update analysis object
    an->SetEigenvectors(mesh_sol);
    an->LoadAllSolutions();
  }
  
  if(simdata.export_vtk_modes){
    PostProcessModes(*an, modalfile, simdata.vtk_res);
  }

  return an;
}

//! Reads sim data from file
SimData ReadSimData(const std::string &dataname){
  using json = nlohmann::json;
  std::ifstream f(dataname);
  json data = json::parse(f);
  SimData sd;
  sd.meshfile = data["meshfile"];
  sd.prefix =  data["prefix"];
  sd.lambda =  data["wavelength"];
  sd.scale = data["scale"];
  sd.mode = wgma::planarwg::mode::TE;
  sd.ncore = data["ncore"];
  sd.nclad = data["nclad"];
  sd.compare_pml = data["compare_pml"];
  std::vector<double> tmpvec_double;
  tmpvec_double = data["alpha_pml_x"].get<std::vector<double>>();
  if(tmpvec_double.size()!=2){
    DebugStop();
  }
  sd.alphaPMLx = {tmpvec_double[0], tmpvec_double[1]};
  tmpvec_double = data["alpha_pml_y"].get<std::vector<double>>();
  if(tmpvec_double.size()!=2){
    DebugStop();
  }
  sd.alphaPMLy = {tmpvec_double[0], tmpvec_double[1]};
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
  sd.print_gmesh=data.value("print_gmes",true);
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
                     wgma::wganalysis::WgmaPlanar &src_an,
                     wgma::wganalysis::WgmaPlanar &match_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const std::map<int64_t,int64_t> &periodic_els,
                     const SimData &simdata)
{

  /*********************
   * solve(scattering) *  
   *********************/  
  TPZSimpleTimer tscatt("Scattering");

  TPZAutoPointer<TPZCompMesh> scatt_mesh_wpbc =
    CreateScattMesh(gmesh,gmshmats,simdata,false);

  const auto n_eigenpairs_left = src_an.GetEigenvalues().size();
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
  
  //set up post processing
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_imag",
    "Field_abs"};

  const std::string suffix = "wpbc";
  const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
  auto vtk = TPZVTKGenerator(scatt_mesh_wpbc, fvars, scatt_file, simdata.vtk_res);
  
  
  src_an.LoadAllSolutions();
  match_an.LoadAllSolutions();

  

  //compute wgbc coefficients
  WpbcData src_data;
  src_data.cmesh = src_an.GetMesh();
  ComputeWpbcCoeffs(src_an,  src_data.wgbc_k,
                    src_data.wgbc_f, false, src_coeffs,
                    simdata.mode,simdata.n_threads);

  WpbcData match_data;
  match_data.cmesh = match_an.GetMesh();
  ComputeWpbcCoeffs(match_an, match_data.wgbc_k,
                    match_data.wgbc_f,true, {},
                    simdata.mode,simdata.n_threads);
    


  

  auto CreateErrorMesh = [&simdata, &gmesh, &gmshmats, &fvars,
                          &periodic_els](TPZAutoPointer<TPZCompMesh> &error_mesh,
                                         TPZFMatrix<CSTATE> &projected_modes,
                                         wgma::wganalysis::WgmaPlanar &modal_an,
                                         TPZAutoPointer<TPZVTKGenerator> &vtk_error,
                                         const std::string &name)
  {
    const auto &nclad = simdata.nclad;
    const auto &ncore = simdata.ncore;
    
    error_mesh = CreateModalAnalysisMesh(gmesh, gmshmats, simdata, nclad*nclad, ncore*ncore, name);
    //now we transfer the modal solution from the WPBC to the error mesh and store it
    TransferSolutionBetweenPeriodicMeshes(error_mesh, modal_an.GetMesh(), periodic_els);
    projected_modes = error_mesh->Solution();
    const int neq = error_mesh->Solution().Rows();
    error_mesh->Solution().Resize(neq, 1);
    if(simdata.export_vtk_error){
      const std::string error_file = simdata.prefix+"_error_"+name;
      vtk_error = new TPZVTKGenerator(error_mesh, fvars, error_file, simdata.vtk_res);
    }
  };


  /*
    we want to analyse the solution close to the ports and see if the
    number of modes is sufficient to represent it.
   */
  TPZAutoPointer<TPZCompMesh> error_mesh_left{nullptr}, error_mesh_right{nullptr};
  TPZFMatrix<CSTATE> projected_modes_left, projected_modes_right;
  TPZAutoPointer<TPZVTKGenerator> vtk_error_left{nullptr}, vtk_error_right{nullptr};
  CreateErrorMesh(error_mesh_left,projected_modes_left,src_an,vtk_error_left,"eval_left");
  CreateErrorMesh(error_mesh_right,projected_modes_right,match_an,vtk_error_right,"eval_right");
  
  

  auto SetupProjMesh = [&simdata, &fvars](TPZAutoPointer<TPZCompMesh> &proj_mesh,
                                          TPZAutoPointer<TPZCompMesh> &error_mesh,
                                          std::set<int64_t> &bound_connects,
                                          TPZAutoPointer<TPZVTKGenerator> &vtk_proj,
                                          const std::string &name){
    // this will be the restricted mesh close to the wg port
    proj_mesh = error_mesh->Clone();
    ReplaceMaterialsForProjection(proj_mesh);
    wgma::cmeshtools::FindDirichletConnects(proj_mesh, bound_connects);
    if(simdata.export_vtk_error){
      const std::string proj_file = simdata.prefix+"_proj_"+name;
      vtk_proj = new TPZVTKGenerator(proj_mesh, {"Solution"}, proj_file, simdata.vtk_res);
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
                           src_coeffs, nm_left,nm_right,simdata);
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

void SolveScatteringWithPML(TPZAutoPointer<TPZGeoMesh> gmesh,
                            wgma::wganalysis::WgmaPlanar &src_an,
                            wgma::wganalysis::WgmaPlanar &match_an,
                            const TPZVec<std::map<std::string, int>> &gmshmats,
                            const std::map<int64_t,int64_t> &periodic_els,
                            const SimData &simdata)
{
  /*********************
   * solve(scattering) *  
   *********************/  
  TPZSimpleTimer tscatt("Scattering");

  const auto n_eigenpairs_left = src_an.GetEigenvalues().size();
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

  //set up post processing vars
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_imag",
    "Field_abs"};
  
  //ok first we solve with a PML
  TPZAutoPointer<TPZCompMesh> scatt_mesh_pml =
    CreateScattMesh(gmesh,gmshmats,simdata,true);
  {
    const std::string suffix = "pml";
    const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
    SolveWithPML(scatt_mesh_pml,src_an,src_coeffs,simdata);
    if(simdata.export_vtk_scatt){
      auto vtk = TPZVTKGenerator(scatt_mesh_pml, fvars, scatt_file, simdata.vtk_res);
      vtk.Do();
    }
  }

  TPZAutoPointer<TPZCompMesh> scatt_mesh_wpbc =
    CreateScattMesh(gmesh,gmshmats,simdata,false);


  const std::string suffix = "wpbc";
  const std::string scatt_file = simdata.prefix+"_scatt_"+suffix;
  auto vtk = TPZVTKGenerator(scatt_mesh_wpbc, fvars, scatt_file, simdata.vtk_res);
  
  
  src_an.LoadAllSolutions();
  match_an.LoadAllSolutions();

  

  //compute wgbc coefficients
  WpbcData src_data;
  src_data.cmesh = src_an.GetMesh();
  ComputeWpbcCoeffs(src_an,  src_data.wgbc_k,
                    src_data.wgbc_f, false, src_coeffs,
                    simdata.mode,simdata.n_threads);

  WpbcData match_data;
  match_data.cmesh = match_an.GetMesh();
  ComputeWpbcCoeffs(match_an, match_data.wgbc_k,
                    match_data.wgbc_f,true, {},
                    simdata.mode,simdata.n_threads);
    


  TPZFMatrix<CSTATE> sol_pml, ref_sol;
  STATE norm_sol_pml{1}, norm_error_pml{1};
  //here we will store the error between pml approx and wgbc approx
  TPZAutoPointer<TPZCompMesh> error_mesh{nullptr};
  TPZAutoPointer<TPZVTKGenerator> vtk_error{nullptr}, vtk_wpbc_eval{nullptr}, vtk_proj_eval{nullptr};
  {
    const auto &nclad = simdata.nclad;
    const auto &ncore = simdata.ncore;
    error_mesh = CreateModalAnalysisMesh(gmesh, gmshmats, simdata, nclad*nclad, ncore*ncore, "eval_left");
    //just to set size
    sol_pml = error_mesh->Solution();

    //now error mesh will contain all the modes
    TransferSolutionBetweenPeriodicMeshes(error_mesh, src_an.GetMesh(), periodic_els);
    //let us get the distance between src and error mesh
    const STATE dist =
      error_mesh->Element(0)->Reference()->Node(0).Coord(0)-
      src_an.GetMesh()->Element(0)->Reference()->Node(0).Coord(0);
    //now we compute the expected (propagated) solution
    ref_sol.Redim(sol_pml.Rows(),1);

    
    for(int i = 0; i < src_coeffs.size();i++){
      CSTATE *ref_sol_ptr = ref_sol.Elem();
      if(src_coeffs[i] != (CSTATE)0){
        TPZFMatrix<CSTATE> &mode = error_mesh->Solution();
        const auto neq = mode.Rows();
        const auto offset = neq*i;
        const auto beta = std::sqrt(src_an.GetEigenvalues()[i]);
        const auto coeff = src_coeffs[i]*std::exp(-1i*beta*dist);
        std::cout<<"beta "<<beta<<" dist "<<dist<<std::endl;
        CSTATE *mode_ptr = mode.Elem() + offset;
        for(int ieq = 0; ieq < neq; ieq++){
          *ref_sol_ptr++ += coeff*(*mode_ptr++);
        }
      }
    }

    
    
    wgma::cmeshtools::ExtractSolFromMesh(error_mesh, scatt_mesh_pml, sol_pml);
    error_mesh->LoadSolution(sol_pml);
    auto normsol =
        wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
    normsol.SetNThreads(simdata.n_threads);
    norm_sol_pml = std::real(normsol.ComputeNorm()[0]);

    //export PML solution at error mesh
    if(simdata.export_vtk_error){
      const std::string pml_eval_file = simdata.prefix+"_pml_eval";
      TPZVTKGenerator vtk_pml_eval(error_mesh, fvars, pml_eval_file, simdata.vtk_res);
      vtk_pml_eval.Do();
    }
    

    const std::string error_file = simdata.prefix+"_pml_error";
    vtk_error = new TPZVTKGenerator(error_mesh, fvars, error_file, simdata.vtk_res);
    const std::string wpbc_eval_file = simdata.prefix+"_wpbc_eval";
    vtk_wpbc_eval = new TPZVTKGenerator(error_mesh, fvars, wpbc_eval_file, simdata.vtk_res);
    
    std::cout<<"norm sol pml: "<<norm_sol_pml<<std::endl;

    TPZFMatrix<CSTATE> diff_sol(sol_pml);
    diff_sol-=ref_sol;
    error_mesh->LoadSolution(diff_sol);
    auto norm_error =
        wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
    norm_error.SetNThreads(simdata.n_threads);
    norm_error_pml = std::real(norm_error.ComputeNorm()[0])/norm_sol_pml;
    std::cout<<"norm error pml: "<<norm_error_pml<<std::endl;

    //export PML solution at error mesh
    if(simdata.export_vtk_error){
      error_mesh->LoadSolution(ref_sol);
      const std::string proj_eval_file = simdata.prefix+"_proj_eval";
      TPZVTKGenerator vtk_proj_eval(error_mesh, fvars, proj_eval_file, simdata.vtk_res);
      vtk_proj_eval.Do();
    }
    
  }

  TPZFMatrix<CSTATE> sol_wpbc = sol_pml;//just to have the same size
  
  std::map<std::pair<int,int>,STATE> wpbc_error_res;

  std::ostringstream errorinfo;
  errorinfo.precision(std::numeric_limits<STATE>::max_digits10);
  
  for(int im_left = 0; im_left < nmodes_left.size(); im_left++){
    const int nm_left = nmodes_left[im_left];
    for(int im_right = 0; im_right < nmodes_right.size(); im_right++){
      const int nm_right = nmodes_right[im_right];
      RestrictDofsAndSolve(scatt_mesh_wpbc, src_data, match_data,
                           src_coeffs, nm_left,nm_right,simdata);
      //plot
      if(simdata.export_vtk_scatt){vtk.Do();}

      wgma::cmeshtools::ExtractSolFromMesh(error_mesh, scatt_mesh_wpbc, sol_wpbc);
      error_mesh->LoadSolution(sol_wpbc);
      if(simdata.export_vtk_error){vtk_wpbc_eval->Do();}
      sol_wpbc -= ref_sol;
      error_mesh->LoadSolution(sol_wpbc);
      if(simdata.export_vtk_error){vtk_error->Do();}
      auto normsol =
        wgma::post::SolutionNorm<wgma::post::SingleSpaceIntegrator>(error_mesh);
      normsol.SetNThreads(simdata.n_threads);
      const auto norm = std::real(normsol.ComputeNorm()[0]);
      const auto error = norm/norm_sol_pml;
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
  
  std::cout<<"+++++++++++++PML COMPARISON+++++++++++++"<<std::endl;
  for(auto [nm,error] : wpbc_error_res){
    std::cout<<"nmodes (left) "<<nm.first<<" nmodes (right) "<<nm.second
             <<" norm rel error: "<<error <<" percentage PML error: "<<(error/norm_error_pml)*100<<std::endl;
  }
  std::cout<<"+++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  
  
}

TPZAutoPointer<TPZCompMesh>
CreateScattMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                const TPZVec<std::map<std::string, int>> &gmshmats,
                const SimData &simdata, bool usingPML)
{

    const CSTATE &ncore = simdata.ncore;
    const CSTATE &nclad = simdata.nclad;
  
    const CSTATE alphaPMLx = simdata.alphaPMLx;
    const CSTATE alphaPMLy = simdata.alphaPMLy;

    const auto &mode = simdata.mode;
    const auto &pOrder = simdata.porder;
    const auto &lambda = simdata.lambda;
    const auto &scale = simdata.scale;
    
    const std::string &prefix = simdata.prefix;

    
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData scatt_data;
    std::map<std::string, std::pair<CSTATE, CSTATE>> scatt_mats;
    scatt_mats["core_left"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    scatt_mats["core_right"] = std::make_pair<CSTATE, CSTATE>(ncore*ncore, 1.);
    scatt_mats["cladding_left"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad,1.);
    scatt_mats["cladding_right"] = std::make_pair<CSTATE, CSTATE>(nclad*nclad,1.);
    std::map<std::string, wgma::bc::type> scatt_bcs;
    scatt_bcs["scatt_bnd_mid"] = wgma::bc::type::PEC;
    if(usingPML){
      scatt_bcs["scatt_bnd_right"] = wgma::bc::type::PEC;
      scatt_bcs["scatt_bnd_left"] = wgma::bc::type::PEC;
    }
    
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, scatt_mats, scatt_bcs,
                                            {alphaPMLx,alphaPMLy}, scatt_data);


    //materials that will represent our source
    std::set<int> src_ids;
    /*
      if the domain is extended, we will create a current source at
      source_core_left and source_clad_left.
      otherwise, the source is modelled as a waveguide port boundary
      condition and there is no need for source mats
     */
    if (usingPML){
      const std::string srcMat[] = {"source_left_core", "source_left_clad"};
      for(const auto &mat : srcMat){
        src_ids.insert(gmshmats[1].at(mat));
      }
      //now we add the 1d pml mats to source mats
      for(const auto &pml : scatt_data.pmlvec){
        const std::string pattern{"source_left_clad"};
        const auto rx = std::regex{pattern, std::regex_constants::icase };
    
        const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
        if(found_pattern){
          const auto matdim = 1;
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
      probeMats.push_back("source_right_clad");
      probeMats.push_back("source_right_core");
      probeMats.push_back("source_left_clad");
      probeMats.push_back("source_left_core");
      probeMats.push_back("eval_left_clad");
      probeMats.push_back("eval_left_core");
      probeMats.push_back("eval_right_clad");
      probeMats.push_back("eval_right_core");
      //now for 1d pml mats
      for(const auto &pml : scatt_data.pmlvec){
        const std::string pattern_left{"source_left_clad"};
        const auto rx_left =
          std::regex{pattern_left, std::regex_constants::icase };
        const std::string pattern_right{"source_right_clad"};
        const auto rx_right =
          std::regex{pattern_right, std::regex_constants::icase };
        const std::string pattern_eval_left{"eval_left_clad"};
        const auto rx_eval_left =
          std::regex{pattern_eval_left, std::regex_constants::icase };
        const std::string pattern_eval_right{"eval_right_clad"};
        const auto rx_eval_right =
          std::regex{pattern_eval_right, std::regex_constants::icase };
        const bool found_pattern =
          std::regex_search(*(pml->names.begin()), rx_left) ||
          std::regex_search(*(pml->names.begin()), rx_right) ||
          std::regex_search(*(pml->names.begin()), rx_eval_left) ||
          std::regex_search(*(pml->names.begin()), rx_eval_right);

        if(found_pattern){
          probeMats.push_back(*pml->names.begin());
        }
      }
      for(auto &mat : probeMats){
        const auto matdim = 1;
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
        const std::string pattern_left{"xm"};
        const auto rx_left = std::regex{pattern_left, std::regex_constants::icase };
        const std::string pattern_right{"xp"};
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

    return wgma::scattering::CMeshScattering2D(gmesh, mode, pOrder, scatt_data,src_ids,
                                               lambda,scale);
 };
void PostProcessModes(wgma::wganalysis::WgmaPlanar &an,
                      std::string filename,
                      const int vtkres){
  TPZSimpleTimer postProc("Post processing");
    
  const std::string file = filename;
  TPZVec<std::string> fvars = {
    "Field_real",
    "Field_imag",
    "Field_abs",
    // "Field_phase",
    "Deriv_real",
    "Deriv_imag",
    "Deriv_abs",
    // "Deriv_phase"
  };

  auto cmesh = an.GetMesh();
  auto vtk = TPZVTKGenerator(cmesh, fvars, file, vtkres);
  const int64_t maxval{100};
  const auto nsol = std::min(maxval,an.GetEigenvectors().Cols());
  std::cout<<"Exporting "<<nsol<<" solutions"<<std::endl;
  for(auto isol = 0; isol < nsol ; isol++){
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

  //we add to vec
  {
    CSTATE *ptr_f = &fvec.g(pos_orig,0);
    CSTATE *ptr_wgbc = wgbc_f.begin();
    for(auto imode = 0; imode < nm; imode++){
      *ptr_f++ += *ptr_wgbc++;
    }
  }
  //now we add to matrix
  for(auto imode = 0; imode < nm; imode++){
    const auto ipos = pos_filt+imode;
    for(auto jmode = 0; jmode < nm; jmode++){
      const auto jpos = pos_filt+jmode;
      const auto kval = mat->Get(ipos,jpos);
      mat->Put(ipos,jpos,kval+wgbc_k.Get(imode,jmode));
    }
  }
}

void ComputeWpbcCoeffs(wgma::wganalysis::WgmaPlanar& an,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       const bool positive_z, const TPZVec<CSTATE> &coeff,
                       const wgma::planarwg::mode mode,
                       const int nthreads){
  auto mesh = an.GetMesh();

  TPZFMatrix<CSTATE>& sol_orig = mesh->Solution();
  const int neq = sol_orig.Rows();
  const int nsol = sol_orig.Cols();
  
  wgma::post::WaveguidePortBC<wgma::post::SingleSpaceIntegrator> wgbc(mesh);
  wgbc.SetNThreads(nthreads);
  const bool is_te = mode == wgma::planarwg::mode::TE;
  TPZManVector<CSTATE,1000> betavec = an.GetEigenvalues();
  for(auto &b : betavec){b = std::sqrt(b);}
  wgbc.SetTE(is_te);
  wgbc.SetPositiveZ(positive_z);
  constexpr bool adj{false};
  wgbc.SetAdjoint(adj);
  if(coeff.size()){
    wgbc.SetSrcCoeff(coeff);
  }
  wgbc.SetBeta(betavec);
  wgbc.ComputeContribution();
  wgbc.GetContribution(wgbc_k,wgbc_f);
}
void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WpbcData& src_data,
                          WpbcData& match_data,
                          const TPZVec<CSTATE> &source_coeffs,
                          const int nmodes_src,
                          const int nmodes_match,
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

  constexpr bool sym{false};
  auto scatt_an = wgma::scattering::Analysis(scatt_mesh, simdata.n_threads,
                                             simdata.optimize_bandwidth,
                                             simdata.filter_bnd_eqs,
                                             sym);


  std::cout<<"nmodes on source boundary: "<<nmodes_src<<std::endl;
  std::cout<<"nmodes on outgoing boundary: "<<nmodes_match<<std::endl;
    
  scatt_an.Assemble();
  //now we must add the waveguide port terms
  AddWaveguidePortContribution(scatt_an, indep_con_id_match,
                               nmodes_match, match_data.wgbc_k, match_data.wgbc_f);
  AddWaveguidePortContribution(scatt_an, indep_con_id_src,
                               nmodes_src, src_data.wgbc_k, src_data.wgbc_f);
  scatt_an.Solve();
}

void
ReplaceMaterialsForProjection(TPZAutoPointer<TPZCompMesh> proj_mesh)
{
  //replace all materials for L2 projection
  TPZMaterialT<CSTATE> *last_mat{nullptr};
  for(auto [id,mat] : proj_mesh->MaterialVec()){
    auto bnd = dynamic_cast<TPZBndCond*>(mat);
    if(!bnd){
      constexpr int dim{1};
      constexpr int soldim{1};
      auto newmat = new wgma::materials::SolutionProjection<CSTATE>(id,dim,soldim);
      delete mat;
      proj_mesh->MaterialVec()[id]=newmat;
      last_mat = newmat;
    }
  }
  for(auto [id,mat] : proj_mesh->MaterialVec()){
    auto bnd = dynamic_cast<TPZBndCond*>(mat);
    if(bnd){
      TPZFMatrix<CSTATE> v1;
      TPZVec<CSTATE> v2;
      auto newmat = last_mat->CreateBC(last_mat,id,0,v1,v2);
      delete bnd;
      proj_mesh->MaterialVec()[id]=dynamic_cast<TPZMaterial*>(newmat);
    }
  }
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

  auto elvec = dest_mesh->ElementVec();
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
                  wgma::wganalysis::WgmaPlanar& src_an,
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

  TPZFMatrix<CSTATE> sol(scatt_cmesh->NEquations(),1);
  //we will add the components of the solution one by one
  const int nsol = src_coeffs.size();
  bool first_assemble{true};
  for(int isol = 0; isol < nsol; isol++){
    if(IsZero(src_coeffs[isol])){continue;}
    src_an.LoadSolution(isol);
    auto beta = std::sqrt(src_an.GetEigenvalues()[isol]);
    wgma::scattering::LoadSource1D(scatt_cmesh, src, src_coeffs[isol]);
    wgma::scattering::SetPropagationConstant(scatt_cmesh, beta);
    if(first_assemble){
      first_assemble = false;
      scatt_an.Assemble();
    }else{
      scatt_an.AssembleRhs(src.id);
    }
    scatt_an.Solve();
    scatt_an.LoadSolution();
    TPZFMatrix<CSTATE> &curr_sol = scatt_cmesh->Solution();
    sol+=curr_sol;
  }
  scatt_an.LoadSolution(sol);
}