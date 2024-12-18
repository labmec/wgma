/**
   arrow.cpp

   This target is used to apply the waveguide port bc in a two-port device,
   it is first being tested with the arrow waveguide (antiresonant reflecting
   optical waveguide)
***/

#include "wganalysis.hpp"          // for Wgma2D, CMeshWgma2D2D
#include "scattering.hpp"          // for CMeshScattering3D
#include "gmeshtools.hpp"
#include "post/waveguidecoupling.hpp"
#include "post/orthowgsol.hpp"
#include "util.hpp"                // for CreatePath, ExtractPath
#include <json_util.hpp>

#include <SPZPeriodicData.h>
#include <TPZSimpleTimer.h>
#include <TPZVTKGenerator.h>
#include <TPZLinearEigenSolver.h>
#include <Electromagnetics/TPZWgma.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzvec.h>
#include <pzfmatrix.h>
#include <tpzautopointer.h>
#include <pzlog.h>                 // for TPZLogger
#include <thread>
#include <regex>                   // for regex_search, match_results<>::_Un...


using namespace std::complex_literals;
//!minimum shared sim data
struct SimData{
  //!.msh mesh
  std::string meshfile;
  //!wavelength
  STATE lambda{4.0};
  //!geometric scaling (floating point precision)
  REAL scale{1};
  //!materials used in the 3D analysis
  TPZVec<std::string> mats_3d;
  //!materials used in the ingoing waveguide port
  TPZVec<std::string> mats_port_in;
  //!materials used in the outgoing waveguide port
  TPZVec<std::string> mats_port_out;
  //!materials used in the probe associated with the ingoing waveguide port
  TPZVec<std::string> mats_probe_in;
  //!materials used in the probe associated with the outgoing waveguide port
  TPZVec<std::string> mats_probe_out;
  //!whether the first two modes of the input port should be in the x and y direction (plane wave)
  bool planewave_in;
  //!whether the first two modes of the output port should be in the x and y direction (plane wave)
  bool planewave_out;
  //!map of refractive indices
  std::map<std::string,CSTATE> refractive_indices;
  //!map of domain regions and number of directional refinement steps
  std::map<std::string,int> refine_regions;
  //!polynomial order
  int porder{-1};
  //!pml attenuation constant in radial direction
  CSTATE alphaPMLr{0};
  //!pml attenuation constant in x-direction
  CSTATE alphaPMLx{0};
  //!pml attenuation constant in y-direction
  CSTATE alphaPMLy{0};
  //!number of eigenvalues computed on left port
  int n_eigenpairs_left;
  //!number of eigenvalues computed on right port
  int n_eigenpairs_right;
  //!perform validation test on waveguide cross section (no discontinuity)
  bool check_mode_propagation;
  //!whether to compute reflectivity/transmittivity
  bool compute_reflection_norm;
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


//!needed data from modal analysis to create waveguide port bc
struct WpbcData{
  TPZAutoPointer<TPZCompMesh> cmesh;
  TPZFMatrix<CSTATE> wgbc_k;
  TPZVec<CSTATE> wgbc_f;
};

//! Reads sim data from file
SimData ReadSimData(const std::string &dataname);
//! Compute modal analysis for a waveguide port
TPZAutoPointer<wgma::wganalysis::Wgma2D>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const SimData& simdata,
  const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &el_map,
  const TPZVec<std::string> &mats,
  const bool planewave,
  int &nEigenpairs,
  const TPZEigenSort sortingRule,
  bool usingSLEPC,
  const std::string &suffix);

std::map<int,int> SplitMaterialsNearWpbc(const TPZAutoPointer<TPZCompMesh> &modal_mesh,
                                         std::set<int> &all_matids);

void
SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
                TPZAutoPointer<wgma::wganalysis::Wgma2D> &src_an,
                TPZAutoPointer<wgma::wganalysis::Wgma2D> &match_an,
                const TPZVec<std::map<std::string, int>> &gmshmats,
                const std::map<int,int> &split_mats,
                const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &periodic_els,
                const SimData &simdata);
void
SolveModePropagation(TPZAutoPointer<TPZGeoMesh> gmesh,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> &src_an,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> &match_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const std::map<int,int> &split_mats,
                     const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &periodic_els,
                     const SimData &simdata);


int main(int argc, char *argv[]) {

  wgma::wganalysis::using_tbb_mat=true;
  wgma::scattering::using_tbb_mat=true;
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif
  if(argc<2){
    PZError<<"Unexpected number of parameters. USAGE: ./arrow param_file"<<std::endl;
    return -1;
  }

  const std::string dataname = argv[1];
  SimData simdata = ReadSimData(dataname);
  
  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(simdata.prefix));
  
  /******************
   * eigensolver options *
   ******************/

  // how to sort eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetMagnitude};
  constexpr bool usingSLEPC {true};

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
  TPZAutoPointer<TPZGeoMesh> gmesh{nullptr};
  TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> periodic_els;
  {
    TPZSimpleTimer timer("ReadMesh",true);
    TPZAutoPointer<SPZPeriodicData> periodic_data{nullptr};
    gmesh = wgma::gmeshtools::ReadPeriodicGmshMesh(simdata.meshfile, simdata.scale,
                                                   gmshmats, periodic_data,
                                                   verbosity_lvl);
    TPZVec<std::pair<int,int>> desired_mats;
    const auto np = periodic_data->dep_mat_ids.size();
    for(auto i = 0; i < np; i++){
      const auto dep = periodic_data->dep_mat_ids[i];
      const auto indep = periodic_data->indep_mat_ids[i];
      desired_mats.push_back({dep,indep});
    }
    wgma::gmeshtools::GetPeriodicElements(gmesh.operator->(),
                                          desired_mats,
                                          periodic_data,
                                          periodic_els);
  }

  //this might break periodicity!
  for(auto [name,nref] : simdata.refine_regions){
    //we need to find the material
    int matid{-1};
    bool found{false};
    for(const auto &mats : gmshmats){
      if(found) break;
      for(const auto &[matname,id] : mats){
        if(found) break;
        if (matname == name){
          matid = id;
          found = true;
        }
      }
    }
    if(found){
      std::cout<<"Refining around "<<name<<" with id "<<matid
               <<" "<<nref<<" times "<<std::endl;
      wgma::gmeshtools::DirectionalRefinement(gmesh,{matid},nref);
    }else{
      std::cout<<"Could not find refinement target "<<name<<std::endl;
    }
  }
  // print wgma_gmesh to .txt and .vtk format
  if (simdata.print_gmesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = simdata.prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }

  //in port modal analysis
  TPZAutoPointer<wgma::wganalysis::Wgma2D>
    modal_an_in{nullptr};
  {
    modal_an_in =  ComputeModalAnalysis(gmesh,gmshmats,simdata, periodic_els,
                                        simdata.mats_port_in,
                                        simdata.planewave_in,
                                        simdata.n_eigenpairs_left,sortingRule,
                                        usingSLEPC,"_port_in");
  }

  //out port modal analysis
  TPZAutoPointer<wgma::wganalysis::Wgma2D>
    modal_an_out{nullptr};
  if(simdata.n_eigenpairs_right){
    modal_an_out = ComputeModalAnalysis(gmesh,gmshmats,simdata, periodic_els,
                                        simdata.mats_port_out,
                                        simdata.planewave_out,
                                        simdata.n_eigenpairs_right,sortingRule,
                                        usingSLEPC,"_port_out");
  }

  std::set<int> all_matids;
  for(auto &mats : gmshmats){//dim
    for(auto &[name,id] : mats){//name,matid
      all_matids.insert(id);
    }
  }
  
  auto modal_l_map =
    SplitMaterialsNearWpbc(modal_an_in->GetMesh(),all_matids);
  std::map<int,int> modal_r_map;
    if(modal_an_out){
      modal_r_map = 
        SplitMaterialsNearWpbc(modal_an_out->GetMesh(),all_matids);
    }

  //now we combine the maps but inverting key->value, so we have new_mat->old_mat
  std::map<int,int> split_mats;
  for(auto [old_mat,new_mat] : modal_l_map){
    split_mats[new_mat] = old_mat;
  }
  for(auto [old_mat,new_mat] : modal_r_map){
    split_mats[new_mat] = old_mat;
  }

  if(!simdata.check_mode_propagation){
    SolveScattering(gmesh, modal_an_in,  modal_an_out, gmshmats,
                    split_mats,periodic_els, simdata);
  }else{
    SolveModePropagation(gmesh, modal_an_in,  modal_an_out, gmshmats,
                         split_mats,periodic_els, simdata);
  }
  return 0;
}


#include "post/solutionnorm.hpp"
#include "post/wgnorm.hpp"
#include "post/waveguideportbc.hpp"
#include "post/planewave.hpp"
#include "precond.hpp"

#include "TPZParallelUtils.h"
#include <TPZMatrixWindow.h>
#include <TPZSpStructMatrix.h>
#include <TPZStructMatrixOMPorTBB.h>
#include "TPZYSMPMatrix.h"
#include <TPZNullMaterial.h>
#include <pzelementgroup.h>
#include <pzvec_extras.h>
#include <pzstepsolver.h>


wgma::cmeshtools::PhysicalData
FillDataForModalAnalysis(const TPZVec<std::map<std::string, int>> &gmshmats,
                         const SimData& simdata,
                         const TPZVec<std::string> &mats,
                         const std::string &suffix);

inline std::string
CheckForBoundary(const TPZVec<std::map<std::string, int>> &gmshmats,
                 const int dim,
                 const std::string &pattern){
  for(auto &[name,id] : gmshmats[dim]){
    const auto rx = std::regex{pattern, std::regex_constants::icase };
    const bool found = std::regex_search(name,rx);
    if(found){return name;}
  }
  return "";
}

void FindPeriodicBoundaries(const TPZVec<std::map<std::string, int>> &gmshmats,
                            const int dim,
                            const std::string &suffix,
                            const std::string &pt1,
                            const std::string &pt2,
                            std::string &s_dep,
                            std::string &s_indep);

void ComputeCouplingMat(wgma::wganalysis::Wgma2D &an,
                        std::string filename,
                        const int nthreads,
                        const bool conj);

void TransformModes(wgma::wganalysis::Wgma2D& an);

void PostProcessModes(wgma::wganalysis::Wgma2D &an,
                      std::string filename,
                      const int vtkres,
                      const int nthreads);

TPZAutoPointer<TPZCompMesh>
CreateScattMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                const TPZVec<std::map<std::string, int>> &gmshmats,
                const std::map<int,int> &split_mats,
                std::set<int> &mats_near_wpbc,
                const SimData &simdata,
                const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &el_map);

void ComputeWpbcCoeffs(wgma::wganalysis::Wgma2D& an,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       const bool positive_z, const TPZVec<CSTATE> &coeff,
                       const int nthreads);

std::set<int>
UpdatePhysicalDataSplittedMats(TPZAutoPointer<TPZGeoMesh> &gmesh,
                               wgma::cmeshtools::PhysicalData& data,
                               const std::map<int,int> &matid_map,
                               const std::set<int> &orig_volids,
                               const int dim);

void
TransferSolutionBetweenPeriodicMeshes(TPZAutoPointer<TPZCompMesh> dest_mesh,
                                      TPZAutoPointer<TPZCompMesh> src_mesh,
                                      const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &periodic_els);

STATE
ComputeProjAndError(TPZAutoPointer<TPZCompMesh> proj_mesh,
                    TPZAutoPointer<TPZCompMesh> error_mesh,
                    TPZAutoPointer<TPZCompMesh> scatt_mesh,
                    const TPZFMatrix<CSTATE> &projected_modes,
                    const int nm,
                    const std::set<int64_t> &bound_connects,
                    TPZAutoPointer<TPZVTKGenerator> vtk_proj,
                    TPZAutoPointer<TPZVTKGenerator> vtk_error,
                    const std::string &name,
                    const SimData &simdata);

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WpbcData& src_data,
                          WpbcData& match_data,
                          const TPZVec<CSTATE> &source_coeffs,
                          const int nmodes_src,
                          const int nmodes_match,
                          const std::set<int> &mats_near_wpbc,
                          const SimData &simdata,
                          int64_t &refl_pos,
                          int64_t &trans_pos);

void CreateElementGroups(TPZCompMesh *cmesh,const std::set<int> &mat_ids);

void AddWaveguidePortContribution(wgma::scattering::Analysis &scatt_an, 
                                  const int64_t indep_con_id,
                                  const int nm,
                                  const TPZFMatrix<CSTATE> &wgbc_k,
                                  const TPZVec<CSTATE> &wgbc_f);

void SetupPrecond(wgma::scattering::Analysis &scatt_an,
                  const std::set<int64_t> &indep_cons);

//! Reads sim data from file
SimData ReadSimData(const std::string &dataname){
  using json = nlohmann::json;
  std::ifstream f(dataname);
  json data = json::parse(f);
  SimData sd;
  
  //so we start by reading the materials
  std::vector<std::string> tmpvec_str;
  tmpvec_str = data["mats_3d"].get<std::vector<std::string>>();
  for(auto mat : tmpvec_str){sd.mats_3d.push_back(mat);}
  tmpvec_str = data["mats_port_in"].get<std::vector<std::string>>();
  for(auto mat : tmpvec_str){sd.mats_port_in.push_back(mat);}
  tmpvec_str = data.value("mats_probe_in",std::vector<std::string>{});
  for(auto mat : tmpvec_str){sd.mats_probe_in.push_back(mat);}
  tmpvec_str = data.value("mats_port_out",std::vector<std::string>{});
  for(auto mat : tmpvec_str){sd.mats_port_out.push_back(mat);}
  tmpvec_str = data.value("mats_probe_out",std::vector<std::string>{});
  for(auto mat : tmpvec_str){sd.mats_probe_out.push_back(mat);}

  sd.planewave_in = data.value("planewave_in",false);
  sd.planewave_out = data.value("planewave_out",false);

  
  auto tmpvec_map =
    data["refractive indices"].get<std::map<std::string,std::vector<double>>>();
  for(const auto &[name,n] : tmpvec_map){
    if(n.size() == 0 || n.size()>2){
      DebugStop();
    }
    if(n.size()==2){
      sd.refractive_indices[name] = {n[0],n[1]};
    }else{
      sd.refractive_indices[name] = {n[0],0};
    }
    std::cout<<"Read material "<<name<<" with refractive index "
             <<sd.refractive_indices[name]<<std::endl;
  }
  //now we check if every material has a refractive index
  for(auto mat : sd.mats_3d){
    if(sd.refractive_indices.count(mat)==0){
      PZError<<__PRETTY_FUNCTION__
             <<"\nCould not find refractive index of material: "<<mat<<std::endl;
      DebugStop();
    }
  }
  //including the port materials
  for(auto mat : sd.mats_port_in){
    const auto suffix_length = std::strlen("_port_in");
    const auto name = mat.substr(0,mat.length()-suffix_length);
    if(sd.refractive_indices.count(name)==0){
      PZError<<__PRETTY_FUNCTION__
             <<"\nCould not find refractive index of material: "<<mat<<std::endl;
      DebugStop();
    }
  }
  for(auto mat : sd.mats_port_out){
    const auto suffix_length = std::strlen("_port_out");
    const auto name = mat.substr(0,mat.length()-suffix_length);
    if(sd.refractive_indices.count(name)==0){
      PZError<<__PRETTY_FUNCTION__
             <<"\nCould not find refractive index of material: "<<mat<<std::endl;
      DebugStop();
    }
  }

  sd.refine_regions = data.value("refine_regions", std::map<std::string,int> {});
  
  sd.meshfile = data["meshfile"];
  sd.prefix =  data["prefix"];
  sd.lambda =  data["wavelength"];
  sd.scale = data["scale"];
  sd.couplingmat = data["export_coupling_mat"];
  sd.check_mode_propagation = data["check_mode_propagation"];
  sd.direct_solver = data.value("direct_solver",false);
  sd.compute_reflection_norm = data.value("compute_reflection_norm",false);
  std::vector<double> tmpvec_double;
  tmpvec_double = data.value("alpha_pml_x", std::vector<double>{0,0});
  if(tmpvec_double.size()!=2){
    DebugStop();
  }
  sd.alphaPMLx = {tmpvec_double[0], tmpvec_double[1]};
  tmpvec_double = data.value("alpha_pml_y", std::vector<double>{0,0});
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

TPZAutoPointer<wgma::wganalysis::Wgma2D>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const SimData& simdata,
  const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &el_map,
  const TPZVec<std::string> &mats,
  const bool planewave,
  int &nEigenpairs,
  const TPZEigenSort sortingRule,
  bool usingSLEPC,
  const std::string &suffix)
{
  auto modal_data =
    FillDataForModalAnalysis(gmshmats,simdata,mats,suffix);

  STATE max_n{0};
  for(auto [matid, er, ur] : modal_data.matinfovec){
    const auto n = std::sqrt(er);
    if(std::real(n) > max_n){max_n = std::real(n);}
  }

  
  const STATE k0scale = (2*M_PI/simdata.lambda)*simdata.scale;
  const STATE target = -1.00001*(k0scale*max_n)*(k0scale*max_n);
  const auto &p_order = simdata.porder;
  const auto &lambda = simdata.lambda;
  const auto &scale = simdata.scale;
  
  auto modal_cmesh = wgma::wganalysis::CMeshWgma2DPeriodic(gmesh,p_order,modal_data,
                                                           el_map,
                                                           lambda, scale,true);

  constexpr bool print_cmesh{false};
  if(print_cmesh){
    wgma::cmeshtools::PrintCompMesh(modal_cmesh[0], simdata.prefix+"_cmesh_mf"+suffix);
    wgma::cmeshtools::PrintCompMesh(modal_cmesh[1], simdata.prefix+"_cmesh_hc"+suffix);
    wgma::cmeshtools::PrintCompMesh(modal_cmesh[2], simdata.prefix+"_cmesh_h1"+suffix);
  }
  /******************************
   * solve(modal analysis left) *
   ******************************/
  constexpr bool verbose{false};
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
                                 simdata.optimize_bandwidth,
                                 simdata.filter_bnd_eqs);
  an->SetSolver(*solver);

  std::string modalfile{simdata.prefix+"_modal"+suffix};

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

    if(planewave){
      using namespace wgma::post;
      auto decomp = PlanewaveDecomposition<SingleSpaceIntegrator>(an->GetHCurlMesh());
      decomp.SetNThreads(simdata.n_threads);
      CSTATE v11{0}, v12{0};
      CSTATE v21{0}, v22{0};
      decomp.ComputeCoefficients(v11, v12, v21, v22);

      int64_t neq_total{0}, neq_h1{0}, neq_hcurl{0};
      an->CountActiveEqs(neq_total, neq_h1, neq_hcurl);
      TPZFMatrix<CSTATE> &ev = an->GetEigenvectors();
      //now we need to change the first two solutions

      TPZMatrixWindow<CSTATE> first_mode(ev.Elem(), neq_total,1,neq_total,neq_total);
      TPZMatrixWindow<CSTATE> second_mode(ev.Elem()+neq_total, neq_total,1,neq_total,neq_total);
      TPZFMatrix<CSTATE> first_mode_cp = first_mode;
      TPZFMatrix<CSTATE> second_mode_cp = second_mode;

      std::cout<<"v11 "<<v11<<" v12 "<<v12<<std::endl;
      std::cout<<"v21 "<<v21<<" v22 "<<v22<<std::endl;

      first_mode = v11*first_mode_cp + v12*second_mode_cp;
      second_mode = v21*first_mode_cp + v22*second_mode_cp;
      an->LoadAllSolutions();
    }

    TPZVec<CSTATE> &betavec = an->GetEigenvalues();
    nEigenpairs = betavec.size();
    std::cout<<nEigenpairs<<" eigenpairs have converged"<<std::endl;
    
    for(CSTATE &b : betavec){
      constexpr STATE tol{1e-6};
      b = std::sqrt(-b);
      if(std::abs(b.real())<tol && b.imag() > 0){
          // std::cout<<"beta transf: ";
        b=-b;
      }
    }
    
    if(simdata.export_csv_modes){
      std::ostringstream eigeninfo;
      typedef std::numeric_limits< double > dbl;
      eigeninfo.precision(dbl::max_digits10);
      for(const auto &b : betavec){
        const auto pos_sign = std::imag(b) > 0 ? "+" : "-";
        eigeninfo<<std::fixed<<std::real(b)<<pos_sign<<std::abs(std::imag(b))<<"j\n";
      }
      std::ofstream eigenfile(simdata.prefix+"_evalues"+suffix+".csv",std::ios::trunc);
      eigenfile<<eigeninfo.str();
      eigenfile.close();
    }


    constexpr bool ortho{true};
    if(ortho){
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

      
      std::string couplingfile{simdata.prefix+"_coupling"+suffix+".csv"};
      ComputeCouplingMat(*an,couplingfile,simdata.n_threads,false);
      couplingfile = simdata.prefix+"_coupling"+suffix+"_conj.csv";
      ComputeCouplingMat(*an,couplingfile,simdata.n_threads,true);
      an->LoadAllSolutions();
    }

    /*
      in the modal analysis we perform a change of variables
      now we transform back the solutions
    */
    {
      TPZSimpleTimer timer("TransformModes",true);
      TransformModes(*an);
    }
    TPZSimpleTimer timer("Normalise",true);
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


void SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> &src_an,
                     TPZAutoPointer<wgma::wganalysis::Wgma2D> &match_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const std::map<int,int> &split_mats,
                     const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &periodic_els,
                     const SimData &simdata)
{

  const bool check_probe_in = simdata.mats_probe_in.size();
  const bool check_probe_out = simdata.mats_probe_out.size();
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
  if(match_an){match_an->LoadAllSolutions();}

  //set up post processing
  TPZVec<std::string> fvars_3d = {
    "Field_real",
    "Field_imag",
    "Field_abs"};
  
  

  std::set<int> mats_near_wpbc;
  TPZAutoPointer<TPZCompMesh> scatt_mesh_wpbc =
    CreateScattMesh(gmesh,gmshmats,split_mats,mats_near_wpbc,simdata,periodic_els);

  const std::string suffix = "wpbc";
  const std::string scatt_file = simdata.prefix+"_scatt"+suffix;
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
    {
      auto h1mesh = src_an->GetH1Mesh();
      auto mfmesh = src_an->GetMesh();
      wgma::cmeshtools::RemovePeriodicity(h1mesh);
      wgma::cmeshtools::RemovePeriodicity(mfmesh);
      src_an=nullptr;
    }
    
    
    if(match_an){
      match_data.cmesh = match_an->GetHCurlMesh();
      ComputeWpbcCoeffs(match_an, match_data.wgbc_k,
                      match_data.wgbc_f,true, {},
                      simdata.n_threads);
      auto h1mesh = match_an->GetH1Mesh();
      auto mfmesh = match_an->GetMesh();
      wgma::cmeshtools::RemovePeriodicity(h1mesh);
      wgma::cmeshtools::RemovePeriodicity(mfmesh);
      match_an=nullptr;
    }
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
                                         const TPZVec<std::string> &mats,
                                         const std::string &suffix)
  {
    auto modal_data = FillDataForModalAnalysis(gmshmats, simdata, mats, suffix);
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
      const std::string error_file = simdata.prefix+"_error"+suffix;
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
  if(check_probe_in){
    CreateErrorMesh(error_mesh_left,src_data.cmesh,projected_modes_left,vtk_error_left,
                    simdata.mats_probe_in,"_probe_in");
  }
  if(check_probe_out){
    CreateErrorMesh(error_mesh_right,match_data.cmesh,projected_modes_right,vtk_error_right,
                    simdata.mats_probe_out,"_probe_out");
  }
  
  

  auto SetupProjMesh = [&simdata, &fvars_2d](TPZAutoPointer<TPZCompMesh> &proj_mesh,
                                             TPZAutoPointer<TPZCompMesh> &error_mesh,
                                             std::set<int64_t> &bound_connects,
                                             TPZAutoPointer<TPZVTKGenerator> &vtk_proj,
                                             const std::string &suffix){
    // this will be the restricted mesh close to the wg port
    proj_mesh = error_mesh->Clone();
    wgma::cmeshtools::FindDirichletConnects(proj_mesh, bound_connects);
    if(simdata.export_vtk_error){
      const std::string proj_file = simdata.prefix+"_proj"+suffix;
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

  if(check_probe_in){
    SetupProjMesh(proj_mesh_left,error_mesh_left,bound_connects_left,vtk_proj_left,"_in");
  }
  if(check_probe_out){
    SetupProjMesh(proj_mesh_right,error_mesh_right,bound_connects_right,vtk_proj_right,"_out");
  }

  std::map<std::pair<int,int>,std::pair<STATE,STATE>> wpbc_error_res;


  std::ostringstream errorinfo;
  errorinfo.precision(std::numeric_limits<STATE>::max_digits10);

  if(simdata.export_csv_error){
    std::ofstream errorfile(simdata.prefix+"_error.csv",std::ios::trunc);
    errorfile<<"nm(left),nm(right),error(left),error(right)"<<std::endl;
  }
  
  for(int im_left = 0; im_left < nmodes_left.size(); im_left++){
    const int nm_left = nmodes_left[im_left];
    for(int im_right = 0; im_right < nmodes_right.size(); im_right++){
      const int nm_right = nmodes_right[im_right];
      //eq num for obtaining reflection and transmittivity
      int64_t refl_pos{-1},trans_pos{-1};
      RestrictDofsAndSolve(scatt_mesh_wpbc, src_data, match_data,
                           src_coeffs, nm_left,nm_right,
                           mats_near_wpbc,simdata,
                           refl_pos,
                           trans_pos
                           );
      //plot
      if(simdata.export_vtk_scatt){vtk.Do();}

      //get reflection and transmission

      TPZFMatrix<CSTATE> &sol = scatt_mesh_wpbc->Solution();
      if(simdata.compute_reflection_norm){
        //for now we assume they are sequential
        const int nm = simdata.source_coeffs.size();
        std::cout<<"wavelength: "<<simdata.lambda;
        for(int i = 0; i < nm; i++){
          const auto ref = sol.GetVal(refl_pos+i,0)-src_coeffs[i];
          const auto trans =  trans_pos >= 0 ? sol.GetVal(trans_pos+i,0) : 0;
          std::cout<<" src "<<src_coeffs[i]
                   <<" ref "<<ref
                   <<" ref norm "<<std::abs(ref)
                   <<" trans "<<trans
                   <<" trans norm "<<std::abs(trans)<<std::endl;
        }
        std::string outputfile = simdata.prefix+"_reflection.csv";
        std::ofstream ost;
        ost.open(outputfile, std::ios_base::app);
        ost << std::setprecision(std::numeric_limits<STATE>::max_digits10);
        ost << simdata.lambda<<',';
        for(int i = 0; i < nm; i++){
          const CSTATE ref = sol.GetVal(refl_pos+i,0)-src_coeffs[i];
          const CSTATE trans = trans_pos >= 0 ? sol.GetVal(trans_pos+i,0) : 0;
          const char ref_sign = ref.imag() > 0 ? '+' : '-';
          const char trans_sign = trans.imag() > 0 ? '+' : '-';
          ost <<ref.real()<<ref_sign<<std::abs(ref.imag())<<'j'<<','
              <<trans.real()<<trans_sign<<std::abs(trans.imag())<<'j';
          if(i == nm - 1){ost <<std::endl;}
          else{ost << ",";}
        }
      }
      const auto error_left =
        check_probe_in ? 
        ComputeProjAndError(proj_mesh_left,error_mesh_left,scatt_mesh_wpbc,
                            projected_modes_left,nm_left,bound_connects_left,
                            vtk_proj_left,vtk_error_left,"left",
                            simdata)
        : 0;
      const auto error_right =
        check_probe_out ?
        ComputeProjAndError(proj_mesh_right,error_mesh_right,scatt_mesh_wpbc,
                            projected_modes_right,nm_right,bound_connects_right,
                            vtk_proj_right,vtk_error_right,"right",
                            simdata)
        : 0;
      wpbc_error_res.insert({{nm_left,nm_right},{error_left,error_right}});

      if(simdata.export_csv_error){
        errorinfo<<nm_left<<','<<nm_right<<','
                 <<std::fixed<<error_left<<','
                 <<std::fixed<<error_right<<'\n';
        std::ofstream errorfile(simdata.prefix+"_error.csv",std::ios::app);
        errorfile<<errorinfo.str();
        errorinfo.clear();
      }
      //removing restrictions
      wgma::cmeshtools::RemovePeriodicity(scatt_mesh_wpbc);
      scatt_mesh_wpbc->ComputeNodElCon();
      scatt_mesh_wpbc->CleanUpUnconnectedNodes();
    }
  }

  if(check_probe_in || check_probe_out){
    std::cout<<"+++++++++++++WPBC COMPARISON+++++++++++++"<<std::endl;
    for(auto [nm,error] : wpbc_error_res){
      std::cout<<"nmodes (left) "<<nm.first<<" nmodes (right) "<<nm.second
               <<" norm error (left): "<<error.first<<" norm error (right): "<<error.second <<std::endl;
    }
    std::cout<<"+++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  }
  
  wgma::cmeshtools::RemovePeriodicity(src_data.cmesh);
  if(match_data.cmesh){wgma::cmeshtools::RemovePeriodicity(match_data.cmesh);}
}

void SolveModePropagation(TPZAutoPointer<TPZGeoMesh> gmesh,
                          TPZAutoPointer<wgma::wganalysis::Wgma2D> &src_an,
                          TPZAutoPointer<wgma::wganalysis::Wgma2D> &match_an,
                          const TPZVec<std::map<std::string, int>> &gmshmats,
                          const std::map<int,int> &split_mats,
                          const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &periodic_els,
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
      auto modal_data = FillDataForModalAnalysis(gmshmats, simdata, simdata.mats_probe_in, "_probe_in");
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
        const auto beta = src_an->GetEigenvalues()[i];
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
    CreateScattMesh(gmesh,gmshmats,split_mats,mats_near_wpbc,simdata,periodic_els);


  const std::string suffix = "_wpbc";
  const std::string scatt_file = simdata.prefix+"_scatt"+suffix;
  auto vtk = TPZVTKGenerator(scatt_mesh_wpbc, fvars_3d, scatt_file, simdata.vtk_res);
  vtk.SetNThreads(simdata.n_threads);

  TPZFMatrix<CSTATE> sol_pml;
  STATE norm_sol_pml{1}, norm_error_pml{1};
  
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
      //eq num for obtaining reflection and transmittivity
      int64_t refl_pos{-1},trans_pos{-1};
      RestrictDofsAndSolve(scatt_mesh_wpbc, src_data, match_data,
                           src_coeffs, nm_left,nm_right,
                           mats_near_wpbc,simdata,
                           refl_pos, trans_pos);
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
}


/**
   AUXILIARY METHODS
**/


wgma::cmeshtools::PhysicalData
FillDataForModalAnalysis(const TPZVec<std::map<std::string, int>> &gmshmats,
                         const SimData& simdata,
                         const TPZVec<std::string> &mats,
                         const std::string &suffix)
{
  // setting up cmesh data
  const auto &alphaPMLx = simdata.alphaPMLx;
  const auto &alphaPMLy = simdata.alphaPMLy;
    
  wgma::cmeshtools::PhysicalData modal_data;
  std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
  for(const auto &matname : mats){
    //now we remove the suffix
    const auto suffix_length = suffix.size();
    const auto name = matname.substr(0,matname.length()-suffix_length);
    const CSTATE n = simdata.refractive_indices.at(name);
    modal_mats[matname] = std::pair<CSTATE,CSTATE>(n*n,1.);
  }
  std::map<std::string, wgma::bc::type> modal_bcs;
  //dimension of the modal analysis 
  constexpr int modal_dim{2};
  //first we check for periodic BCs
  constexpr int bcdim{1};
  std::string depbc, indepbc;
  FindPeriodicBoundaries(gmshmats,bcdim,suffix,"xm","xp",depbc,indepbc);
  modal_bcs[depbc] = wgma::bc::type::PERIODIC;
  modal_bcs[indepbc] = wgma::bc::type::PERIODIC;
  FindPeriodicBoundaries(gmshmats,bcdim,suffix,"ym","yp",depbc,indepbc);
  modal_bcs[depbc] = wgma::bc::type::PERIODIC;
  modal_bcs[indepbc] = wgma::bc::type::PERIODIC;
  // auto pec_bnd = CheckForBoundary(gmshmats,bcdim,"bound"+suffix);
  
  // if(pec_bnd.size() > 0){
  //   modal_bcs[pec_bnd] = wgma::bc::type::PEC;
  // }
  
  wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                          {alphaPMLx,alphaPMLy}, modal_data, modal_dim);
  //we must now filter the 2D PMLs
  std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
  for(const auto &pml : modal_data.pmlvec){
    const std::string pattern{suffix};
    const auto rx = std::regex{pattern, std::regex_constants::icase };
    
    const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
    if(found_pattern){pmlvec.push_back(pml);}
  }
  modal_data.pmlvec = pmlvec;
  return modal_data;
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

  //is transformed already
  TPZVec<CSTATE> betavec = an.GetEigenvalues();
    
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
  const auto neq_full = mesh_sol.Rows();
  const auto neq_indep = mf_mesh->NEquations();
  TPZMatrixWindow<CSTATE> reduced_sol(mesh_sol,0,0,neq_indep,nsol);
  an.SetEigenvectors(reduced_sol);
  an.LoadAllSolutions();
}

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
    "Et_abs",
    "Material"};
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

TPZAutoPointer<TPZCompMesh>
CreateScattMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                const TPZVec<std::map<std::string, int>> &gmshmats,
                const std::map<int,int> &split_mats,
                std::set<int> &mats_near_wpbc,
                const SimData &simdata,
                const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &el_map
                )
{

  const TPZVec<std::string> &mats = simdata.mats_3d;

  const CSTATE alphaPMLx = simdata.alphaPMLx;
  const CSTATE alphaPMLy = simdata.alphaPMLy;

  const auto &pOrder = simdata.porder;
  const auto &lambda = simdata.lambda;
  const auto &scale = simdata.scale;
    
  const std::string &prefix = simdata.prefix;

    
  // setting up cmesh data
  wgma::cmeshtools::PhysicalData scatt_data;
  std::map<std::string, std::pair<CSTATE, CSTATE>> scatt_mats;
  for(const auto &name : mats){
    const CSTATE n = simdata.refractive_indices.at(name);
    scatt_mats[name] = std::pair<CSTATE,CSTATE>(n*n,1.);
  }
  std::map<std::string, wgma::bc::type> scatt_bcs;

  constexpr int dim{3};
  constexpr int bcdim{2};
  std::string depbc, indepbc;
  FindPeriodicBoundaries(gmshmats,bcdim,"","xm","xp",depbc,indepbc);
  scatt_bcs[depbc] = wgma::bc::type::PERIODIC;
  scatt_bcs[indepbc] = wgma::bc::type::PERIODIC;
  FindPeriodicBoundaries(gmshmats,bcdim,"","ym","yp",depbc,indepbc);
  scatt_bcs[depbc] = wgma::bc::type::PERIODIC;
  scatt_bcs[indepbc] = wgma::bc::type::PERIODIC;
  auto pec_bnd = CheckForBoundary(gmshmats,bcdim,"bound_vol");
  
  if(pec_bnd.size() > 0){
    scatt_bcs[pec_bnd] = wgma::bc::type::PEC;
  }
    
    
  wgma::cmeshtools::SetupGmshMaterialData(gmshmats, scatt_mats, scatt_bcs,
                                          {alphaPMLx,alphaPMLy,0}, scatt_data);


  //materials that will represent our source
  std::set<int> src_ids;
  /*
    probe mats are regions of the domain in which we want to be able
    to evaluate our solution
    they are also used to ensure that the computational elements are created
    so every region that will be used in a waveguide port bc must be
    also inserted as a probe mat
  */
  std::vector<std::string> probeMats;
  {
    const std::string pattern_left{"port_in"};
    const auto rx_left =
      std::regex{pattern_left, std::regex_constants::icase };
    const std::string pattern_right{"port_out"};
    const auto rx_right =
      std::regex{pattern_right, std::regex_constants::icase };
    const std::string pattern_probe_left{"probe_in"};
    const auto rx_probe_left =
      std::regex{pattern_probe_left, std::regex_constants::icase };
    const std::string pattern_probe_right{"probe_out"};
    const auto rx_probe_right =
      std::regex{pattern_probe_right, std::regex_constants::icase };

    constexpr int dim{3};
    constexpr int probedim{dim-1};
    for(const auto &[name,id] : gmshmats[probedim]){
      const bool found_pattern =
        std::regex_search(name, rx_left) ||
        std::regex_search(name, rx_right) ||
        std::regex_search(name, rx_probe_left) ||
        std::regex_search(name, rx_probe_right);
      if(found_pattern){
        scatt_data.probevec.push_back({id,probedim});
      }
    }
  }
    
    

    
    
  //due to the waveguide port bc we need to filter out some PML regions
  {
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

  std::set<int> volids;
  for(auto [id,dummy1,dummy2] : scatt_data.matinfovec){
    volids.insert(id);
  }
    
  mats_near_wpbc =
    UpdatePhysicalDataSplittedMats(gmesh, scatt_data, split_mats,
                                   volids, dim);
  auto cmesh =
    wgma::scattering::CMeshScattering3DPeriodic(gmesh, pOrder, scatt_data,
                                                el_map,
                                                src_ids,
                                                lambda,scale,false);
  return cmesh;
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
  if(coeff.size()){
    wgbc.SetSrcCoeff(coeff);
  }
  wgbc.SetBeta(betavec);
  wgbc.ComputeContribution();
  wgbc.GetContribution(wgbc_k,wgbc_f);
}

void
TransferSolutionBetweenPeriodicMeshes(TPZAutoPointer<TPZCompMesh> dest_mesh,
                                      TPZAutoPointer<TPZCompMesh> src_mesh,
                                      const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &periodic_els)
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
    int64_t src_gel_index{-1};
    for(auto periodic_map : periodic_els){
      auto periodic_it = periodic_map->find(dest_gel_index);
      if(periodic_it!=periodic_map->end()){
        src_gel_index = periodic_it->second;
        break;
      }
    }
    if(src_gel_index < 0){
      std::cout<<__PRETTY_FUNCTION__
               <<"\nCould not find periodic pair of element "
               <<dest_gel_index<<std::endl;
      DebugStop();
    }
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

void RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WpbcData& src_data,
                          WpbcData& match_data,
                          const TPZVec<CSTATE> &source_coeffs,
                          const int nmodes_src,
                          const int nmodes_match,
                          const std::set<int> &mats_near_wpbc,
                          const SimData &simdata,
                          int64_t &refl_pos,
                          int64_t &trans_pos)
{

  auto gmesh = scatt_mesh->Reference();
  

  auto match_mesh = match_data.cmesh;
  auto src_mesh = src_data.cmesh;

  /**
     no dirichlet connects must be restricted!!
  **/
  std::set<int64_t> boundConnects;
  wgma::cmeshtools::FindDirichletConnects(scatt_mesh, boundConnects);



  int64_t indep_con_id_match = -1;
  if(match_mesh){
    indep_con_id_match =
      wgma::cmeshtools::RestrictDofs(scatt_mesh, match_mesh, nmodes_match, boundConnects);
  }
  
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
    AddWaveguidePortContribution(scatt_an, indep_con_id_src,
                                 nmodes_src, src_data.wgbc_k, src_data.wgbc_f);
    if(match_mesh){
      AddWaveguidePortContribution(scatt_an, indep_con_id_match,
                                 nmodes_match, match_data.wgbc_k, match_data.wgbc_f);
    }
  }
  
  
  if(!simdata.direct_solver){
    std::set<int64_t> indices = {indep_con_id_src};
    if(match_mesh){indices.insert(indep_con_id_match);}
    SetupPrecond(scatt_an, indices);
  }
  TPZSimpleTimer tsolve("Solve",true);
  scatt_an.Solve();

  const auto &block = scatt_mesh->Block();
  {
    const auto &indep_con = scatt_mesh->ConnectVec()[indep_con_id_src];
    const auto seqnum = indep_con.SequenceNumber();
    refl_pos = block.Position(seqnum);
  }

  if(match_mesh){
    const auto &indep_con = scatt_mesh->ConnectVec()[indep_con_id_match];
    const auto seqnum = indep_con.SequenceNumber();
    trans_pos = block.Position(seqnum);
  }else{
    trans_pos = -1;
  }
  
}

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

  
void FindPeriodicBoundaries(const TPZVec<std::map<std::string, int>> &gmshmats,
                            const int dim,
                            const std::string &suffix,
                            const std::string &pt1,
                            const std::string &pt2,
                            std::string &s_dep,
                            std::string &s_indep){

  const std::string pattern_1 = suffix+"_periodic_"+pt1;
  const auto indepname = CheckForBoundary(gmshmats, dim, pattern_1);
  if(indepname.size() > 0){
    //found
    const std::string pattern_2 = suffix+"_periodic_"+pt2;
    const auto depname = CheckForBoundary(gmshmats, dim, pattern_2);
    if(depname.size() > 0){
      //found
      s_dep = depname;
      s_indep = indepname;
    }else{
      PZError<<"ERROR:\n"
             <<"Setting up modal analysis for "<<suffix
             <<" and found periodic boundary "<<indepname
             <<" but did not find corresponding "<<pattern_2
             <<"\nAborting.."<<std::endl;
      DebugStop();
    }
  }
}