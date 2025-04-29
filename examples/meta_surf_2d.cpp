/**
   arrow.cpp

   This target is used to apply the waveguide port bc in a two-port device,
   it is first being tested with the arrow waveguide (antiresonant reflecting
   optical waveguide)
***/

#include "wganalysis.hpp"          // for Wgma2D, CMeshWgma2D2D
#include "planewaveanalysis.hpp"
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
  //!wavelength vector
  TPZVec<STATE> wl_vec;
  //!map of refractive indices
  TPZVec<std::map<std::string,CSTATE>> ref_index_vec;
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
  //!map of refractive indices
  std::map<std::string,CSTATE> refractive_indices;
  //!map of domain regions and number of directional refinement steps
  std::map<std::string,int> refine_regions;
  //!whether curved regions are described in csv file
  bool curved_els{false};
  //!polynomial order
  int porder{-1};
  //! maximum integer used on in port (nmodes = 2*(2*k+1)(2*k+1)
  int max_k_in;
  //! maximum integer used on out port (nmodes = 2*(2*k+1)(2*k+1)
  int max_k_out;
  //! whether to use direct solver
  bool direct_solver;
  //! tolerance for the iterative solver
  REAL solver_tol;
  //!pairs of mode index/coefficient to be used as source
  std::vector<std::pair<int,double>> source_coeffs;
  //!whether to filter dirichlet eqs
  bool filter_bnd_eqs{true};
  //!renumber equations (for modal analysis, scattering is always false)
  bool optimize_bandwidth{true};
  //!output geometric mesh in .txt and .vtk files
  bool print_gmesh{false};
  //!post process modal fields
  bool export_vtk_modes{false};
  //!post process scatt fields
  bool export_vtk_scatt{false};
  //!whether eigensolver is verbose
  bool eigen_verbose{false};
  //!vtk resolution
  int vtk_res{0};
  //!number of threads
  int n_threads{(int)std::thread::hardware_concurrency()};
  //!prefix for both meshes and output files
  std::string prefix{""};
};


//!needed data to be exported from modal analysis
struct ModalData{
  TPZAutoPointer<TPZCompMesh> cmesh_mf;
  TPZAutoPointer<TPZCompMesh> cmesh_hcurl;
  TPZAutoPointer<TPZCompMesh> cmesh_h1;
  TPZVec<CSTATE> eigenvalues;
};

//!needed data from modal analysis to create waveguide port bc
struct WpbcData{
  TPZAutoPointer<TPZCompMesh> cmesh;
  TPZFMatrix<CSTATE> wgbc_k;
  TPZVec<CSTATE> wgbc_f;
};

//! Reads sim data from file
SimData ReadSimData(const std::string &dataname);
//! Refine geometric mesh
void RefineRegions(TPZAutoPointer<TPZGeoMesh> &gmesh,
                   const TPZVec<std::map<std::string, int>> &gmshmats,
                   const std::map<std::string,int> &refine_regions);

//! Replace elements for curved elements in cylinder/spherical regions
void ReplaceForCurvedEls(const std::string & meshfile,
                         TPZAutoPointer<TPZGeoMesh> &gmesh,
                         const REAL scale);
//! Compute modal analysis for a waveguide port
TPZAutoPointer<ModalData>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const SimData& simdata,
  const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &el_map,
  const TPZVec<std::string> &mats,
  const int max_k,
  const std::string &suffix);

std::map<int,int> SplitMaterialsNearWpbc(const TPZAutoPointer<TPZCompMesh> &modal_mesh,
                                         std::set<int> &all_matids);

REAL
SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
                TPZAutoPointer<ModalData> &src_an,
                TPZAutoPointer<ModalData> &match_an,
                const TPZVec<std::map<std::string, int>> &gmshmats,
                const std::map<int,int> &split_mats,
                const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &periodic_els,
                const SimData &simdata,
                TPZFMatrix<CSTATE> &last_sol);

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


    if(simdata.curved_els){
      ReplaceForCurvedEls(simdata.meshfile, gmesh, simdata.scale);
    }
  }


  //now we refine towards any given entities in refined_regions
  RefineRegions(gmesh, gmshmats, simdata.refine_regions);
  
  // print wgma_gmesh to .txt and .vtk format
  if (simdata.print_gmesh) {
    // prefix for the wgma_gmesh files
    const std::string filename = simdata.prefix + "_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh, filename);
  }


  //this map will be filled in the first iteration and
  //allows for a faster assembly of the scatt matrix
  std::map<int,int> split_mats;
  //this matrix will store the solution from the previous iteration
  TPZFMatrix<CSTATE> last_sol;
  //number of wavelength points
  const int nwl_pts = simdata.wl_vec.size();
  for(int iwl = 0; iwl < nwl_pts; iwl++){
    auto timer_begin = std::chrono::high_resolution_clock::now();
    simdata.lambda = simdata.wl_vec[iwl];
    simdata.refractive_indices = simdata.ref_index_vec[iwl];
    //in port modal analysis
    TPZAutoPointer<ModalData> modal_an_in =
      ComputeModalAnalysis(gmesh,gmshmats,simdata, periodic_els,
                           simdata.mats_port_in,
                           simdata.max_k_in,"_port_in");

    //out port modal analysis
    TPZAutoPointer<ModalData>
      modal_an_out{nullptr};
    if(simdata.mats_port_out.size()){
      modal_an_out = ComputeModalAnalysis(gmesh,gmshmats,simdata, periodic_els,
                                          simdata.mats_port_out,
                                          simdata.max_k_out,"_port_out");
    }

    

    if(iwl == 0){
      std::set<int> all_matids;
      for(auto &mats : gmshmats){//dim
        for(auto &[name,id] : mats){//name,matid
          all_matids.insert(id);
        }
      }
      //we just need to do it on the first iteration
      auto modal_l_map =
        SplitMaterialsNearWpbc(modal_an_in->cmesh_mf,all_matids);
      std::map<int,int> modal_r_map;
      if(modal_an_out){
        modal_r_map = 
          SplitMaterialsNearWpbc(modal_an_out->cmesh_mf,all_matids);
      }

      //now we combine the maps but inverting key->value, so we have new_mat->old_mat
      for(auto [old_mat,new_mat] : modal_l_map){
        split_mats[new_mat] = old_mat;
      }
      for(auto [old_mat,new_mat] : modal_r_map){
        split_mats[new_mat] = old_mat;
      }

      std::cout<<"split mats"<<std::endl;
    }
    const auto res = SolveScattering(gmesh, modal_an_in,  modal_an_out, gmshmats,
                                     split_mats,periodic_els, simdata, last_sol);
    auto timer_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = timer_end-timer_begin;
    std::cout<<"wavelength "<<simdata.lambda<<" took "<<duration.count()<<" ms"<<std::endl;
    if(res > simdata.solver_tol){
      //fall back to direct solver
      simdata.direct_solver = true;
    }
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

TPZAutoPointer<TPZCompMesh>
CreateScattMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                const TPZVec<std::map<std::string, int>> &gmshmats,
                const std::map<int,int> &split_mats,
                std::set<int> &mats_near_wpbc,
                const SimData &simdata,
                const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &el_map);

void ComputeWpbcCoeffs(ModalData & modal_data,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       const bool positive_z, const TPZVec<CSTATE> &coeff,
                       const int nthreads);

std::set<int>
UpdatePhysicalDataSplittedMats(TPZAutoPointer<TPZGeoMesh> &gmesh,
                               wgma::cmeshtools::PhysicalData& data,
                               const std::map<int,int> &matid_map,
                               const std::set<int> &orig_volids,
                               const int dim);


REAL RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WpbcData& src_data,
                          WpbcData& match_data,
                          const TPZVec<CSTATE> &source_coeffs,
                          const int nmodes_src,
                          const int nmodes_match,
                          const std::set<int> &mats_near_wpbc,
                          const SimData &simdata,
                          TPZFMatrix<CSTATE> &sol_vec,
                          int64_t &refl_pos,
                          int64_t &trans_pos);

void CreateElementGroups(TPZCompMesh *cmesh,const std::set<int> &mat_ids);

void AddWaveguidePortContribution(wgma::scattering::Analysis &scatt_an, 
                                  const int64_t indep_con_id,
                                  const int nm,
                                  const TPZFMatrix<CSTATE> &wgbc_k,
                                  const TPZVec<CSTATE> &wgbc_f);

void SetupPrecond(wgma::scattering::Analysis &scatt_an,
                  const std::set<int64_t> &indep_cons,
                  const REAL tol,
                  int from_current);

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
  tmpvec_str = data.value("mats_port_out",std::vector<std::string>{});
  for(auto mat : tmpvec_str){sd.mats_port_out.push_back(mat);}

  //we check if every port material has a corresponding 3d mat

  for(auto mat : sd.mats_port_in){
    const auto suffix_length = std::strlen("_port_in");
    const auto name = mat.substr(0,mat.length()-suffix_length);
    if ( std::find(sd.mats_3d.begin(), sd.mats_3d.end(), name) == sd.mats_3d.end() ){
      PZError<<__PRETTY_FUNCTION__
             <<"\nCould not find corresponding material of "<<mat<<std::endl;
      DebugStop();
    }
  }

  for(auto mat : sd.mats_port_out){
    const auto suffix_length = std::strlen("_port_out");
    const auto name = mat.substr(0,mat.length()-suffix_length);
    if ( std::find(sd.mats_3d.begin(), sd.mats_3d.end(), name) == sd.mats_3d.end() ){
      PZError<<__PRETTY_FUNCTION__
             <<"\nCould not find corresponding material of "<<mat<<std::endl;
      DebugStop();
    }
  }
  
  auto test_str_ad =
    data["test_str"].get<std::vector<std::tuple<double,std::map<std::string,std::vector<double>>>>>();

  for(auto [wl, matinfo] : test_str_ad){
    sd.wl_vec.push_back(wl);

    std::map<std::string,CSTATE> refractive_indices;
    
    for(const auto &[name,n] : matinfo){
      if(n.size() == 0 || n.size()>2){
        DebugStop();
      }
      if(n.size()==2){
        refractive_indices[name] = {n[0],n[1]};
      }else{
        refractive_indices[name] = {n[0],0};
      }
    }
    //now we check if every material has a refractive index
    for(auto mat : sd.mats_3d){
      if(refractive_indices.count(mat)==0){
        PZError<<__PRETTY_FUNCTION__
               <<"\nCould not find refractive index of material: "<<mat<<std::endl;
        DebugStop();
      }
    }
    sd.ref_index_vec.push_back(refractive_indices);
  }

  sd.refine_regions = data.value("refine_regions", std::map<std::string,int> {});

  sd.curved_els = data["curved_els"];
  sd.meshfile = data["meshfile"];
  sd.prefix =  data["prefix"];
  sd.scale = data["scale"];
  sd.direct_solver = data.value("direct_solver",false);
  sd.solver_tol = data.value("solver_tol",(REAL)5e-5);
  
  sd.porder = data["porder"];
  sd.source_coeffs = data["source_coeffs"].get<std::vector<std::pair<int,double>>>();
  sd.max_k_in = data["max_k_in"];
  sd.max_k_out = data["max_k_out"];
  sd.filter_bnd_eqs = data.value("filter_bnd_eqs",true);
  sd.print_gmesh=data.value("print_gmesh",true);
  sd.optimize_bandwidth = data.value("optimize_bandwidth",true);
  sd.export_vtk_modes = data.value("export_vtk_modes",false);
  sd.export_vtk_scatt = data.value("export_vtk_scatt",true);
  sd.eigen_verbose = data.value("eigen_verbose",false);
  sd.vtk_res = data.value("vtk_res",(int)0);
  sd.n_threads = data.value("n_threads",(int)std::thread::hardware_concurrency());
  return sd;
}

//! Refine geometric mesh
void
RefineRegions(TPZAutoPointer<TPZGeoMesh> &gmesh,
              const TPZVec<std::map<std::string, int>> &gmshmats,
              const std::map<std::string,int> &refine_regions)
{

  TPZVec<std::map<std::string,int>> arranged_mats(4,{});
  //this might break periodicity!
  for(auto [name,nref] : refine_regions){
    //we need to find the material
    int matid{-1};
    bool found{false};
    for(int idim = 0; idim < 4; idim++){
      const auto &mats = gmshmats[idim];
      if(found) break;
      for(const auto &[matname,id] : mats){
        if(found) break;
        if (matname == name){
          matid = id;
          arranged_mats[idim].insert({name,matid});
          found = true;
        }
      }
    }
    if(!found){
      std::cout<<"Could not find refinement target "<<name<<std::endl;
    }
  }
  for(auto allmats : arranged_mats){
    for(auto [name,matid] : allmats){
      const auto nref = refine_regions.at(name);
      //maybe we only want to set porder==0 next to it
      if(nref > 0){
        std::cout<<"Refining around "<<name<<" with id "<<matid
                 <<" "<<nref<<" times "<<std::endl;
        wgma::gmeshtools::DirectionalRefinement(gmesh,{matid},nref);
      }
    }
  }
}

//! Replace elements for curved elements in cylinder/spherical regions
void ReplaceForCurvedEls(const std::string & meshfile, TPZAutoPointer<TPZGeoMesh> &gmesh,
                         const REAL scale)
{
  //useful lambda for reading csv file
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

  TPZVec<wgma::gmeshtools::ArcData> arcs;
  TPZVec<wgma::gmeshtools::CylinderData> cyls;
  TPZVec<wgma::gmeshtools::SphereData> spheres;
  TPZVec<wgma::gmeshtools::TorusData> toruses;
  //let us remove the .msh extension
  const std::string file_prefix = meshfile.substr(0, meshfile.length() - 4);
  //first arcs
  {
    const std::string arc_suffix = "_arcdata.csv";
    const std::string arc_file = file_prefix + arc_suffix;
    std::ifstream read(arc_file);
    if (!read) {
      std::cout << "Couldn't find the arc data file " << arc_file << std::endl;
    }
    auto line = getNextLineAndSplitIntoTokens(read); // header
    line = getNextLineAndSplitIntoTokens(read);
    // we expect xc, yc, zc, r (in um), and matid
    const auto factor = 1./scale;
    while (line.size() == 5) {
      wgma::gmeshtools::ArcData arc;

      arc.m_xc = std::stod(line[0]) * factor;
      arc.m_yc = std::stod(line[1]) * factor;
      arc.m_zc = std::stod(line[2]) * factor;
      arc.m_radius = std::stod(line[3]) * factor;
      arc.m_matid = std::stoi(line[4]);
      const int narcs = arcs.size();
      arcs.Resize(narcs + 1);
      arcs[narcs] = arc;

      line = getNextLineAndSplitIntoTokens(read);
    }
  }
  //then cylinders
  {
    const std::string cyl_suffix = "_cyldata.csv";
    const std::string cyl_file = file_prefix + cyl_suffix;
    std::ifstream read(cyl_file);
    if (!read) {
      std::cout << "Couldn't find the cylinder data file " << cyl_file << std::endl;
    }
    auto line = getNextLineAndSplitIntoTokens(read); // header
    line = getNextLineAndSplitIntoTokens(read);
    // we expect xc, yc, zc, xaxis, yaxis, zaxis, r, and matid
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
      cyls.push_back(cyl);

      line = getNextLineAndSplitIntoTokens(read);
    }
  }
  
  //then spheres
  {
    const std::string sphere_suffix = "_spheredata.csv";
    const std::string sphere_file = file_prefix + sphere_suffix;
    std::ifstream read(sphere_file);
    if (!read) {
      std::cout << "Couldn't find the sphere data file " << sphere_file << std::endl;
    }
    auto line = getNextLineAndSplitIntoTokens(read); // header
    line = getNextLineAndSplitIntoTokens(read);
    // we expect xc, yc, zc, xaxis, yaxis, zaxis, r, and matid
    const auto factor = 1./scale;
    while (line.size() == 5) {
      wgma::gmeshtools::SphereData sphere;
      sphere.m_xc = std::stod(line[0]) * factor;
      sphere.m_yc = std::stod(line[1]) * factor;
      sphere.m_zc = std::stod(line[2]) * factor;
      sphere.m_radius = std::stod(line[3]) * factor;
      sphere.m_matid = std::stoi(line[4]);
      spheres.push_back(sphere);

      line = getNextLineAndSplitIntoTokens(read);
    }
  }
  //then toruses
  {
    const std::string torus_suffix = "_torusdata.csv";
    const std::string torus_file = file_prefix + torus_suffix;
    std::ifstream read(torus_file);
    if (!read) {
      std::cout << "Couldn't find the torus data file " << torus_file << std::endl;
    }
    auto line = getNextLineAndSplitIntoTokens(read); // header
    line = getNextLineAndSplitIntoTokens(read);
    // we expect xc, yc, zc, r_small, r_large (in um), and matid
    const auto factor = 1./scale;
    while (line.size() == 6) {
      wgma::gmeshtools::TorusData torus;

      torus.m_xc = std::stod(line[0]) * factor;
      torus.m_yc = std::stod(line[1]) * factor;
      torus.m_zc = std::stod(line[2]) * factor;
      torus.m_r_small = std::stod(line[3]) * factor;
      torus.m_r_large = std::stod(line[4]) * factor;
      torus.m_matid = std::stoi(line[5]);
      toruses.push_back(torus);

      line = getNextLineAndSplitIntoTokens(read);
    }
  }
  
  wgma::gmeshtools::SetExactArcRepresentation(gmesh, arcs, false);
  wgma::gmeshtools::SetExactCylinderRepresentation(gmesh, cyls, false);
  wgma::gmeshtools::SetExactSphereRepresentation(gmesh, spheres, false);
  wgma::gmeshtools::SetExactTorusRepresentation(gmesh, toruses, false);
  wgma::gmeshtools::ReplaceNeighsWithBlend(*gmesh);
}

TPZAutoPointer<ModalData>
ComputeModalAnalysis(
  TPZAutoPointer<TPZGeoMesh> gmesh,
  const TPZVec<std::map<std::string, int>> &gmshmats,
  const SimData& simdata,
  const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &el_map,
  const TPZVec<std::string> &mats,
  const int max_k,
  const std::string &suffix)
{

  TPZSimpleTimer analysis("Modal analysis",true);
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
  const bool verbose = true;//simdata.eigen_verbose;

  //now we find the coordinates of the boundaries
  REAL xMin{0},xMax{0},yMin{0},yMax{0},zMin{0},zMax{0};
  std::set<int> mat_ids;
  for(auto &matinfo : modal_data.matinfovec){
    mat_ids.insert(std::get<0>(matinfo));
  }
  wgma::gmeshtools::FindRegionLimits(gmesh, mat_ids,
                                     xMin, xMax,
                                     yMin, yMax,
                                     zMin, zMax);
  
  const STATE lx = xMax-xMin;
  const STATE ly = yMax-yMin;
  
  auto modal_cmesh = wgma::planewaveanalysis::CMeshPlaneWave2D(gmesh,p_order,modal_data,
                                                               el_map, lambda,
                                                               lx, ly, max_k,
                                                               scale,verbose);

  constexpr bool print_cmesh{false};
  if(print_cmesh){
    wgma::cmeshtools::PrintCompMesh(modal_cmesh[0], simdata.prefix+"_cmesh_mf"+suffix);
    wgma::cmeshtools::PrintCompMesh(modal_cmesh[1], simdata.prefix+"_cmesh_hc"+suffix);
    wgma::cmeshtools::PrintCompMesh(modal_cmesh[2], simdata.prefix+"_cmesh_h1"+suffix);
  }

  TPZAutoPointer<wgma::planewaveanalysis::Analysis>
    an = new wgma::planewaveanalysis::Analysis(modal_cmesh, simdata.n_threads,
                                               simdata.optimize_bandwidth,
                                               simdata.filter_bnd_eqs);
  
  std::string modalfile{simdata.prefix+"_modal"+suffix};

  an->Run();
  //load all obtained modes into the mesh
  an->LoadSolution();

  TPZVec<CSTATE> betavec = an->GetEigenvalues();


  // if(simdata.export_vtk_modes){
  //   const auto beta_sqrt = betavec;
  //   an->SetEigenvalues(betavec_cp);
  //   PostProcessModes(*an, modalfile, simdata.vtk_res,simdata.n_threads);
  //   an->LoadAllSolutions();
  //   an->SetEigenvalues(beta_sqrt);
  // }
    
  
  TPZSimpleTimer timer("Normalise");
  //now we normalise them
  auto cmesh = an->GetMesh();
  //leave empty for all valid matids
  std::set<int> matids {};
  constexpr bool conj{true};
  auto norm =
    wgma::post::WgNorm<wgma::post::MultiphysicsIntegrator,1>(cmesh,matids,
                                                             conj,simdata.n_threads);
  norm.SetNThreads(simdata.n_threads);    
  norm.SetBeta(betavec);
  norm.SetWavelength(simdata.lambda/simdata.scale);
  norm.Normalise();
  TPZFMatrix<CSTATE> &mesh_sol=cmesh->Solution();
  an->LoadSolution(mesh_sol);
  auto normvec = norm.ComputeNorm();
  
  // for(auto iev = 0; iev < betavec.size(); iev++){
  //   std::cout<<"iev "<<iev<<" beta "<<betavec[iev]<<" norm "<<normvec[iev]<<std::endl;
  // }
  //we dont need them anymore, let us free up memory
  an->GetSolver().SetMatrix(nullptr);
  TPZAutoPointer<ModalData> data = new ModalData;
  data->cmesh_h1 = an->GetH1Mesh();
  data->cmesh_hcurl = an->GetHCurlMesh();
  data->cmesh_mf = an->GetMesh();
  data->eigenvalues = an->GetEigenvalues();
  return data;
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

#include "TPZCompMeshTools.h"
#include "pzsubcmesh.h"
REAL SolveScattering(TPZAutoPointer<TPZGeoMesh> gmesh,
                     TPZAutoPointer<ModalData> &src_an,
                     TPZAutoPointer<ModalData> &match_an,
                     const TPZVec<std::map<std::string, int>> &gmshmats,
                     const std::map<int,int> &split_mats,
                     const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &periodic_els,
                     const SimData &simdata,
                     TPZFMatrix<CSTATE> &last_sol)
{
  /*********************
   * solve(scattering) *  
   *********************/  
  TPZSimpleTimer tscatt("Scattering");

  const auto nmodes_left = src_an->cmesh_mf->Solution().Cols();
  //maybe we have dirichlet on the out boundary?
  const auto nmodes_right =
    match_an ? match_an->cmesh_mf->Solution().Cols() : 0;
  /*
    the source is written as a linear combination of the modes
    this vector contains the coefficients of such combination
  */
  TPZVec<CSTATE> src_coeffs(nmodes_left,0);
  for(auto [i, alpha] : simdata.source_coeffs){
    if(i >= src_coeffs.size()){
      std::cout<<"ERROR: src coefficient bigger than computed number of modes\n"
               <<"i: "<<i<<" alpha "<<alpha<<std::endl;
      DebugStop();
    }
    src_coeffs[i] = alpha;
  }


  //set up post processing
  TPZVec<std::string> fvars_3d = {
    "Field_real",
    "Field_imag",
    "Field_abs",
    "Material"};
  
  

  std::set<int> mats_near_wpbc;
  TPZAutoPointer<TPZCompMesh> scatt_mesh_wpbc =
    CreateScattMesh(gmesh,gmshmats,split_mats,mats_near_wpbc,simdata,periodic_els);

  //now we reduce the polynomial order on refined edges
  {
    //just to avoid iterating through the same elemnet over and over
    std::set<int64_t> refined_els;
    
    for(auto [name,nref] : simdata.refine_regions){
      bool found{false};
      int matid{-1};
      for(auto idim = 0; idim < 3 && found==false; idim++){
        auto &mats_dim = gmshmats[idim];
        if(mats_dim.count(name)){
          found=true;
          matid = mats_dim.at(name);
        }
      }
      if(found==false){
        DebugStop();
      }
      //now we iterate through the mesh
      for(auto gel : gmesh->ElementVec()){
        if(!gel){continue;}
        if(gel->MaterialId() !=matid){continue;}
        //we want the most refined subelements
        if(gel->NSubElements()){continue;}
        //now we found an element that resulted from refinement, we need to see
        //its neighbouring compels
        const auto nsides = gel->NSides();
        for(auto is = 0; is < nsides-1; is++){
          TPZGeoElSide gelside(gel,is);
          TPZGeoElSide neigh = gelside.Neighbour();
          while(neigh && neigh != gelside){
            auto neigh_gel = neigh.Element();
            if(!neigh_gel){DebugStop();}
            auto cel = neigh_gel->Reference();
            if(cel){
              if(cel->Mesh() == scatt_mesh_wpbc.operator->()){
                if(refined_els.count(cel->Index())==0){
                  //ok, now we found someone that needs their p-order reduced
                  refined_els.insert(cel->Index());
                  //we need to distinguish between edge and face/interior connects
                  const int nedges = neigh_gel->NSides(1);
                  const int ncon = cel->NConnects();
                  for(auto icon = 0; icon < ncon; icon++){
                    TPZConnect &c = cel->Connect(icon);
                    if(c.Order()>0 && c.NShape() > 0){
                      if(c.HasDependency()) {DebugStop();}
                      const auto cindex = cel->ConnectIndex(icon);
                      const int64_t seq = c.SequenceNumber();
                      if(seq < 0){
                        DebugStop();
                      }
                      c.SetOrder(0,cindex);

                      //edges will have one function, faces and interior 0
                      const int nshape = icon < nedges ?  1 : 0;
                      c.SetNShape(nshape);
                      // reset the size of the block of the connect
                      scatt_mesh_wpbc->Block().Set(seq,nshape);
                    }
                  }//for connects
                }//if first time checking el
              }//if same mesh
            }//if cel
            neigh = neigh.Neighbour();
          }//while neigh
        }//for sides
      }//for gel
    }//for mat

    
    if(refined_els.size()){
      scatt_mesh_wpbc->ComputeNodElCon();
      scatt_mesh_wpbc->CleanUpUnconnectedNodes();
      scatt_mesh_wpbc->ExpandSolution();
    }


    if constexpr (0)//now we will create the subcmesh
    {
      std::cout<<"creating subcmesh with "<<refined_els.size()<<" elements"<<std::endl;
      TPZSubCompMesh * subcmesh = new TPZSubCompMesh(*scatt_mesh_wpbc);
      for(auto el_index : refined_els){
        subcmesh->TransferElement(scatt_mesh_wpbc.operator->(), el_index);
      }


      scatt_mesh_wpbc->ComputeNodElCon();
      subcmesh->MakeAllInternal();
      
      scatt_mesh_wpbc->ComputeNodElCon();
      {
        int64_t nc = scatt_mesh_wpbc->NConnects();
        for (int64_t ic = 0; ic<nc; ic++) {
            TPZConnect &c = scatt_mesh_wpbc->ConnectVec()[ic];
            if(c.NElConnected() == 0 && c.HasDependency())
            {
                c.RemoveDepend();
            }
        }
      }
      scatt_mesh_wpbc->ComputeNodElCon();
      scatt_mesh_wpbc->CleanUpUnconnectedNodes();



      subcmesh->ComputeNodElCon();
      subcmesh->CleanUpUnconnectedNodes();

      
      TPZAutoPointer<TPZCompMesh> subcmeshptr(subcmesh);
      wgma::cmeshtools::PrintCompMesh(subcmeshptr, simdata.prefix+"subcmesh.txt");
      subcmesh->Print();
      subcmesh->InitializeBlock();
      TPZCompEl *cel = dynamic_cast<TPZCompEl*>(subcmesh);
      int64_t numeq2=0;
      //??
      int64_t ic;
      for (ic=0; ic<cel->NConnects(); ic++) {
        int64_t conindex = cel->ConnectIndex(ic);
        TPZConnect &cn = cel->Mesh()->ConnectVec()[conindex];
        if (cn.SequenceNumber()<0 || cn.HasDependency()) {
          DebugStop();
        }
        int64_t seqnum = cn.SequenceNumber();
        int blsize = cel->Mesh()->Block().Size(seqnum);
        numeq2 += blsize;
      }
      std::cout<<"neq "<<numeq2
               <<" n int "<<subcmesh->NumInternalEquations()<<std::endl;
    }
    
    TPZCompMeshTools::CreatedCondensedElements(scatt_mesh_wpbc.operator->(),
                                               false, false);
  }
  const std::string suffix = "wpbc";
  const std::string scatt_file = simdata.prefix+"_scatt"+suffix;
  auto vtk = TPZVTKGenerator(scatt_mesh_wpbc, fvars_3d, scatt_file, simdata.vtk_res);
  vtk.SetNThreads(simdata.n_threads);
  

  
  //compute wgbc coefficients
  WpbcData src_data;
  WpbcData match_data;

  {
    TPZSimpleTimer timer("wpbc coeffs",true);
    src_data.cmesh = src_an->cmesh_hcurl;
    ComputeWpbcCoeffs(src_an,  src_data.wgbc_k,
                      src_data.wgbc_f, false, src_coeffs,
                      simdata.n_threads);

    //only hcurl mesh is needed from now on, we can delete the  analysis object
    {
      auto h1mesh = src_an->cmesh_h1;
      auto mfmesh = src_an->cmesh_mf;
      wgma::cmeshtools::RemovePeriodicity(h1mesh);
      wgma::cmeshtools::RemovePeriodicity(mfmesh);
      src_an=nullptr;
    }
    
    
    if(match_an){
      match_data.cmesh = match_an->cmesh_hcurl;
      ComputeWpbcCoeffs(match_an, match_data.wgbc_k,
                      match_data.wgbc_f,true, {},
                      simdata.n_threads);
      auto h1mesh = match_an->cmesh_h1;
      auto mfmesh = match_an->cmesh_mf;
      wgma::cmeshtools::RemovePeriodicity(h1mesh);
      wgma::cmeshtools::RemovePeriodicity(mfmesh);
      match_an=nullptr;
    }
  }


  /*
    dirichlet boundary connects should not be restricted, otherwise
    this will result in all the equations on the same dependency
    being removed as well
  */
  std::set<int64_t> bound_connects_left, bound_connects_right;

  
  
  //eq num for obtaining reflection and transmittivity
  int64_t refl_pos{-1},trans_pos{-1};
    
  const REAL computed_res =
    RestrictDofsAndSolve(scatt_mesh_wpbc, src_data, match_data,
                         src_coeffs, nmodes_left,nmodes_right,
                         mats_near_wpbc,simdata,
                         last_sol,
                         refl_pos,
                         trans_pos
                         );
  //plot
  if(simdata.export_vtk_scatt){vtk.Do();}
  //get reflection and transmission

  TPZFMatrix<CSTATE> &sol = scatt_mesh_wpbc->Solution();
    
  {
    //for now we assume they are sequential
    const int nm = simdata.source_coeffs.size();
    std::cout<<"wavelength: "<<simdata.lambda;
    for(int i = 0; i < nm; i++){
      const auto s11 = sol.GetVal(refl_pos+i,0)-src_coeffs[i];
      const auto s21 =  trans_pos >= 0 ? sol.GetVal(trans_pos+i,0) : 0;
      const auto ref = std::abs(s11)*std::abs(s11);
      const auto trans = std::abs(s21)*std::abs(s21);
      std::cout<<" src "<<src_coeffs[i]
               <<" s11 "<<s11
               <<" ref "<<ref
               <<" s21 "<<s21
               <<" trans "<<trans
               <<" t + r "<<trans+ref<<std::endl;
    }
    std::string outputfile = simdata.prefix+"_reflection.csv";
    std::ofstream ost;
    ost.open(outputfile, std::ios_base::app);
    ost << std::setprecision(std::numeric_limits<STATE>::max_digits10);
    ost << simdata.lambda<<',';
    for(int i = 0; i < nm; i++){
      const CSTATE s11 = sol.GetVal(refl_pos+i,0)-src_coeffs[i];
      const CSTATE s21 = trans_pos >= 0 ? sol.GetVal(trans_pos+i,0) : 0;
      const char s11_sign = s11.imag() > 0 ? '+' : '-';
      const char s21_sign = s21.imag() > 0 ? '+' : '-';
      ost <<s11.real()<<s11_sign<<std::abs(s11.imag())<<'j'<<','
          <<s21.real()<<s21_sign<<std::abs(s21.imag())<<'j'<<',';
    }
    //useful for debuggin weird results
    ost <<computed_res<<std::endl;
  }
  //removing restrictions
  wgma::cmeshtools::RemovePeriodicity(scatt_mesh_wpbc);
  scatt_mesh_wpbc->ComputeNodElCon();
  scatt_mesh_wpbc->CleanUpUnconnectedNodes();
  
  wgma::cmeshtools::RemovePeriodicity(src_data.cmesh);
  if(match_data.cmesh){wgma::cmeshtools::RemovePeriodicity(match_data.cmesh);}
  return computed_res;
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
  wgma::cmeshtools::PhysicalData modal_data;
  std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
  for(const auto &matname : mats){
    //now we remove the suffix
    const auto suffix_length = suffix.size();
    const auto name = matname.substr(0,matname.length()-suffix_length);
    const CSTATE n = simdata.refractive_indices.at(name);
    modal_mats[matname] = std::pair<CSTATE,CSTATE>(n,1.);
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
                                          {0,0,0}, modal_data, modal_dim);
  return modal_data;
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
                                          {0,0,0}, scatt_data);


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
    constexpr int dim{3};
    constexpr int probedim{dim-1};
    for(const auto &[name,id] : gmshmats[probedim]){
      const bool found_pattern =
        std::regex_search(name, rx_left) ||
        std::regex_search(name, rx_right);
      if(found_pattern){
        scatt_data.probevec.push_back({id,probedim});
      }
    }
  }
  
  std::set<int> volids;
  for(auto [id,dummy1,dummy2] : scatt_data.matinfovec){
    volids.insert(id);
  }
    
  mats_near_wpbc =
    UpdatePhysicalDataSplittedMats(gmesh, scatt_data, split_mats,
                                   volids, dim);
  const bool verbose{true};
  //we must not condense since we will change p order
  const bool condense{false};
  auto cmesh =
    wgma::scattering::CMeshScattering3DPeriodic(gmesh, pOrder, scatt_data,
                                                el_map,
                                                src_ids,
                                                lambda,scale,verbose,condense);
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

void ComputeWpbcCoeffs(ModalData& an,
                       TPZFMatrix<CSTATE> &wgbc_k, TPZVec<CSTATE> &wgbc_f,
                       const bool positive_z, const TPZVec<CSTATE> &coeff,
                       const int nthreads){
  auto mesh = an.cmesh_mf;

  TPZFMatrix<CSTATE>& sol_orig = mesh->Solution();
  const int neq = sol_orig.Rows();
  const int nsol = sol_orig.Cols();
  
  wgma::post::WaveguidePortBC<wgma::post::MultiphysicsIntegrator,1> wgbc(mesh);
  wgbc.SetNThreads(nthreads);
  TPZManVector<CSTATE,1000> betavec = an.eigenvalues;
  if(coeff.size()){
    wgbc.SetSrcCoeff(coeff);
  }
  wgbc.SetBeta(betavec);
  wgbc.ComputeContribution();
  wgbc.GetContribution(wgbc_k,wgbc_f);
}

REAL RestrictDofsAndSolve(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                          WpbcData& src_data,
                          WpbcData& match_data,
                          const TPZVec<CSTATE> &source_coeffs,
                          const int nmodes_src,
                          const int nmodes_match,
                          const std::set<int> &mats_near_wpbc,
                          const SimData &simdata,
                          TPZFMatrix<CSTATE> &sol_vec,
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
  if(strmtrx){
    //300 is the number of shape functions in a k4 hexahedron
    const int maxsz = std::max(nmodes_match,nmodes_src) + 540;
    const int bufsz = maxsz*maxsz;
    strmtrx->BufferSizeForUserMatrix(bufsz);
  }
  TPZSimpleTimer timer("WPBC:Assemble+solve",true);
  {
    TPZSimpleTimer tassemble("Assemble",true);
    if(sol_vec.Rows() > 0){
      std::cout<<"running with custom init vec"<<std::endl;
      const auto eqfilt = scatt_an.StructMatrix()->EquationFilter();
      int64_t neq {0};
      if(eqfilt.IsActive()){
        neq = eqfilt.NActiveEquations();
      }else{
        neq = scatt_mesh->NEquations();
      }
      if(sol_vec.Rows() != neq){
        sol_vec.Redim(neq,1);
      }
      scatt_an.SetInitVecCustom(sol_vec);
    }
    std::cout<<"Assembling..."<<std::endl;
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
      dynamic_cast<TPZFYsmpMatrix<CSTATE>*>(scatt_an.StructMatrix()->Create());
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
    int from_current = sol_vec.Rows() > 0 ? 1 : 0;
    SetupPrecond(scatt_an, indices, simdata.solver_tol, from_current);
  }
  TPZSimpleTimer tsolve("Solve",true);
  scatt_an.Solve();
  auto sol = scatt_an.Solution();
  const auto eqfilt = scatt_an.StructMatrix()->EquationFilter();
  if(eqfilt.IsActive()){
    const auto neq = eqfilt.NActiveEquations();
    sol_vec.Resize(neq,1);
    eqfilt.Gather(sol, sol_vec);
  }else{
    sol_vec = sol;
  }

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

  return scatt_an.GetResidual();
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
                  const std::set<int64_t> &indep_cons,
                  const REAL tol,
                  int from_current) {
  TPZSimpleTimer solve("SetupPrecond", true);
      
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