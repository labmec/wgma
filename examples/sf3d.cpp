/**
stepfiber_3d.cpp

This target performs the modal analysis of a step-index optical
and then the subsequent scattering analysis.
***/

// wgma includes
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include <wganalysis.hpp>
#include <scattering.hpp>
#include <util.hpp>
#include <precond.hpp>
#include <slepcepshandler.hpp>
#include <post/solutionnorm.hpp>
// pz includes
#include <MMeshType.h>      //for MMeshType
#include <TPZSimpleTimer.h> //for TPZSimpleTimer
#include <pzlog.h>          //for TPZLogger
#include <TPZVTKGenerator.h>
#include <TPZPardisoSolver.h>
#include <TPZYSMPMatrix.h>
#include <TPZSYSMPMatrix.h>
#include <regex>//for string search
#include <pzcheckgeom.h>
#include <tpzchangeel.h>
#include <pzstepsolver.h>
#include <pzseqsolver.h>
#include <pzmgsolver.h>

#include <thread>

// Sets geometric info regarding all the cylinders in the mesh
TPZVec<wgma::gmeshtools::CylinderData> SetUpCylData(std::string_view filename,
                                                    const REAL scale);


void ComputeModes(wgma::wganalysis::Wgma2D &an,
                  TPZAutoPointer<TPZCompMesh> cmesh,
                  const bool orthogonalise,
                  const int nThreads);

void PostProcessModes(wgma::wganalysis::Wgma2D &an,
                      TPZAutoPointer<TPZCompMesh> cmesh,
                      std::string filename,
                      const int vtkres,
                      std::set<int> sols = {});

int main(int argc, char *argv[]) {
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif
  /***********************
   * setting the problem *
   ***********************/
  
  //refractive index of the fibers core
  constexpr STATE n_core{1.4457};
  //refractive index of the fibers cladding
  constexpr STATE n_clad{1.4378};
  //magnetic permeability of the core
  constexpr STATE ur_core{1};
  //magnetic permeability of the cladding
  constexpr STATE ur_clad{1};
  // operational wavelength
  constexpr STATE lambda{4.0};
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
  // polynomial order to be used in the modal analysis
  constexpr int pOrder2D{1};
  // polynomial order to be used in the scattering analysis
  constexpr int pOrder3D{1};
  //PML attenuation constant
  /*//working for lambda = 30
  constexpr STATE modal_alphaPMLx{0.5};
  constexpr STATE modal_alphaPMLy{0.5};
  
  constexpr STATE scatt_alphaPMLx{1.0};
  constexpr STATE scatt_alphaPMLy{1.0};
  constexpr STATE scatt_alphaPMLz{5.0};*/

  /*
    0.009: good ez, bad et
    0.015: better, et still weird
    0.05: ez already spreading to PML
   */
  constexpr STATE modal_alphaPMLx{0.4};
  constexpr STATE modal_alphaPMLy{0.4};
  
  constexpr STATE scatt_alphaPMLx{0.4};
  constexpr STATE scatt_alphaPMLy{0.4};
  constexpr STATE scatt_alphaPMLz{0.4};
  /******************
   * solver options *
   ******************/

  // number of threads to use
  const int nThreads = std::thread::hardware_concurrency();
  const int nThreadsDebug = std::thread::hardware_concurrency();
  // how to sort eigenvaluesn
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  bool usingSLEPC {true};
  constexpr int nEigenvalues{1};
  constexpr CSTATE target = -2.09;//n_core*n_core;
  constexpr bool symMat{false};
  /*********************
   * exporting options *
   *********************/

  // whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  // path for output files
  const std::string path {"res_sf3d/"};
  // common prefix for both meshes and output files
  const std::string basisName{"sf3d"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(prefix));

  
  // resolution of the .vtk file in which the solution will be exported
  constexpr int vtkRes{0};

  /********************
   * advanced options *
   ********************/

  // reorder the equations in order to optimize bandwidth
  constexpr bool optimizeBandwidth2D{true};
  // reorder the equations in order to optimize bandwidth
  constexpr bool optimizeBandwidth3D{true};
  /*
    The equations corresponding to homogeneous dirichlet boundary condition(PEC)
    can be filtered out of the global system in order to achieve better
    conditioning.
   */
  constexpr bool filterBoundaryEqs2D{true};
  constexpr bool filterBoundaryEqs3D{true};
  /*
    Whether to use a non-linear representation for cylinders
   */
  constexpr bool cyl3D{true};

  //! Whether to use a direct solver
  constexpr bool direct{false};

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
  constexpr bool verbosity_lvl{false};
  auto gmesh = wgma::gmeshtools::ReadGmshMesh(meshfile, scale,
                                              gmshmats,verbosity_lvl);

  if(cyl3D){
    /**
     in order to exactly represent all the circles in the mesh,
     the python script that generated the mesh also generated this .csv file
    **/
    const std::string cylfile{"meshes/"+basisName+"_cyldata.csv"};
    auto cyldata = SetUpCylData(cylfile, scale);
    wgma::gmeshtools::SetExactCylinderRepresentation(gmesh, cyldata);
  }
  

  /********************************
   * cmesh(modal analysis)        *
   ********************************/

  wgma::cmeshtools::PhysicalData modal_data;
  auto modal_l_cmesh = [gmesh,&gmshmats,&modal_data, scale](){
    // setting up cmesh data

    std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
    modal_mats["src_clad"] = std::make_pair<CSTATE, CSTATE>(n_clad*n_clad, 1.);
    modal_mats["src_core"] = std::make_pair<CSTATE, CSTATE>(n_core*n_core, 1.);
    std::map<std::string, wgma::bc::type> modal_bcs;
    modal_bcs["modal_bnd"] = wgma::bc::type::PEC;
    //dimension of the modal analysis 
    constexpr int modal_dim{2};
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                            {modal_alphaPMLx,modal_alphaPMLy},
                                            modal_data, modal_dim);
    //we must now filter the 2D PMLs
    std::vector<TPZAutoPointer<wgma::pml::data>>  pmlvec;
    for(const auto &pml : modal_data.pmlvec){
      const std::string pattern{"src_clad"};
      const auto rx = std::regex{pattern, std::regex_constants::icase };
    
      const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
      if(found_pattern){pmlvec.push_back(pml);}
    }
    modal_data.pmlvec = pmlvec;
    return wgma::wganalysis::CMeshWgma2D(gmesh,pOrder2D, modal_data,lambda, scale);
  }();

  /******************************
   * solve(modal analysis left) *
   ******************************/
  TPZAutoPointer<wgma::wganalysis::Wgma2D>
    modal_l_an = new wgma::wganalysis::Wgma2D(modal_l_cmesh, nThreads,
                                              optimizeBandwidth2D, filterBoundaryEqs2D);
  {
    auto solver_left =
      wgma::wganalysis::SetupSolver(target, nEigenvalues, sortingRule, usingSLEPC);
    modal_l_an->SetSolver(*solver_left);
  }
  std::string modal_left_file{prefix+"_modal"};
  if(usingSLEPC){
    modal_left_file += "_slepc";
  }else{
    modal_left_file += "_krylov";
  }
  constexpr bool ortho{true};
  ComputeModes(modal_l_an, modal_l_cmesh[0], ortho, nThreads);
  if(exportVtk){
    PostProcessModes(modal_l_an, modal_l_cmesh[0], modal_left_file, vtkRes);
  }
  
  //materials that will represent our source
  std::set<int> src_ids;
  {
    const std::string srcMat[] = {"src_core", "src_clad"};
    for(const auto &mat : srcMat){
      if(gmshmats[2].count(mat)){
        src_ids.insert(gmshmats[2].at(mat));
      }
    }
    for(const auto &pml : modal_data.pmlvec){
      const std::string pattern{"src_clad"};
      const auto rx = std::regex{pattern, std::regex_constants::icase };
    
      const bool found_pattern = std::regex_search(*(pml->names.begin()), rx);
      if(found_pattern){
        for (auto p : pml->ids){
          src_ids.insert(p);
        }
      }
    }
  }
  
  
  auto scatt_cmesh = [gmesh,&gmshmats, src_ids](){
    
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData scatt_data;
    std::map<std::string, std::pair<CSTATE, CSTATE>> scatt_mats;
    scatt_mats["core"] = std::make_pair<CSTATE, CSTATE>(n_core*n_core, 1.);
    scatt_mats["cladding"] = std::make_pair<CSTATE, CSTATE>(n_clad*n_clad,1.);
    std::map<std::string, wgma::bc::type> scatt_bcs;
    scatt_bcs["scatt_bnd"] = wgma::bc::type::PEC;

    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, scatt_mats, scatt_bcs,
                                            {scatt_alphaPMLx,
                                             scatt_alphaPMLy,
                                             scatt_alphaPMLz},
                                            scatt_data);
    
    //materials in which we would like to evaluate the solution
    const std::string probeMats[] = {};

    TPZSimpleTimer scatt("Scattering mesh", true);
    return wgma::scattering::CMeshScattering3D(gmesh, pOrder3D, scatt_data,src_ids,
                                               lambda,scale);
  }();

  
  constexpr int which_src {0};
  modal_l_an->LoadSolution(which_src);
  const auto beta = sqrt(-modal_l_an->GetEigenvalues()[which_src]);

  {
    TPZSimpleTimer tscatt("Load source",true);
    wgma::scattering::SourceWgma src;
    src.id = src_ids;
    src.modal_cmesh = modal_l_cmesh[0];
    wgma::scattering::LoadSource2D(scatt_cmesh, src);
    wgma::scattering::SetPropagationConstant(scatt_cmesh, beta);
  }

  modal_l_an = nullptr;
  modal_l_cmesh.resize(0);

  /*********************
   * solve(scattering) *  
   *********************/  
  TPZSimpleTimer tscatt("Scattering", true);

  auto scatt_an = wgma::scattering::Analysis(scatt_cmesh, nThreadsDebug,
                                             optimizeBandwidth3D, filterBoundaryEqs3D,
                                             symMat);
  
  {
    TPZSimpleTimer tscatt("Assemble",true);
    scatt_an.Assemble();

    auto mat = scatt_an.GetSolver().Matrix();
    mat->SetSymmetry(SymProp::Sym);
    mat->SetDefPositive(false);
  }

  if(direct){
    //get pardiso control
    auto *pardiso = scatt_an.GetSolver().GetPardisoControl();
    pardiso->SetMessageLevel(1);
  
    pardiso->ResetParam();
    constexpr auto sys_type = SymProp::Sym;
    constexpr auto prop = TPZPardisoSolver<CSTATE>::MProperty::EIndefinite;
    pardiso->SetMatrixType(sys_type,prop);
    TPZSimpleTimer tscatt("Solve",true);
    scatt_an.Solve();
    pardiso->FreePardisoMemory();
  }
  else{
    TPZSimpleTimer solve("precond+solve", true);
    constexpr REAL tol = 5e-8;
      
    auto &solver = dynamic_cast<TPZStepSolver<CSTATE>&>(scatt_an.GetSolver());

    TPZAutoPointer<TPZMatrixSolver<CSTATE>> precond;
    {
      TPZSimpleTimer pre("precond",true);
      const auto eqfilt = scatt_an.StructMatrix()->EquationFilter();
      std::set<int> bnd_ids;
      bnd_ids.insert(gmshmats[2].at("scatt_bnd"));
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
    //solves the system
    scatt_an.Solve();
  }

  const std::string scatt_file = prefix+"_scatt";
    
  {
    TPZSimpleTimer tpostprocess("Post processing");
    TPZVec<std::string> fvars = {
      "Field_real",
      "Field_imag",
      "Field_abs"};
    auto vtk = TPZVTKGenerator(scatt_cmesh, fvars, scatt_file, vtkRes);
    vtk.Do();
  }
  
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
                  TPZAutoPointer<TPZCompMesh> cmesh,
                  const bool orthogonalise,
                  const int nThreads)
{
  
  TPZSimpleTimer analysis("Modal analysis");
  an.Assemble();

  auto matA = an.GetSolver().MatrixA();
  auto matB = an.GetSolver().MatrixB();

  matA->SetDefPositive(false);
  matB->SetDefPositive(false);
  static constexpr bool computeVectors{true};
  an.Solve(computeVectors);
  //let us discard all solutions with negative real part of beta

  auto &ev = an.GetEigenvalues();
  auto &evec = an.GetEigenvectors();

  //load all obtained modes into the mesh
  an.LoadAllSolutions();
  using namespace wgma::post;
  if(orthogonalise){
    //leave empty for all valid matids
    std::set<int> matids {};
    constexpr bool conj = true;
    auto ortho = SolutionNorm<MultiphysicsIntegrator>(cmesh, matids, conj, nThreads);
    //orthogonalise the modes
    ortho.Normalise();
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