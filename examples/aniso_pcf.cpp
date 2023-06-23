/**
aniso_pcf.cpp

This target performs the modal analysis of a step-index optical fiber.

It also illustrates the usage of the SLEPc solver and how to create
complex curved structured meshes in NeoPZ.
***/

//wgma includes
#include "wganalysis.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include "twisted_wgma.hpp"
#include <Electromagnetics/TPZCylindricalPML.h>
#include <TPZYSMPPardiso.h>
#include "util.hpp"
//pz includes
#include <pzcmesh.h>                     //for TPZCompMesh
#include <pzgmesh.h>                     //for TPZGeoMesh
#include <pzlog.h>                       //for TPZLogger
#include <TPZElectromagneticConstants.h> //for pzelectromag::cZero
#include <TPZQuadEigenSolver.h>        //for TPZQuadEigenSolver
#include <TPZVTKGenerator.h>
#include <TPZPardisoSolver.h>
#include <post/solutionnorm.hpp>
#include <TPZSimpleTimer.h>              //for TPZSimpleTimer

// Sets geometric info regarding all the circles in the mesh
TPZVec<wgma::gmeshtools::ArcData> SetUpArcData(std::string_view filename,
                                               const REAL scale);

int main(int argc, char *argv[]) {
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif
  /***********************
   * setting the problem *
   ***********************/
  
  //refractive index of the holes in fiber
  constexpr STATE n_air{1};
  //refractive index of the cladding
  constexpr STATE n_clad{1.444024};
  // operational wavelength
  constexpr STATE lambda{1.55};
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
  constexpr int pOrder2D{2};
  
  constexpr STATE modal_alphaPMLx{0.006};
  constexpr STATE modal_alphaPMLy{0.006};
  /******************
   * solver options *
   ******************/

  // number of threads to use
  const int nThreads = std::thread::hardware_concurrency();
  const int nThreadsDebug = std::thread::hardware_concurrency();
  // how to sort eigenvaluesn
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetMagnitude};
  constexpr int nEigenpairs{3};
  constexpr int krylovDim{30};
  constexpr CSTATE target = 1.43464883801848*1i;

  constexpr bool computeVectors{true};
  /*********************
   * exporting options *
   *********************/

  // whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  // path for output files
  const std::string path {"res_aniso_pcf/"};
  // common prefix for both meshes and output files
  const std::string basisName{"pcf"};
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
  constexpr bool optimizeBandwidth{true};
  /*
    The equations corresponding to homogeneous dirichlet boundary condition(PEC)
    can be filtered out of the global system in order to achieve better
    conditioning.
   */
  constexpr bool filterBoundaryEqs{true};
  /*
    Whether to use a non-linear representation for cylinders
   */
  constexpr bool arc3D{true};

  constexpr bool printGMesh{false};

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

  if(arc3D){
    /**
     in order to exactly represent all the circles in the mesh,
     the python script that generated the mesh also generated this .csv file
    **/
    const std::string arcfile{"meshes/"+basisName+"_circdata.csv"};
    auto arcdata = SetUpArcData(arcfile, scale);
    wgma::gmeshtools::SetExactArcRepresentation(gmesh, arcdata);
  }
  
  //print gmesh to .txt and .vtk format
  if(printGMesh)
  {
    //prefix for the gmesh files
    const auto filename = prefix+"_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh,filename);
  }

  //setting up cmesh data
  wgma::cmeshtools::PhysicalData modal_data;
  auto modal_cmesh = [gmesh,&gmshmats, &modal_data, scale](){
    // setting up cmesh data
    wgma::cmeshtools::PhysicalData modal_data;

    std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
    modal_mats["air"] = std::make_pair<CSTATE, CSTATE>(n_air*n_air, 1.);
    modal_mats["cladding"] = std::make_pair<CSTATE, CSTATE>(n_clad*n_clad, 1.);
    std::map<std::string, wgma::bc::type> modal_bcs;
    modal_bcs["modal_bnd"] = wgma::bc::type::PEC;

    //dimension of the modal analysis 
    constexpr int modal_dim{2};
    wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                            {modal_alphaPMLx,modal_alphaPMLy},
                                            modal_data, modal_dim);
    return wgma::wganalysis::CMeshWgmaAniso2D(gmesh,pOrder2D, modal_data,lambda, scale);
  }();

  //now we replace all the materials by their torsioned version

  {
    const STATE alpha{0.05};
    auto mfmesh = modal_cmesh[0];
    auto &matvec = mfmesh->MaterialVec();
    TPZVec<TPZBndCond*> bndmats;
    //first we get all bc mats
    for(auto &matpair : matvec){
      auto *mat = dynamic_cast<TPZBndCond*>(matpair.second);
      if(mat){
        bndmats.push_back(mat);
      }
    }

    std::set<TPZMaterial*> mats_to_delete;
    for(auto &matpair : matvec){
      const auto id = matpair.first;
      auto *mat = dynamic_cast<TPZAnisoWgma*>(matpair.second);
      auto *pml = dynamic_cast<TPZCombinedSpacesCylindricalPML<TPZAnisoWgma>*>(matpair.second);
      if(mat){
        auto twisted_mat = new wgma::materials::TwistedWgma(*mat);
        mats_to_delete.insert(mat);
        twisted_mat->SetAlpha(alpha);
        if(pml){
          REAL rMin{0},dR{0};
          STATE alphapml{0};
          auto twisted_pml = new wgma::materials::TwistedWgmaPML(id,*twisted_mat);
          pml->GetAttR(rMin,alphapml,dR);
          twisted_pml->SetAttR(rMin,alphapml,dR);
          delete twisted_mat;
          twisted_mat = twisted_pml;
        }
        matpair.second = twisted_mat;
        for(auto bnd : bndmats){
          if(bnd->Material()->Id() == id){
            bnd->SetMaterial(twisted_mat);
            break;
          }
        }
      }
    }
    for (auto mat : mats_to_delete){delete mat;}
  }

  
  //WGAnalysis class is responsible for managing the modal analysis
  wgma::wganalysis::WgmaAniso2D analysis(modal_cmesh,nThreads,optimizeBandwidth,filterBoundaryEqs);

  TPZAutoPointer<TPZQuadEigenSolver<CSTATE>> krylov_solver =
    new TPZQuadEigenSolver<CSTATE>();
  krylov_solver->SetKrylovDim(krylovDim);

  krylov_solver->SetTarget(target);
  krylov_solver->SetNEigenpairs(nEigenpairs);
  krylov_solver->SetEigenSorting(sortingRule);

  analysis.SetSolver(*krylov_solver);

  analysis.Assemble();

  {
    auto &solv = analysis.GetSolver();

    auto sparse_k =
      TPZAutoPointerDynamicCast<TPZFYsmpMatrixPardiso<CSTATE>>(solv.MatrixK());
    auto &pardiso_control = sparse_k->GetPardisoControl();
    auto sys = SymProp::NonSym;
    auto prop = TPZPardisoSolver<CSTATE>::MProperty::EIndefinite;
    pardiso_control.SetMatrixType(sys,prop);
    auto param = pardiso_control.GetParam();
      //fParam[0] No default values
    param[0] = 1;
    //param[1]  use Metis for the ordering
    param[1] = 2;
    /*param[3]  Preconditioned CGS/CG. 
      0 = // No iterative-direct algorithm
      10*L+K
      L = stoppping criterion: 10^-L
      K = 
      0: The factorization is always computed as required by phase
      1: CGS iteration replaces the computation of LU. 
      The preconditioner is LU that was computed at a previous step
      (the first step or last step with a failure) in a sequence of
      solutions needed for identical sparsity patterns.
      2: CGS iteration for symmetric positive definite matrices
      Replaces the computation of LLt. The preconditioner is LLT
      that was computed at a previous step
      (the first step or last step with a failure)
      in a sequence of solutions needed for identical sparsity patterns. 
    */
    //param[4]  No user fill-in reducing permutation
    param[3] = 0;
    param[4] = 0;

    //param[7]  Maximum number of iterative refinement steps that the solver performs when perturbed pivots are obtained during the numerical factorization. 
    param[7] = 8;
      	
    //param[8]  Tolerance level for the relative residual in the iterative refinement process. (10^-{param[8]})
    param[8] = 23;
    //param[9]  Perturb the pivot elements with 1E-param[9]
    param[9] = 23;
    //param[10]  Use nonsymmetric permutation and scaling MPS
    param[10] = 0;

      
    //param[12]  Maximum weighted matching algorithm is switched-off (default for symmetric).
    param[12] = 0;
    //param[26] Whether to check matrix data
    param[26] = 1;
    //param[59]  Do not use OOC
    param[59] = 0;
    pardiso_control.SetParam(param);
    
  }
  analysis.Solve(computeVectors);

  if (!computeVectors && !exportVtk) return 0;

  const std::string plotfile = prefix+"_modal";

  
  {
    
    TPZSimpleTimer tpostprocess("Post processing");
    TPZVec<std::string> fvars = {
      "Ez_real",
      "Ez_abs",
      "Et_real",
      "Et_abs",
      "Material"};
    auto vtk = TPZVTKGenerator(modal_cmesh[0], fvars, plotfile, vtkRes);
    auto ev = analysis.GetEigenvalues();
    for (int isol = 0; isol < ev.size(); isol++) {
      auto currentKz = CSTATE(-1i)*ev[isol];
      std::cout<<"\rPost processing step "<<isol+1<<" out of "<<ev.size()
               <<"(kz = "<<currentKz<<")"<<std::endl;
      analysis.LoadSolution(isol);
      vtk.Do();
    }
  }
  return 0;
}

// Sets geometric info regarding all the circles in the mesh
TPZVec<wgma::gmeshtools::ArcData> SetUpArcData(std::string_view filename,
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
  // we expect xc, yc, zc, r (in um), and matid
  TPZVec<wgma::gmeshtools::ArcData> arcs;
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
  return arcs;
}