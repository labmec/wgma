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
#include "materials/acousticmodesomega.hpp"
#include "slepcepshandler.hpp"
#include <TPZYSMPPardiso.h>
#include "util.hpp"
//pz includes
#include <pzcmesh.h>                     //for TPZCompMesh
#include <pzgmesh.h>                     //for TPZGeoMesh
#include <pzlog.h>                       //for TPZLogger
#include <TPZElectromagneticConstants.h> //for pzelectromag::cZero
#include <TPZKrylovEigenSolver.h>        //for TPZKrylovEigenSolver
#include <TPZVTKGenerator.h>
#include <TPZPardisoSolver.h>
#include <post/solutionnorm.hpp>
#include <TPZSimpleTimer.h>              //for TPZSimpleTimer

// Sets geometric info regarding all the circles in the mesh
TPZVec<wgma::gmeshtools::ArcData> SetUpArcData(std::string_view filename,
                                               const REAL scale);


TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const int neigenpairs, const CSTATE target,
            const int krylovDim,
            TPZEigenSort sorting, bool usingSLEPC);

TPZAutoPointer<TPZCompMesh>
CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
            const int pOrder, const CSTATE beta, const STATE scale,
            const STATE torsion,
            const TPZVec<std::map<std::string, int>> &gmshmats,
            const std::map<std::string, std::tuple<STATE, STATE, STATE>> &modal_mats,
            const std::map<int, wgma::bc::type> &modal_bcs);

TPZEigenAnalysis
CreateAnalysis(TPZAutoPointer<TPZCompMesh> modal_cmesh, const int nThreads,
               const bool optimizeBandwidth);

constexpr STATE ComputeMu(const STATE E, const STATE v){
    return E/(2*(1+v));
}
constexpr STATE ComputeLambda(const STATE E, const STATE v){
    return E*v/((1+v)*(1-2*v));
}

using namespace std::complex_literals;
int main(int argc, char *argv[]) {
#ifdef PZ_LOG
    /**if the NeoPZ library was configured with log4cxx,
     * the log should be initialised as:*/
    TPZLogger::InitializePZLOG();
#endif
    /***********************
     * setting the problem *
     ***********************/

    //the mesh is described in micrometers
    constexpr STATE char_length = {1/1e6};
    
    //taken from https://pubs.aip.org/aip/app/article/4/7/071101/1024179/Brillouin-optomechanics-in-nanophotonic-structures
    //value in tonne/m^3
    constexpr STATE rho_clad{2.203};
    //this value is in GPa
    constexpr STATE young_clad{73.1};
    constexpr STATE poisson_clad{0.17};
    constexpr STATE lambda_clad = ComputeLambda(young_clad, poisson_clad);
    constexpr STATE mu_clad = ComputeMu(young_clad, poisson_clad);
    constexpr STATE torsion{0.0};
    /*
      we set \beta_{mech} = 2\beta_{electromagnetic}
      using the value from
      Leaky Modes In Twisted Microstructured Optical Fibers
      A. Nicolet, F. Zolla, Y. O. Agha & S. Guenneau,
      2007,
      neff = 1.43877448428430
      thus
      beta = 5832314 rad/m
      
     */
    constexpr CSTATE beta{2*5.832314};

    constexpr REAL scale_geom{1};
    constexpr REAL scale_mat{1};

    /******************
     *  fem options   *
     ******************/
    // polynomial order to be used in the modal analysis
    constexpr int pOrder{2};
  
    /******************
     * solver options *
     ******************/

    // number of threads to use
    const int nThreads = std::thread::hardware_concurrency();
    const int nThreadsDebug = 0;
    // how to sort eigenvaluesn
    constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
    constexpr int nEigenvalues{10};
    constexpr bool usingSLEPC{true};
    constexpr int krylovDim{150};
    //we search for dominant modes
    const CSTATE target = 0;

    constexpr bool computeVectors{true};
    /*********************
     * exporting options *
     *********************/

    // whether to export the solution as a .vtk file
    constexpr bool exportVtk{true};
    // path for output files
    const std::string path {"res_aniso_acoustics_omega/"};
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
    auto gmesh = wgma::gmeshtools::ReadGmshMesh(meshfile, scale_geom,
                                                gmshmats,verbosity_lvl);

    const std::string arcfile{"meshes/"+basisName+"_circdata.csv"};
    auto arcdata = SetUpArcData(arcfile, scale_geom);
    if(arc3D){
        /**
           in order to exactly represent all the circles in the mesh,
           the python script that generated the mesh also generated this .csv file
        **/
        wgma::gmeshtools::SetExactArcRepresentation(gmesh, arcdata);
    }
  
    //stores the mat ids of all circles
    std::set<int> circles_ids;
    {
        for(auto &arc : arcdata){
            circles_ids.insert(arc.m_matid);
        }
        //now we find the outer boundary (outside PML) and remove it
        constexpr int bcdim{1};
        bool found{false};
        for(auto [name,id] : gmshmats[bcdim]){
            if(name == "modal_bnd"){
              circles_ids.erase(id);
              found = true;
              break;
            }
        }
        if(!found){
            DebugStop();
        }
    }

    //let us refine near the circles
    constexpr int nrefdir{2};
    wgma::gmeshtools::DirectionalRefinement(gmesh, circles_ids, nrefdir);
    
    //print gmesh to .txt and .vtk format
    if(printGMesh)
    {
        //prefix for the gmesh files
        const auto filename = prefix+"_gmesh";
        wgma::gmeshtools::PrintGeoMesh(gmesh,filename);
    }

    //setting up cmesh data
    //rho, lambda, mu
    std::map<std::string, std::tuple<STATE, STATE, STATE>> modal_mats;
    modal_mats["cladding"] = std::make_tuple(rho_clad,lambda_clad,mu_clad);
  
    std::map<int, wgma::bc::type> modal_bcs;
    //all air-si interfaces are modelled as free surfaces
    for(auto id : circles_ids){
        modal_bcs[id] = wgma::bc::type::PMC;
    }

    auto modal_cmesh = CreateCMesh(gmesh,pOrder,beta,scale_mat,torsion,gmshmats,modal_mats,modal_bcs);

  
    //WGAnalysis class is responsible for managing the modal analysis
    auto analysis = CreateAnalysis(modal_cmesh,nThreads,optimizeBandwidth);

    auto solver = SetupSolver(nEigenvalues,target, krylovDim, sortingRule, usingSLEPC);
  
    analysis.SetSolver(solver);
    
    {
        TPZSimpleTimer timer("Assemble",true);
        analysis.Assemble();
    }
    
    
    analysis.SetComputeEigenvectors(computeVectors);
    analysis.Solve();

    if (!computeVectors && !exportVtk) return 0;

    const std::string plotfile = prefix+"_modal";

  
    {
    
        TPZSimpleTimer tpostprocess("Post processing");
        TPZVec<std::string> fvars = {
            "u_real",
            "u_abs"};
        auto vtk = TPZVTKGenerator(modal_cmesh, fvars, plotfile, vtkRes);
        auto ev = analysis.GetEigenvalues();
        const auto &eigenvectors = analysis.GetEigenvectors();
        const auto neq = eigenvectors.Rows();
        for (int isol = 0; isol < ev.size(); isol++) {
            auto currentOmega = sqrt(ev[isol])/scale_mat;
            std::cout<<"\rPost processing step "<<isol+1<<" out of "<<ev.size()
                     <<ev[isol]<<" (w = "<<currentOmega<<")"<<std::endl;
            
            TPZFMatrix<CSTATE> sol(neq, 1, (CSTATE*)eigenvectors.Elem() + neq*isol, neq);
            modal_cmesh->LoadSolution(sol);
            vtk.Do();
        }
    }
    return 0;
}

TPZAutoPointer<TPZCompMesh>
CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
            const int pOrder, 
            const CSTATE beta, const STATE scale,
            const STATE torsion,
            const TPZVec<std::map<std::string, int>> &gmshmats,
            const std::map<std::string, std::tuple<STATE, STATE, STATE>> &modal_mats,
            const std::map<int, wgma::bc::type> &modal_bcs)
{
    ///Defines the computational mesh based on the geometric mesh
    constexpr bool is_complex{true};
    TPZAutoPointer<TPZCompMesh>  cmesh = new TPZCompMesh(gmesh,is_complex);
    //using traditional H1 elements
    cmesh->SetAllCreateFunctionsContinuous();

    std::set<int> vol_ids;
    constexpr int modaldim{2};
    for(auto [name,id] : gmshmats[modaldim]){
        if(modal_mats.count(name)){
            vol_ids.insert(id);
            auto [rho, lambda, mu] = modal_mats.at(name);
            auto *mat = new wgma::materials::AcousticModesOmega(id,mu,lambda,rho,beta,scale);
            mat->SetTorsion(torsion);
            cmesh->InsertMaterialObject(mat);
        }
    }

    TPZFNMatrix<1, CSTATE> val1(1, 1, 0);
    TPZManVector<CSTATE,1> val2(1, 0.);

    std::set<int> bc_ids;
    for(const auto [id, type] : modal_bcs){
        auto res = wgma::gmeshtools::FindBCNeighbourMat(gmesh, id, vol_ids);

            
        if(!res.has_value()){
            std::cout<<__PRETTY_FUNCTION__
                     <<"\nwarning: could not find neighbour of bc "<<id<<std::endl;
            continue;
        }
        bc_ids.insert(id);
        const int bctype = wgma::bc::to_int(type);
        const int volid = res.value();
        auto *matWG =
            dynamic_cast<TPZMaterialT<CSTATE>*>(cmesh->FindMaterial(volid));
        auto *bcmat = matWG->CreateBC(matWG,id,bctype,val1,val2);
        cmesh->InsertMaterialObject(bcmat);
    }

    std::set<int> all_ids;

    std::set_union(vol_ids.begin(), vol_ids.end(),
                   bc_ids.begin(), bc_ids.end(),
                   std::inserter(all_ids, all_ids.begin()));

    std::cout<<"ids: "<<std::endl;
    for(auto id : all_ids){
      std::cout<<id<<std::endl;
    }
    //seting the polynomial order in the computational mesh
    cmesh->SetDefaultOrder(pOrder);
    //actually creates the computational elements
    cmesh->AutoBuild(all_ids);
    return cmesh;
}


#include <TPZCutHillMcKee.h>
#include <TPZSpStructMatrix.h>

TPZEigenAnalysis
CreateAnalysis(TPZAutoPointer<TPZCompMesh> cmesh, const int nThreads,
               const bool optimizeBandwidth)
{
    TPZEigenAnalysis analysis(cmesh, RenumType::ECutHillMcKee);
    std::cout<<"neq: "<<cmesh->NEquations()<<std::endl;
    TPZAutoPointer<TPZStructMatrix> strmtrx =
      new TPZSpStructMatrix<CSTATE>(cmesh);

    strmtrx->SetNumThreads(nThreads);
    analysis.SetStructuralMatrix(strmtrx);
    return analysis;
}

TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const int neigenpairs, const CSTATE target,
            const int krylovDim,
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
    krylov_solver->SetNEigenpairs(neigenpairs);
    krylov_solver->SetAsGeneralised(true);
    krylov_solver->SetKrylovDim(krylovDim);
    
    solver = krylov_solver;
  }
  
  solver->SetEigenSorting(sorting);
  return solver;
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
