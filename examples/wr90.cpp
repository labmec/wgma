/**
wr90.cpp

This target is the simples example on how to use NeoPZ for the modal analysis
of a metal-backed electromagnetic waveguide using the WR90 waveguide as a model.
***/

//wgma includes
#include "gmeshtools.hpp"
#include "cmeshtools.hpp"
#include "wganalysis.hpp"
#include "util.hpp"
//pz includes
#include <MMeshType.h>                   //for MMeshType
#include <pzcmesh.h>                     //for TPZCompMesh
#include <pzgmesh.h>                     //for TPZGeoMesh
#include <TPZGeoMeshTools.h>             //for TPZGeoMeshTools::CreateGeoMeshOnGrid
#include <pzlog.h>                       //for TPZLogger
#include <TPZElectromagneticConstants.h> //for pzelectromag::cZero
#include <TPZKrylovEigenSolver.h>        //for TPZKrylovEigenSolver
#include <TPZVTKGenerator.h>             //for TPZVTKGenerator
#include <TPZSimpleTimer.h>              //for TPZSimpleTimer
/**
   @brief Creates the geometrical mesh associated with a rectangular metallic waveguide. 
The domain can be divided in half (width-wise) for exploring symmetry. */
TPZAutoPointer<TPZGeoMesh>
CreateGMeshRectangularWaveguide(TPZVec<int> &matIdVec, const REAL &scale,
                                const REAL wDomain, const REAL hDomain,
                                TPZVec<int> &nDivs,
                                const MMeshType meshType,
                                const bool usingSymmetry);

int main(int argc, char *argv[]) {
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif

  /***********************
   * setting the problem *
   ***********************/
  
  // width of the waveguide
  constexpr REAL wDomain{0.02286};
  // height of the waveguide
  constexpr REAL hDomain{0.01016};
  //whether to split the domain in two for exploring symmetries
  constexpr bool usingSymmetry{false};
  //which BC to apply if symmetry is used
  constexpr wgma::bc::type sym{wgma::bc::type::PEC};
  // n divisions in x direction
  constexpr int nDivX{7};
  // n divisions in y direction
  constexpr int nDivY{7};
  // type of elements
  constexpr MMeshType meshType{MMeshType::ETriangular};
  // TPZManVector<Type,N> is a vector container with static + dynamic storage.
  // one can also use TPZVec<Type> for dynamic storage
  TPZManVector<int, 2> nDivs = {nDivX, nDivY};
  // polynomial order to be used in the approximation
  constexpr int pOrder{3};
  // magnetic permeability
  constexpr CSTATE ur{1};
  // electric permittivity
  constexpr CSTATE er{1};
  // operational frequency
  constexpr STATE freq = 25e9;
  // wavelength
  constexpr STATE lambda{pzeletromag::cZero/freq};
  /*
    Given the small dimensions of the domain, scaling it can help in 
    achieving good precision. Using 1./k0 as a scale factor results in 
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  constexpr REAL scale{lambda/(2*M_PI)};


  /******************
   * solver options *
   ******************/
  
  //number of threads to use
  constexpr int nThreads{8};
  //number of genvalues to be computed
  constexpr int nEigenpairs{10};
  //whether to compute eigenvectors (instead of just eigenvalues)
  constexpr bool computeVectors{true};
  //how to sort the computed eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  /*
   The simulation uses a Krylov-based Arnoldi solver for solving the
   generalised EVP. A shift-and-inverse spectral transform is applied in 
   the system for better convergence of the eigenvalues.
   The target variable should be close to the desired eigenvalue (in this case,
   the effective index neff). We are looking for neff=1, i.e., the modes with
  lowest cutoff frequency*/
  constexpr CSTATE target = -1.01;
  // Dimension of the krylov space to be used. Suggested to be at least nev * 10
  constexpr int krylovDim{-1};


  /*********************
   * exporting options *
   *********************/

  //whether to print the geometric mesh in .txt and .vtk formats
  constexpr bool printGMesh{false};
  //whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  //resolution of the .vtk file in which the solution will be exported
  constexpr int vtkRes{1};
  //if true, the real part of the electric fields is exported. otherwise, the magnitude
  constexpr bool printRealPart{true};
  //path for exporting file
  const std::string path{"res_wr90/"};
  wgma::util::CreatePath(path);

  /********************
   * advanced options *
   ********************/
  
  //reorder the equations in order to optimize bandwidth
  constexpr bool optimizeBandwidth{true};
  /*
    The equations corresponding to homogeneous dirichlet boundary condition(PEC)
    can be filtered out of the global system in order to achieve better conditioning.
   */
  constexpr bool filterBoundaryEqs{true};

  /*********
   * begin *
   *********/

  
  //scoped-timer 
  TPZSimpleTimer total("Total");

  /* Vector for storing materials(regions) identifiers for the 
     TPZGeoMeshTools::CreateGeoMeshOnGrid function.
     In this example, we have one material for the interior of the waveguide
     and one for each BC*/
  TPZManVector<int, 5> matIdVec({1,-1,-2,-3,-4});
  //creates geometric mesh
  auto gmesh = CreateGMeshRectangularWaveguide(matIdVec,scale,wDomain,hDomain,nDivs,
                                               meshType,usingSymmetry);

  
  //print gmesh to .txt and .vtk format
  if(printGMesh)
  {
    const auto filename = path+"wr90_gmesh_"+
      std::to_string(nDivX) + " _" + std::to_string(nDivY);
    wgma::gmeshtools::PrintGeoMesh(gmesh,filename);
  }
  //Let us fill the struct needed to CMeshWgma2D
  wgma::cmeshtools::PhysicalData data;

  //setting material info
  data.matinfovec.push_back(std::make_tuple(1,er,ur));
  //setting BC information
  data.bcvec.resize(4);
  for(int i = 0; i < 4; i++){
    data.bcvec[i].id = matIdVec[i+1];
    data.bcvec[i].t = i == 1 && usingSymmetry ?
      wgma::bc::type::PMC : wgma::bc::type::PEC;
    
  }
  /*
   The problem uses an H1 approximation space for the longitudinal component 
   and a HCurl approximation space for the transversal one. Therefore, three
  // computational meshes are generated. One for each space and a multiphysics mesh*/
  auto meshVec = wgma::wganalysis::CMeshWgma2D(gmesh,pOrder,data, lambda,scale);

  //Analysis class is responsible for managing the modal analysis
  wgma::wganalysis::Wgma2D analysis(meshVec,nThreads,optimizeBandwidth,filterBoundaryEqs);

  TPZAutoPointer<TPZEigenSolver<CSTATE>>solver{nullptr};

  TPZAutoPointer<TPZKrylovEigenSolver<CSTATE>> krylov_solver =
    new TPZKrylovEigenSolver<CSTATE>;
  TPZSTShiftAndInvert<CSTATE> st;
  krylov_solver->SetSpectralTransform(st);
  krylov_solver->SetKrylovDim(krylovDim);

  krylov_solver->SetTarget(target);
  krylov_solver->SetNEigenpairs(nEigenpairs);
  krylov_solver->SetAsGeneralised(true);
  krylov_solver->SetEigenSorting(sortingRule);

  analysis.SetSolver(*krylov_solver);


  analysis.Run(computeVectors);
  
  if (!computeVectors && !exportVtk) return 0;

  const std::string plotfile = path+"wr90_field";

  TPZSimpleTimer postProc("Post processing");

  TPZVec<std::string> fields = {
    "Ez_real",
    "Ez_abs",
    "Et_real",
    "Et_abs"};
  auto vtk = TPZVTKGenerator(meshVec[0], fields, plotfile, vtkRes);
  auto ev = analysis.GetEigenvalues();
  for (int isol = 0; isol < ev.size(); isol++) {
    auto currentKz = std::sqrt(-1.0*ev[isol]);
    std::cout<<"\rPost processing step "<<isol+1<<" out of "<<ev.size()
             <<"(kz = "<<currentKz<<")"<<std::endl;
    analysis.LoadSolution(isol);
    vtk.Do();
  }
  
  return 0;
}






TPZAutoPointer<TPZGeoMesh> CreateGMeshRectangularWaveguide(
    TPZVec<int> &matIdVec, const REAL &scale, const REAL wDomain,
    const REAL hDomain, TPZVec<int> &nDivs, const MMeshType meshType,
    const bool usingSymmetry)
{
  TPZSimpleTimer timer("Create gmesh");
  // dimension of the problem
  constexpr int dim{2};
  // lower left corner of the domain
  TPZManVector<REAL, 3> minX = {0, 0, 0};
  
  // upper right corner of the domain
  TPZManVector<REAL, 3> maxX = {
    usingSymmetry ? wDomain /( scale * 2) : wDomain / scale,
    hDomain / scale,
    0};
  
  // whether to create boundary elements
  constexpr bool genBoundEls{true};
  
  return TPZGeoMeshTools::CreateGeoMeshOnGrid(dim,minX,maxX,matIdVec,nDivs,
                                              meshType,genBoundEls);
}