/**
stepfiber.cpp

This target performs the modal analysis of a double-ridged waveguide.

It also illustrates how to import a .gmsh mesh and how to perform
directional mesh refinement.
***/

//wgma includes
#include "wganalysis.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include "slepcepshandler.hpp"
#include "util.hpp"
//pz includes
#include <MMeshType.h>                   //for MMeshType
#include <pzcmesh.h>                     //for TPZCompMesh
#include <pzgmesh.h>                     //for TPZGeoMesh
#include <pzlog.h>                       //for TPZLogger
#include <TPZElectromagneticConstants.h> //for pzelectromag::cZero
#include <TPZKrylovEigenSolver.h>        //for TPZKrylovEigenSolver
#include <pzintel.h>
#include <pzmultiphysicselement.h>
#include <Electromagnetics/TPZWgma.h>
#include <pzbuildmultiphysicsmesh.h>
#include <TPZVTKGenerator.h>
#include <TPZSimpleTimer.h>              //for TPZSimpleTimer



int main(int argc, char *argv[]) {
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif

  /***********************
   * setting the problem *
   ***********************/

  std::map<std::string,std::pair<CSTATE,CSTATE>> matmap;
  matmap["air"] = std::make_pair<CSTATE,CSTATE>(1.,1.);
  matmap["core"] = std::make_pair<CSTATE,CSTATE>(3.44*3.44,1.);
  matmap["substrate"] = std::make_pair<CSTATE,CSTATE>(3.40*3.40,1.);
  std::map<std::string,wgma::bc::type> bcmap;
  bcmap["bound"] = wgma::bc::type::PEC;
  // operational wavelength
  constexpr STATE lambda{1.15e-6};
  /*
    Given the small dimensions of the domain, scaling it can help in 
    achieving good precision. Using 1./k0 as a scale factor results in 
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  constexpr REAL scale{lambda/(2*M_PI)};

  /******************
   *  fem options   *
   ******************/
  // polynomial order to be used in the approximation
  constexpr int pOrder{3};

  /******************
   * solver options *
   ******************/
  
  //number of threads to use
  constexpr int nThreads{8};
  //number of genvalues to be computed
  constexpr int nEigenpairs{4};
  //whether to compute eigenvectors (instead of just eigenvalues)
  constexpr bool computeVectors{true};
  //how to sort the computed eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  /*
   The simulation uses a Krylov-based Arnoldi solver for solving the
   generalised EVP. A shift-and-inverse spectral transform is applied in 
   the system for better convergence of the eigenvalues.
   The target variable should be close to the desired eigenvalue (in this case,
   the effective index neff).*/
  constexpr CSTATE target = -12;

  //PML attenuation constant
  constexpr STATE alphaPML{0.04};

  /*********************
   * exporting options *
   *********************/

  //whether to print the geometric mesh in .txt and .vtk formats
  constexpr bool printGMesh{true};
  //whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  //prefix for exported files
  const std::string prefix{"res_ribwg/ribwg"};
  wgma::util::CreatePath(wgma::util::ExtractPath(prefix));
  //resolution of the .vtk file in which the solution will be exported
  constexpr int vtkRes{3};
  //if true, the real part of the electric fields is exported. otherwise, the magnitude
  constexpr bool printRealPart{true};

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
  

  bool usingSLEPC{true};
  //creates gmesh
  //filename
  const std::string filename{"meshes/ribWG.msh"};
  TPZVec<std::map<std::string,int>> gmshmats;
  auto gmesh = wgma::gmeshtools::ReadGmshMesh(filename, scale, gmshmats);

  /*
    Since we know that there is a singularity at the metallic corner,
    we will directionally refine all the elements with a node in this corner.
   */
  std::set<int> refids;
  refids.insert(gmshmats[0].at("corner"));
  constexpr int nrefdir{5};
  wgma::gmeshtools::DirectionalRefinement(gmesh, refids, nrefdir);

  
  //print gmesh to .txt and .vtk format
  if(printGMesh)
  {
    //prefix for the gmesh files
    const std::string filename = prefix +"_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh,filename);
  }

  
  //setting up cmesh data
  wgma::cmeshtools::PhysicalData data;
  wgma::cmeshtools::SetupGmshMaterialData(gmshmats, matmap, bcmap,
                                          {alphaPML, alphaPML},
                                          data);
  /*
   The problem uses an H1 approximation space for the longitudinal component 
   and a HCurl approximation space for the transversal one. Therefore, three
  // computational meshes are generated. One for each space and a multiphysics mesh*/
  auto meshVec = wgma::wganalysis::CMeshWgma2D(gmesh,pOrder,data, lambda,scale);

  /**
     Increasing the polynomial order of an element can only improve
     the quality of the FEM solution as long as the solution is regular
     enough. Therefore, we lower the polynomial order of the elements
     adjacent to the singularities.
   */
  [](TPZVec<TPZAutoPointer<TPZCompMesh>> &meshVec, std::set<int> singids){
    auto mfmesh = meshVec[0];
    TPZVec<TPZCompMesh*> atomicmeshes = {meshVec[1].operator->(),
      meshVec[2].operator->()};

    constexpr int h1index = TPZWgma::H1Index();
    constexpr int hcurlindex = TPZWgma::HCurlIndex();
    /*
      In order to ensure De Rham compatibility, we need to have
      h1 order = hcurl order +1.
      given our choice of design in NeoPZ, minimum hcurl order is 1.
    **/
    constexpr int ordh1 = 2;
    constexpr int ordhcurl = 1;


    for(auto cel : mfmesh->ElementVec()){
      auto gel = cel->Reference();//geometric element
      /*
        Lets check the neighbours for each vertex
      */
      bool shouldrefine = false;
      for(int in=0; in<gel->NCornerNodes() && !shouldrefine; in++){
        TPZGeoElSide gels(gel,in);
        TPZGeoElSide neigh(gels.Neighbour());
        while(gels != neigh && !shouldrefine){
          if(singids.count(neigh.Element()->MaterialId())){
            shouldrefine = true;
          }
          neigh = neigh.Neighbour();
        }
      }
      if(shouldrefine){
        auto celmf = dynamic_cast<TPZMultiphysicsElement*>(cel);
        
        auto intelh1 =
            dynamic_cast<TPZInterpolatedElement*>(celmf->Element(h1index));
        intelh1->PRefine(ordh1);

        auto intelhcurl =
            dynamic_cast<TPZInterpolatedElement*>(celmf->Element(hcurlindex));
        intelhcurl->PRefine(ordhcurl);
      }
    }
    for (auto mesh : atomicmeshes){
      mesh->ExpandSolution();
    }
    TPZBuildMultiphysicsMesh::AddConnects(atomicmeshes, mfmesh.operator->());
    TPZBuildMultiphysicsMesh::TransferFromMeshes(atomicmeshes, mfmesh.operator->());
    mfmesh->ExpandSolution();
  }(meshVec, refids);
  
  //WGAnalysis class is responsible for managing the modal analysis
  wgma::wganalysis::Wgma2D analysis(meshVec,nThreads,optimizeBandwidth,filterBoundaryEqs);
  
  auto solver = wgma::wganalysis::SetupSolver(target,nEigenpairs, sortingRule, usingSLEPC);

  analysis.SetSolver(solver);
  analysis.Run(computeVectors);
  
  if (!computeVectors && !exportVtk) return 0;

  const std::string plotfile = prefix+"_field";

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

