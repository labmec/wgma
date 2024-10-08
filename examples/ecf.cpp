/**
ecf.cpp

This target performs the modal analysis of an Exposed Core Fiber (ECF)
and generates dispersion curves

***/

//wgma includes
#include "wganalysis.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include "slepcepshandler.hpp"
#include <util.hpp>
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
#include <TPZSimpleTimer.h>              //for TPZSimpleTimer
#include <TPZVTKGenerator.h>


//Sets geometric info regarding all the circles in the mesh
TPZVec<wgma::gmeshtools::ArcData> SetUpArcData(const REAL scale);

void RunSimulation(const STATE lambda, const int nEigenpairs,
                   const int pOrder, const std::string& prefix,
                   const bool printGMesh, const bool exportVtk,
                   const bool exportCsv, bool &usingSLEPC);

int main(int argc, char *argv[]) {
  
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif
  
  // polynomial order to be used in the approximation
  constexpr int pOrder{1};
  //number of genvalues to be computed
  constexpr int nEigenpairs{10};
  //whether to print the geometric mesh in .txt and .vtk formats
  constexpr bool printGMesh{false};
  //whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  //whether to export the eigenvalues in .csv format
  constexpr bool exportCsv{true};
  // path for output files
  const std::string path {"res_ecf/"};
  // common prefix for both meshes and output files
  const std::string basisName{"ecf"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(prefix));
  //whether to use SLEPC
  bool usingSLEPC{true};

  
  constexpr int nlambdas{1};
  constexpr STATE minlambda{5e-7};
  constexpr STATE maxlambda{5e-6};
  constexpr STATE deltalambda{(maxlambda-minlambda)/nlambdas};
  TPZSimpleTimer total("Dispersion curve");
  for(int il = 0; il < nlambdas; il++){
    std::cout<<"Running round "<<il+1<<" out of "<<nlambdas<<std::endl;
    const STATE lambda = minlambda + il * deltalambda;
    RunSimulation(lambda, nEigenpairs, pOrder, prefix, printGMesh, exportVtk,
                  exportCsv, usingSLEPC);
  }
  
  if(usingSLEPC){
    wgma::slepc::EPSHandler<CSTATE>::FinalizeSLEPc();
  }
  return 0;
}

void RunSimulation(const STATE lambda, const int nEigenpairs, const int pOrder, const std::string&prefix,
                   const bool printGMesh, const bool exportVtk,
                   const bool exportCsv, bool &usingSLEPC){

  
  std::map<std::string,std::pair<CSTATE,CSTATE>> matmap;
  matmap["air"] = std::make_pair<CSTATE,CSTATE>(1.,1.);
  matmap["core"] = std::make_pair<CSTATE,CSTATE>(1.45*1.45,1.);
  std::map<std::string,wgma::bc::type> bcmap;
  bcmap["bound"] = wgma::bc::type::PEC;
  // operational wavelength
  
  /*
    Given the small dimensions of the domain, scaling it can help in 
    achieving good precision. Using 1./k0 as a scale factor results in 
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  const REAL scale{lambda/(2*M_PI)};

  /******************
   * solver options *
   ******************/
  
  //number of threads to use
  constexpr int nThreads{8};
  //whether to compute eigenvectors (instead of just eigenvalues)
  constexpr bool computeVectors{true};
  //how to sort the computed eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::RealAscending};
  /*
   The simulation uses a Krylov-based Arnoldi solver for solving the
   generalised EVP. A shift-and-inverse spectral transform is applied in 
   the system for better convergence of the eigenvalues.
   The target variable should be close to the desired eigenvalue (in this case,
   the effective index neff).*/
  constexpr CSTATE target = -2.1025;

  //PML attenuation constant
  constexpr STATE alphaPML{0.02};
  //resolution of the .vtk file in which the solution will be exported
  constexpr int vtkRes{1};
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
  
  //creates gmesh
  //filename
  const std::string filename{"meshes/ecf.msh"};
  TPZVec<std::map<std::string,int>> gmshmats;
  auto gmesh = wgma::gmeshtools::ReadGmshMesh(filename, scale, gmshmats);

  /*
    In order to represent curved geometries correctly, we need to use
    TPZArc3D, which uses an exact mapping to represent circumference arcs.
    The neighbouring elements are also modified to take this deformation into account.
   */
  auto arcdata = SetUpArcData(scale);

  wgma::gmeshtools::SetExactArcRepresentation(gmesh, arcdata);
  //print gmesh to .txt and .vtk format
  if(printGMesh)
  {
    //prefix for the gmesh files
    const std::string filename = prefix +"_gmesh";
    wgma::gmeshtools::PrintGeoMesh(gmesh,filename);
  }

  wgma::cmeshtools::PhysicalData data;
  wgma::cmeshtools::SetupGmshMaterialData(gmshmats, matmap, bcmap,
                                          {alphaPML, alphaPML},
                                          data);
  /*
   The problem uses an H1 approximation space for the longitudinal component 
   and a HCurl approximation space for the transversal one. Therefore, three
  // computational meshes are generated. One for each space and a multiphysics mesh*/
  auto meshVec = wgma::wganalysis::CMeshWgma2D(gmesh,pOrder,data, lambda,scale);

  
  //Analysis class is responsible for managing the modal analysis
  wgma::wganalysis::Wgma2D analysis(meshVec,nThreads,
                                            optimizeBandwidth,filterBoundaryEqs);
  
  auto solver = wgma::wganalysis::SetupSolver(target,nEigenpairs, sortingRule, usingSLEPC);

  analysis.SetSolver(*solver);
  analysis.Run(computeVectors);

  if(exportCsv){
    const std::string csvfile = prefix+"_eigenvalues.csv";
    analysis.WriteToCsv(csvfile, lambda);
  }
  
  if (!computeVectors || !exportVtk) return;

  const std::string plotfile = prefix+"_field_";

  {
    
    TPZSimpleTimer tpostprocess("Post processing(new)");
    TPZVec<std::string> fvars = {
      "Ez_real",
      "Ez_abs",
      "Et_real",
      "Et_abs"};
    auto vtk = TPZVTKGenerator(meshVec[0], fvars, plotfile, vtkRes);
    auto ev = analysis.GetEigenvalues();
    for (int isol = 0; isol < ev.size(); isol++) {
      auto currentKz = std::sqrt(-1.0*ev[isol]);
      std::cout<<"\rPost processing step "<<isol+1<<" out of "<<ev.size()
               <<"(kz = "<<currentKz<<")"<<std::flush;
      analysis.LoadSolution(isol);
      //since only the solution has changed (no element/geometry change)
      //we dont need to compute the vtk points again, just the field values
      vtk.Do();
    }
  }
  return;
}



TPZVec<wgma::gmeshtools::ArcData> SetUpArcData(const REAL scale)
{
  TPZVec<wgma::gmeshtools::ArcData> arcdata(7);

  constexpr REAL um = 1e-6;
  const REAL r_core{14*um / scale};
  const REAL r_clad{(12.5*um) / scale}, r_clad_out{(100*um) / scale};

  //rotates a pt around the origin in the xy-plane
  auto RotatePt = [] (const REAL theta, const REAL x, const REAL y){
    const REAL xrot = cos(theta) * x - sin(theta) * y;
    const REAL yrot = sin(theta) * x + cos(theta) * y;
    return std::pair(xrot, yrot);
  };

  const REAL xini{1.05*r_clad+r_core}, yini{0};
  constexpr int ncircs{6};
  for(int ic = 0; ic < ncircs; ic++){
    const REAL theta = ic * M_PI/3;
    REAL xc, yc;
    std::tie(xc, yc) = RotatePt(theta, xini,yini);
    arcdata[ic].m_matid = 20 + ic+1;
    arcdata[ic].m_radius = r_clad;
    arcdata[ic].m_xc = xc;
    arcdata[ic].m_yc = yc;
    arcdata[ic].m_zc = 0;
  }

  arcdata[6].m_matid = 27;
  arcdata[6].m_radius = r_clad_out;
  arcdata[6].m_xc = 0;
  arcdata[6].m_yc = 0;
  arcdata[6].m_zc = 0;
  
  return arcdata;
}