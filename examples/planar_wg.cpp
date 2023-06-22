/**
stepfiber.cpp

This target performs the modal analysis of a double-ridged waveguide.

It also illustrates how to import a .gmsh mesh and how to perform
directional mesh refinement.
***/

//wgma includes
#include <scattering.hpp>
#include "cmeshtools.hpp"
#include "cmeshtools_impl.hpp"//for custom pml regions
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include <util.hpp>
#include "slepcepshandler.hpp"
//pz includes
#include <MMeshType.h>                   //for MMeshType
#include <pzcmesh.h>                     //for TPZCompMesh
#include <pzgmesh.h>                     //for TPZGeoMesh
#include <pzlog.h>                       //for TPZLogger
#include <TPZSimpleTimer.h>              //for TPZSimpleTimer
#include <TPZVTKGenerator.h>


/**
   @brief Implements a source function to be inserted in the problem.
   @note This function corresponds to the fundamental TE mode calculated at
   lambda = 1.55um
   strip width = 3um
   n strip = 1.55
   n air = 1
*/
std::pair<wgma::scattering::srcfunc1d, CSTATE>
GetSourceFunc(const REAL scale, const REAL lambda);

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
  matmap["core"] = std::make_pair<CSTATE,CSTATE>(1.55*1.55,1.);
  std::map<std::string,wgma::bc::type> bcmap;
  bcmap["bound"] = wgma::bc::type::PEC;
  std::set<std::string> srclist;
  srclist.insert("gamma_1");
  srclist.insert("gamma_2");
  srclist.insert("gamma_3");
  srclist.insert("gamma_4");
  srclist.insert("gamma_5");
  // operational wavelength
  constexpr STATE lambda{1.55e-6};
  /*
    Given the small dimensions of the domain, scaling it can help in 
    achieving good precision. Using 1./k0 as a scale factor results in 
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2.
    This scale factor is often referred to as characteristic length
    of the domain.
  */
  constexpr REAL scale{lambda/(2*M_PI)};

  constexpr wgma::planarwg::mode mode{wgma::planarwg::mode::TE};

  /******************
   *  fem options   *
   ******************/
  // polynomial order to be used in the approximation
  constexpr int pOrder{1};

  /******************
   * solver options *
   ******************/
  
  //number of threads to use
  constexpr int nThreads{8};

  //PML attenuation constant
  constexpr STATE alphaxPML {10.00};
  constexpr STATE alphayPML = alphaxPML/2;

  /*********************
   * exporting options *
   *********************/

  //whether to print the geometric mesh in .txt and .vtk formats
  constexpr bool printGMesh{true};
  //whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  // path for output files
  const std::string path {"res_planar_wg/"};
  // common prefix for both meshes and output files
  const std::string basisName{"planar_wg"};
  // prefix for exported files
  const std::string prefix{path+basisName};
  //just to make sure we will output results
  wgma::util::CreatePath(wgma::util::ExtractPath(prefix));
  //resolution of the .vtk file in which the solution will be exported
  constexpr int vtkRes{0};


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
  //check meshes/planar_wg.msh for a simpler example of the TE mode
  const std::string filename{"meshes/wg_disc.msh"};
  TPZVec<std::map<std::string,int>> gmshmats;
  auto gmesh = wgma::gmeshtools::ReadGmshMesh(filename, scale, gmshmats);

  
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
                                          {alphaxPML, alphayPML},
                                          data);


  
  
  std::set<int> source_ids;
  {
    //let us iterate through all 1D materials read by gmsh
    for(auto gmshbc : gmshmats[1]){
      const auto bcname = gmshbc.first;
      if(srclist.find(bcname) != srclist.end()){
        source_ids.insert(gmshbc.second);
      }
    }
  }

  
    
  auto cmesh =
    wgma::scattering::CMeshScattering2D(gmesh, mode, pOrder,
                                        data, source_ids, lambda,scale);

  auto mysource = GetSourceFunc(scale,lambda);

  auto beta = mysource.second;
  
  wgma::scattering::Source1D sources;
  sources.func = mysource.first;
  sources.id = source_ids;
  wgma::scattering::LoadSource1D(cmesh, sources);
  wgma::scattering::SetPropagationConstant(cmesh, beta);

  auto an = wgma::scattering::Analysis(cmesh, nThreads,
                                       optimizeBandwidth, filterBoundaryEqs);

  an.Run();

  const std::string plotfile = prefix+".vtk";
  TPZSimpleTimer postProc("Post processing");

  TPZVec<std::string> fields = {
    "Field_real",
    "Field_imag",
    "Field_abs"};
  auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
  vtk.Do();
}


std::pair<wgma::scattering::srcfunc1d, CSTATE>
GetSourceFunc(const REAL scale, const REAL lambda)
{
  /*
    The following source was calculated for
    lambda = 1.55um
    strip width = 3um
    n strip = 1.55
    n air = 1
  */
  const STATE k0 = scale * 2*M_PI/lambda;
  //calculated effective index
  const STATE nev = 1.20078029226310301069702290988;
  const STATE n1 = 1.55;
  const STATE n2 = 1;
  const STATE kx = k0*sqrt(n1*n1 - nev*nev);
  const STATE gammax = k0*sqrt(nev*nev- n2*n2);
  const REAL d = 0.15 * 1e-6 / scale;
  const STATE beta = k0 * nev;
  
  
  wgma::scattering::srcfunc1d sourcefunc =
    [d, kx, gammax, beta](const TPZVec<REAL> &loc,
                          CSTATE &val,
                          TPZVec<CSTATE> &deriv){
      const auto x = loc[1];
      const auto xabs = std::abs(x);
      if(xabs < d){
        val = cos(kx*x);
      }else{
        val = cos(kx*d)*exp(-gammax*(xabs-d));
      }
      deriv = {0,0,0};
    };
  return std::make_pair(sourcefunc,beta);
}