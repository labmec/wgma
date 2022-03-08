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
#include "slepcepshandler.hpp"
//pz includes
#include <MMeshType.h>                   //for MMeshType
#include <pzcmesh.h>                     //for TPZCompMesh
#include <pzgmesh.h>                     //for TPZGeoMesh
#include <pzlog.h>                       //for TPZLogger
#include <TPZSimpleTimer.h>              //for TPZSimpleTimer



/**
   @brief Implements a source function to be inserted in the problem.
   @note This function corresponds to the fundamental TE mode calculated at
   lambda = 1.55um
   strip width = 3um
   n strip = 1.55
   n air = 1
*/
ForcingFunctionBCType<CSTATE> GetSourceFunc(const REAL scale, const REAL lambda);

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
  bcmap["gamma_1"] = wgma::bc::type::SOURCE;
  bcmap["gamma_2"] = wgma::bc::type::SOURCE;
  bcmap["gamma_3"] = wgma::bc::type::SOURCE;
  bcmap["gamma_4"] = wgma::bc::type::SOURCE;
  bcmap["gamma_5"] = wgma::bc::type::SOURCE;
  bcmap["bound"] = wgma::bc::type::PEC;
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
  //prefix for exported files
  const std::string prefix{"planar_wg"};
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
                                          alphaxPML, alphayPML,
                                          data);


  std::vector<wgma::bc::source> sources;

  {
    auto mysource = GetSourceFunc(scale,lambda);
    
    //let us iterate through all 1D materials read by gmsh
    for(auto gmshbc : gmshmats[1]){
      const auto bcname = gmshbc.first;
      if(bcmap.find(bcname) != bcmap.end()){
        const auto type = bcmap.at(bcname);
        if(type == wgma::bc::type::SOURCE){
          wgma::bc::source src;
          src.id = gmshbc.second;
          src.func = mysource;
          sources.push_back(src);
        }
      }
    }
  }
  
  auto cmesh =
    wgma::cmeshtools::CMeshScattering2D(gmesh, mode, pOrder,
                                        data, sources, lambda,scale);


  auto an = wgma::ScatteringAnalysis(cmesh, nThreads,
                                     optimizeBandwidth, filterBoundaryEqs);

  an.Run();

  const std::string plotfile = prefix+".vtk";
  an.PostProcess(plotfile, vtkRes);
}


ForcingFunctionBCType<CSTATE> GetSourceFunc(const REAL scale, const REAL lambda)
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
  
  
  ForcingFunctionBCType<CSTATE> sourcefunc =
    [d, kx, gammax, beta](const TPZVec<REAL> &loc,
                          TPZVec<CSTATE> &rhsVal,
                          TPZFMatrix<CSTATE> &matVal){
    
      const auto x = loc[1];
      const auto xabs = std::abs(x);
      if(xabs < d){
        rhsVal[0] = cos(kx*x);
      }else{
        rhsVal[0] = cos(kx*d)*exp(-gammax*(xabs-d));
      }
      rhsVal[1] = beta ;
    };
  return sourcefunc;
}