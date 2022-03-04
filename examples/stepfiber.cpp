/**
stepfiber.cpp

This target performs the modal analysis of a step-index optical fiber.

It also illustrates the usage of the SLEPc solver and how to create
complex curved structured meshes in NeoPZ.
***/

//wgma includes
#include "wganalysis.hpp"
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
#include "slepcepshandler.hpp"
//pz includes
#include <MMeshType.h>                   //for MMeshType
#include <pzcmesh.h>                     //for TPZCompMesh
#include <pzgmesh.h>                     //for TPZGeoMesh
#include <pzlog.h>                       //for TPZLogger
#include <TPZElectromagneticConstants.h> //for pzelectromag::cZero
#include <TPZKrylovEigenSolver.h>        //for TPZKrylovEigenSolver





#include <TPZSimpleTimer.h>              //for TPZSimpleTimer
/**
   @brief Creates the geometrical mesh associated with a step-index optical fiber.
   @param[in] rCore radius of the core
   @param[in] boundDist distance from the core from which the PML begins
   @param[in] dPML pml width
   @param[in] nLayersPml number of layers in the pml
   @param[in] factor number of subdivisions used in which quadrilateral region of the mesh
   @param[in] refine whether to refine the elements near the core/cladding interface
   @param[in] scale scale used in the domain (for better floating point precision)
   @param[out] matIdVec vector with the material identifiers: core/cladding/pmls/bcs
   @param[out] pmlTypeVec vector with the pml types (in the same order as `matIdVec`)
*/
TPZAutoPointer<TPZGeoMesh>
CreateStepFiberMesh(
  const REAL rCore, const REAL boundDist, const int factor,
  bool refine, const REAL scale, TPZVec<int> &matIdVec,
  const bool usingPML, const REAL dPML, const int nLayersPML,
  TPZVec<wgma::pml::type> &pmlTypeVec
  );

TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const int neigenpairs, const CSTATE target,
            TPZEigenSort sorting, bool usingSLEPC);

int main(int argc, char *argv[]) {
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif

  /***********************
   * setting the problem *
   ***********************/

  //see CreateStepFiberMesh for details on parameters
  
  //refractive index of the fibers core
  constexpr STATE core_n{1.4457};
  //refractive index of the fibers cladding
  constexpr STATE cladding_n{1.4378};
  //magnetic permeability of the core
  constexpr STATE coreUr{1};
  //magnetic permeability of the cladding
  constexpr STATE claddingUr{1};
  //radius of the core
  constexpr REAL rCore{8.0e-6};
  // operational wavelength
  constexpr STATE lambda{1.55e-6};
  /*both distances from the core and pml width are measured
   in wavelengths in the cladding material*/
  constexpr REAL lambdaCladding = lambda*cladding_n;
  //whether to surround the domain by a PML
  constexpr bool usingPML{true};
  //distance from the core from which the PML begins
  constexpr REAL boundDist{2*lambdaCladding};
  //pml width
  constexpr REAL dPML{4*lambdaCladding};
  //number of layers in the pml
  constexpr int nLayersPML{2};
  /**PML attenuation constant
     for a quadratic PML, one suggestion is to calculate it as:
     \frac{-3 \ln (R)}{4 \omega d n}
     where
     R = desired theoretical reflection coefficient
     \omega = operational frequency
     d = pml thickness
     n = reffraction index of region next to pml.

     for this value:
     R = 10^{-25}
     \omega = 2 \pi c_0 / \lambda
     d = \lambda * 4
     n = 1.4378
   */
  constexpr STATE alphaPML{2.7e-9};
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
  constexpr int pOrder{2};

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
  constexpr CSTATE target = -2.09;


  /*********************
   * exporting options *
   *********************/

  //whether to print the geometric mesh in .txt and .vtk formats
  constexpr bool printGMesh{true};
  //whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  //resolution of the .vtk file in which the solution will be exported
  constexpr int vtkRes{0};
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

  /* Vector for storing materials(regions) identifiers for the 
     TPZGeoMeshTools::CreateGeoMeshOnGrid function.
     In this example, we have one material for the interior of the waveguide
     and one for each BC*/
  TPZManVector<int, 11> matIdVec;
  TPZManVector<wgma::pml::type,8> pmlTypeVec;
  constexpr int factor{3};
  constexpr bool refine{false};
  

  bool usingSLEPC{true};
  //creates gmesh
  auto gmesh = CreateStepFiberMesh(rCore,boundDist,factor,refine,scale,
                                   matIdVec,usingPML,dPML,nLayersPML,pmlTypeVec);

  //print gmesh to .txt and .vtk format
  if(printGMesh)
  {
    //prefix for the gmesh files
    const std::string prefix = refine ? "ref_" : "noref_";
    const auto suffix = usingPML? "wpml" : "";
    const auto filename = "stepfiber_gmesh_"+
      std::to_string(factor) + " _" + prefix + suffix;
    wgma::gmeshtools::PrintGeoMesh(gmesh,filename);
  }

  //setting up cmesh data
  //Let us fill the struct needed to CMeshWgma2D
  wgma::cmeshtools::PhysicalData data;
  
  {
    int matCount = 0;
    data.matinfovec.push_back(
      std::make_tuple(matIdVec[matCount++],
                      core_n*core_n,
                      coreUr*coreUr));
    data.matinfovec.push_back(
      std::make_tuple(matIdVec[matCount++],
                      cladding_n*cladding_n,
                      claddingUr*claddingUr));
    
    constexpr int nPML = usingPML ? 8 : 0;

    for(int ipml = 0; ipml < nPML; ipml++){
      wgma::pml::data pml;
      pml.id = matIdVec[matCount++];
      pml.alphax = alphaPML;
      pml.alphay = alphaPML;
      pml.t = pmlTypeVec[ipml];
      data.pmlvec.push_back(pml);
    }
    
    wgma::bc::data bc;
    bc.id = matIdVec[matCount];
    bc.t = wgma::bc::type::PEC;

    data.bcvec.push_back(bc);
  }

  /*
   The problem uses an H1 approximation space for the longitudinal component 
   and a HCurl approximation space for the transversal one. Therefore, three
  // computational meshes are generated. One for each space and a multiphysics mesh*/
  auto meshVec = wgma::cmeshtools::CMeshWgma2D(gmesh,pOrder,data, lambda,scale);

  
  //WGAnalysis class is responsible for managing the modal analysis
  wgma::WGAnalysis analysis(meshVec,nThreads,optimizeBandwidth,filterBoundaryEqs);
  
  auto solver = SetupSolver(nEigenpairs, target, sortingRule, usingSLEPC);

  analysis.SetSolver(solver);
  analysis.Run(computeVectors);
  
  if (!computeVectors && !exportVtk) return 0;

  const std::string plotfile = "stepfiber_field_";

  TPZSimpleTimer postProc("Post processing");
  
  analysis.PostProcess(plotfile, vtkRes, printRealPart);
  return 0;
}




TPZAutoPointer<TPZGeoMesh>
CreateStepFiberMesh(
  const REAL realRCore, const REAL boundDist, const int factor,
  bool refine, const REAL scale,  TPZVec<int> &matIdVec,
  const bool usingPML, const REAL realDPML, const int nLayersPML,
  TPZVec<wgma::pml::type> &pmlTypeVec
  )
{

  
  using wgma::gmeshtools::EdgeData;
  using wgma::gmeshtools::QuadData;
  using wgma::gmeshtools::CircArcMap;
  
  const int nElsCoreR = factor * 3,//number of points in core (radial direction)
    nElsCoreT = factor * 4,//number of points in core (tangential direction)
    nElsCladdingR = factor * 2,
    nElsPml = factor * nLayersPML + 1;

  if(std::min<int>({nElsCoreR,nElsCoreT,nElsCladdingR,nElsPml}) < 2 ) {
    std::cout<<"Mesh has not sufficient divisions."<<std::endl;
    std::cout<<"Minimum is 2."<<std::endl;
    DebugStop();
  }

  constexpr int nPtsCore = 8;
  constexpr int nPtsCladding = 4;
  const int nPtsPML = usingPML ? 12 : 0;
  const int nPoints = nPtsCore + nPtsCladding + nPtsPML;
  constexpr int maxNPts = 24;

  constexpr int nQuadsCore = 5;
  constexpr int nQuadsCladding = 4;
  const int nQuadsPML = usingPML ? 8 : 0;
  const int nQuads = nQuadsCore + nQuadsCladding + nQuadsPML;
  constexpr int maxNQuads = 17;

  
  constexpr int matIdCore = 1, matIdCladding = 2, matIdBC= 18;
  constexpr int matIdPMLxp = 10,
    matIdPMLyp = 11,
    matIdPMLxm = 12,
    matIdPMLym = 13,
    matIdPMLxpym = 14,
    matIdPMLxpyp = 15,
    matIdPMLxmyp = 16,
    matIdPMLxmym = 17;

  const REAL rCore = realRCore / scale;
  const REAL bound = rCore + boundDist / scale;
  const REAL dPML = realDPML / scale;
  TPZManVector<REAL,2> xc(3, 0.);
  xc[0] = 0.;
  xc[1] = 0.;



  ///all points
  TPZManVector<TPZVec<REAL>,maxNPts> pointsVec(nPoints,TPZVec<REAL>(3,0.));
  //inner part of the core (lower right, upper right, upper left, lower left)
  pointsVec[0][0]= xc[0] + M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[0][1]= xc[1] - M_SQRT1_2 * std::sin(M_PI_4)*rCore;
  pointsVec[1][0]= xc[0] + M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[1][1]= xc[1] + M_SQRT1_2 * std::sin(M_PI_4)*rCore;
  pointsVec[2][0]= xc[0] - M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[2][1]= xc[1] + M_SQRT1_2 * std::sin(M_PI_4)*rCore;
  pointsVec[3][0]= xc[0] - M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[3][1]= xc[1] - M_SQRT1_2 * std::sin(M_PI_4)*rCore;
  //ok
  //outer part of the core (lower right, upper right, upper left, lower left)
  pointsVec[4][0]= xc[0] + std::cos(M_PI_4)*rCore; pointsVec[4][1]= xc[1] - std::sin(M_PI_4)*rCore;
  pointsVec[5][0]= xc[0] + std::cos(M_PI_4)*rCore; pointsVec[5][1]= xc[1] + std::sin(M_PI_4)*rCore;
  pointsVec[6][0]= xc[0] - std::cos(M_PI_4)*rCore; pointsVec[6][1]= xc[1] + std::sin(M_PI_4)*rCore;
  pointsVec[7][0]= xc[0] - std::cos(M_PI_4)*rCore; pointsVec[7][1]= xc[1] - std::sin(M_PI_4)*rCore;
  //ok

  //cladding (lower right, upper right, upper left, lower left)
  pointsVec[8][0] =  1 * bound; pointsVec[8][1] = -1 * bound;
  pointsVec[9][0] =  1 * bound; pointsVec[9][1] =  1 * bound;
  pointsVec[10][0] = -1 * bound; pointsVec[10][1] =  1 * bound;
  pointsVec[11][0] = -1 * bound; pointsVec[11][1] = -1 * bound;
  //ok
  if(usingPML){
    //xp
    pointsVec[12][0] =  1 * (bound+dPML); pointsVec[12][1] = -1 * bound;
    pointsVec[13][0] =  1 * (bound+dPML); pointsVec[13][1] =  1 * bound;
    //yp
    pointsVec[14][0] =  1 * bound; pointsVec[14][1] =  1 * (bound+dPML);
    pointsVec[15][0] = -1 * bound; pointsVec[15][1] =  1 * (bound+dPML);
    //xm
    pointsVec[16][0] = -1 * (bound+dPML); pointsVec[16][1] =  1 * bound;
    pointsVec[17][0] = -1 * (bound+dPML); pointsVec[17][1] = -1 * bound;
    //ym
    pointsVec[18][0] = -1 * bound; pointsVec[18][1] = -1 * (bound+dPML);
    pointsVec[19][0] = 1 * bound; pointsVec[19][1] = -1 * (bound+dPML);
    //xpym
    pointsVec[20][0] =  1 * (bound+dPML); pointsVec[20][1] = -1 * (bound+dPML);
    //xpyp
    pointsVec[21][0] =  1 * (bound+dPML); pointsVec[21][1] =  1 * (bound+dPML);
    //xmyp
    pointsVec[22][0] = -1 * (bound+dPML); pointsVec[22][1] =  1 * (bound+dPML);
    //xmym
    pointsVec[23][0] = -1 * (bound+dPML); pointsVec[23][1] = -1 * (bound+dPML);
  }

  //CREATING QUADS
  TPZManVector<QuadData,maxNQuads> quadVec(nQuads);
  //core
  quadVec[0].m_nodes ={ 0, 1, 2, 3};
  quadVec[1].m_nodes ={ 0, 4, 5, 1};
  quadVec[2].m_nodes ={ 2, 1, 5, 6};
  quadVec[3].m_nodes ={ 7, 3, 2, 6};
  quadVec[4].m_nodes ={ 7, 4, 0, 3};
  for(int i = 0; i < nQuadsCore; i++){
    quadVec[i].m_matid = matIdCore;
  }
  //cladding
  quadVec[5].m_nodes ={ 4, 8, 9, 5};
  quadVec[6].m_nodes ={ 6, 5, 9, 10};
  quadVec[7].m_nodes ={ 11, 7, 6, 10};
  quadVec[8].m_nodes ={ 11, 8, 4, 7};
  for(int i = 0; i < nQuadsCladding; i++){
    quadVec[5+i].m_matid = matIdCladding;
  }
  //pml
  if(usingPML){
    quadVec[9].m_nodes ={ 13, 9, 8, 12};
    quadVec[9].m_matid = matIdPMLxp;
    quadVec[10].m_nodes ={ 14, 15, 10, 9};
    quadVec[10].m_matid = matIdPMLyp;
    quadVec[11].m_nodes ={ 10, 16, 17, 11};
    quadVec[11].m_matid = matIdPMLxm;
    quadVec[12].m_nodes ={ 8, 11, 18, 19};
    quadVec[12].m_matid = matIdPMLym;
    quadVec[13].m_nodes ={ 21, 14,  9, 13};
    quadVec[13].m_matid = matIdPMLxpyp;
    quadVec[14].m_nodes ={ 15, 22, 16, 10};
    quadVec[14].m_matid = matIdPMLxmyp;
    quadVec[15].m_nodes ={ 11, 17, 23, 18};
    quadVec[15].m_matid = matIdPMLxmym;
    quadVec[16].m_nodes ={ 12,  8, 19, 20};
    quadVec[16].m_matid = matIdPMLxpym;
  }

  for(int i = 0; i < nQuads; i++){
    quadVec[i].m_type = wgma::gmeshtools::ElType::Tri;
  }
  
  //CREATING EDGES
  constexpr int nEdgesCore = 12;
  constexpr int nEdgesCladding = 8;
  const int nEdgesPml = usingPML ? 20 : 0;
  const int nEdges = nEdgesCore + nEdgesCladding + nEdgesPml;
  constexpr int maxNEdges = 40;

  constexpr int matIdInterface{-15};
  TPZManVector<EdgeData,maxNEdges> edgeVec(nEdges);

  //inner edges of the core
  int iedge = 0;
  for(int i = 0; i < 4; i++){
    edgeVec[iedge].m_nel = nElsCoreT;
    edgeVec[iedge].m_create_el = false;
    edgeVec[iedge].m_map = nullptr;
    edgeVec[iedge].m_nodes = {i, (i+1)%4};
    iedge++;
  }
  //radial edges of the core
  for(int i = 0; i < 4; i++){
    edgeVec[iedge].m_nel = nElsCoreR;
    edgeVec[iedge].m_create_el = false;
    edgeVec[iedge].m_map = nullptr;
    edgeVec[iedge].m_nodes = {i, i+4};
    iedge++;
  }
  //outer edges of the core, curved edges
  for(int i = 0; i<4;i++){
    const auto theta_i = -M_PI/4 + i * M_PI/2;
    const auto theta_f = M_PI/4 + i * M_PI/2;
    edgeVec[iedge].m_matid = matIdInterface;//arbitrary
    edgeVec[iedge].m_nel = nElsCoreT;
    edgeVec[iedge].m_create_el = true;
    edgeVec[iedge].m_nodes = {4+i,4+(i+1)%4};//({7,4},{4,5},{5,6},{6,7})
    edgeVec[iedge].m_map = [rCore, xc,theta_i,theta_f](const REAL s){
      return CircArcMap({theta_i,theta_f},xc,rCore,s);};
    iedge++;
  }
  //edges connecting the core to the boundary
  for(int i = 0; i < 4; i++){
    edgeVec[iedge].m_nel = nElsCladdingR;
    edgeVec[iedge].m_nodes = {4+i, 8+i};
    iedge++;
  }

  //interface cladding-pml(or cladding-boundary)
  for(int i = 0; i < 4; i++){
    edgeVec[iedge].m_nel = nElsCoreT;
    edgeVec[iedge].m_create_el = usingPML ? false : true;
    edgeVec[iedge].m_matid = matIdBC;//ignored if usingPML==true
    edgeVec[iedge].m_nodes = {8+i, 8+(i+1)%4};
    iedge++;
  }

  if(usingPML){
    //pml xp,yp,xm,ym
    for(int i = 0; i < 4; i++){
      edgeVec[iedge].m_nel = nElsPml;
      edgeVec[iedge].m_nodes = {8+i, 12+2*i};//ok
      iedge++;
      edgeVec[iedge].m_nel = nElsPml;
      edgeVec[iedge].m_nodes = {8+(i+1)%4,13+2*i};//ok
      iedge++;
      edgeVec[iedge].m_nel = nElsCoreT;
      edgeVec[iedge].m_nodes = {12+2*i,13+2*i};//ok
      edgeVec[iedge].m_create_el = true;
      edgeVec[iedge].m_matid = matIdBC;
      iedge++;
    }
    //pml xpym, xpyp, xmyp, xmym
    for(int i = 0; i < 4; i++){
      edgeVec[iedge].m_nel = nElsPml;
      //({19,20},{13,21},{15,22},{17,23})
      edgeVec[iedge].m_nodes = {13+2*((i+3)%4), 20+i};
      edgeVec[iedge].m_create_el = true;
      edgeVec[iedge].m_matid = matIdBC;
      iedge++;
      edgeVec[iedge].m_nel = nElsPml;
      //({20,12},{21,14},{22,16},{23,18})
      edgeVec[iedge].m_nodes = {20+i,12+2*i};
      edgeVec[iedge].m_create_el = true;
      edgeVec[iedge].m_matid = matIdBC;
      iedge++;
    }
  }
  
  
  const bool nonLinearMapping{true};
  auto gmesh = wgma::gmeshtools::CreateStructuredMesh(
    pointsVec, quadVec, edgeVec, nonLinearMapping);

  //the materials are: core, cladding, pml regions (8) and bc
  const int nMats = usingPML ? 11 : 3;
  matIdVec.Resize(nMats);
  matIdVec[0]=matIdCore;
  matIdVec[1]=matIdCladding;
  if(usingPML){
    pmlTypeVec.Resize(8);
    pmlTypeVec[0]=wgma::pml::type::xp;
    matIdVec[2]=matIdPMLxp;
    pmlTypeVec[1]=wgma::pml::type::yp;
    matIdVec[3]=matIdPMLyp;
    pmlTypeVec[2]=wgma::pml::type::xm;
    matIdVec[4]=matIdPMLxm;
    pmlTypeVec[3]=wgma::pml::type::ym;
    matIdVec[5]=matIdPMLym;
    pmlTypeVec[4]=wgma::pml::type::xpym;
    matIdVec[6]=matIdPMLxpym;
    pmlTypeVec[5]=wgma::pml::type::xpyp;
    matIdVec[7]=matIdPMLxpyp;
    pmlTypeVec[6]=wgma::pml::type::xmyp;
    matIdVec[8]=matIdPMLxmyp;
    pmlTypeVec[7]=wgma::pml::type::xmym;
    matIdVec[9]=matIdPMLxmym;
  }
  
  matIdVec[nMats-1]=matIdBC;

  gmesh->BuildConnectivity();

  if(refine){
    std::set<int> elsToRefine;
    for(const auto gel : gmesh->ElementVec()){
      if(gel->MaterialId() == matIdInterface){
        elsToRefine.insert(gel->Id());
        for(int iv = 0; iv < 2; iv++){
          TPZGeoElSide gelside(gel,iv);
          TPZGeoElSide neighbour = gelside.Neighbour();
          while (neighbour != gelside) {
            const auto neighgel = neighbour.Element();
            if(neighgel->Dimension() == 2){
              elsToRefine.insert(neighgel->Id());
            }
            neighbour = neighbour.Neighbour();
          }
        }
      }
    }
    TPZManVector<TPZGeoEl *,4>sons;
    for(const auto id : elsToRefine){
      auto gel = gmesh->Element(id);
      //just to be on the safe side
      if(!gel->HasSubElement()){
        gel->Divide(sons);
      }
    }
  }
  gmesh->BuildConnectivity();

  return gmesh;
}


TPZAutoPointer<TPZEigenSolver<CSTATE>>
SetupSolver(const int neigenpairs, const CSTATE target,
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

  constexpr int krylovDim{-1};
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
    
    constexpr Precond pc = Precond::LU;
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