/**
stepfiber.cpp

This target performs the modal analysis of a step-index optical fiber.

A HCurl-conforming approximation space is used for the transverse field components and
a H1-conforming approximation space for the axial component.
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
  constexpr STATE coreReffIndex{1.4457};
  //refractive index of the fibers cladding
  constexpr STATE claddingReffIndex{1.4378};
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
  constexpr REAL lambdaCladding = lambda*claddingReffIndex;
  //whether to surround the domain by a PML
  constexpr bool usingPML{true};
  //distance from the core from which the PML begins
  constexpr REAL boundDist{1*lambdaCladding};
  //pml width
  constexpr REAL dPML{4*lambdaCladding};
  //number of layers in the pml
  constexpr int nLayersPML{2};
  //PML attenuation constant
  constexpr STATE alphaPML{0.5};
  /*Given the small dimensions of the domain, scaling it can help in 
    achieving good precision. Uing k0 as a scale factor results in 
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2*/
  constexpr REAL scale{2*M_PI/lambda};

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
  constexpr bool printGMesh{false};
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
  TPZVec<int> volMatIdVec(2,-1);
  TPZVec<wgma::pml::data> pmlDataVec(usingPML ? 8 : 0);
  TPZVec<wgma::bc::data> bcDataVec(1);
  {
    int matCount = 0;
    for(auto &matid: volMatIdVec){
      matid = matIdVec[matCount++];
    }

    int pmlCount{0};
    for(auto &pml : pmlDataVec){
      pml.id = matIdVec[matCount++];
      pml.alpha = alphaPML;
      pml.t = pmlTypeVec[pmlCount++];
    }

    bcDataVec[0].id = matIdVec[matCount];
    bcDataVec[0].t = wgma::bc::type::PEC;
  }
  
  TPZVec<CSTATE> urVec({coreUr, claddingUr});
  TPZVec<CSTATE> erVec({coreReffIndex*coreReffIndex,
      claddingReffIndex*claddingReffIndex});

  /*
   The problem uses an H1 approximation space for the longitudinal component 
   and a HCurl approximation space for the transversal one. Therefore, three
  // computational meshes are generated. One for each space and a multiphysics mesh*/
  auto meshVec = wgma::cmeshtools::CreateCMesh(gmesh,pOrder,volMatIdVec,
                                               urVec, erVec, pmlDataVec,
                                               bcDataVec, lambda,scale);

  
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
  using wgma::gmeshtools::QuadrilateralData;
  const int nPtsCoreR = factor * 2,//number of points in core (radial direction)
    nPtsCoreT = factor * 5,//number of points in core (tangential direction)
    nPtsCladdingR = factor * 4,
    nDivPml = factor * nLayersPML + 1;

  if(std::min<int>({nPtsCoreR,nPtsCoreT,nPtsCladdingR,nDivPml}) < 2 ) {
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

  const REAL rCore = realRCore * scale;
  const REAL bound = rCore + boundDist * scale;
  const REAL dPML = realDPML * scale;
  TPZManVector<REAL,2> xc(2, 0.);
  xc[0] = 0.;
  xc[1] = 0.;



  ///all points
  TPZManVector<TPZVec<REAL>,maxNPts> pointsVec(nPoints,TPZVec<REAL>(2,0.));
  //hole 1
  pointsVec[0][0]= xc[0] + M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[0][1]= xc[1] + M_SQRT1_2 * std::sin(M_PI_4)*rCore;
  pointsVec[1][0]= xc[0] - M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[1][1]= xc[1] + M_SQRT1_2 * std::sin(M_PI_4)*rCore;
  pointsVec[2][0]= xc[0] - M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[2][1]= xc[1] - M_SQRT1_2 * std::sin(M_PI_4)*rCore;
  pointsVec[3][0]= xc[0] + M_SQRT1_2 * std::cos(M_PI_4)*rCore; pointsVec[3][1]= xc[1] - M_SQRT1_2 * std::sin(M_PI_4)*rCore;

  pointsVec[4][0]= xc[0] + std::cos(M_PI_4)*rCore; pointsVec[4][1]= xc[1] + std::sin(M_PI_4)*rCore;
  pointsVec[5][0]= xc[0] - std::cos(M_PI_4)*rCore; pointsVec[5][1]= xc[1] + std::sin(M_PI_4)*rCore;
  pointsVec[6][0]= xc[0] - std::cos(M_PI_4)*rCore; pointsVec[6][1]= xc[1] - std::sin(M_PI_4)*rCore;
  pointsVec[7][0]= xc[0] + std::cos(M_PI_4)*rCore; pointsVec[7][1]= xc[1] - std::sin(M_PI_4)*rCore;


  //cladding
  pointsVec[8][0] =  1 * bound; pointsVec[8][1] = -1 * bound;
  pointsVec[9][0] =  1 * bound; pointsVec[9][1] =  1 * bound;
  pointsVec[10][0] = -1 * bound; pointsVec[10][1] =  1 * bound;
  pointsVec[11][0] = -1 * bound; pointsVec[11][1] = -1 * bound;
  if(usingPML){
    pointsVec[12][0] =  1 * (bound+dPML); pointsVec[12][1] = -1 * bound;
    pointsVec[13][0] =  1 * (bound+dPML); pointsVec[13][1] =  1 * bound;

    pointsVec[14][0] =  1 * bound; pointsVec[14][1] =  1 * (bound+dPML);
    pointsVec[15][0] = -1 * bound; pointsVec[15][1] =  1 * (bound+dPML);

    pointsVec[16][0] = -1 * (bound+dPML); pointsVec[16][1] =  1 * bound;
    pointsVec[17][0] = -1 * (bound+dPML); pointsVec[17][1] = -1 * bound;

    pointsVec[18][0] = -1 * bound; pointsVec[18][1] = -1 * (bound+dPML);
    pointsVec[19][0] = 1 * bound; pointsVec[19][1] = -1 * (bound+dPML);

    pointsVec[20][0] =  1 * (bound+dPML); pointsVec[20][1] = -1 * (bound+dPML);
    pointsVec[21][0] =  1 * (bound+dPML); pointsVec[21][1] =  1 * (bound+dPML);
    pointsVec[22][0] = -1 * (bound+dPML); pointsVec[22][1] =  1 * (bound+dPML);
    pointsVec[23][0] = -1 * (bound+dPML); pointsVec[23][1] = -1 * (bound+dPML);
  }
  TPZFNMatrix<maxNQuads*3,int> quadPointsVec(nQuads,4);
  quadPointsVec(0,0) = 0; quadPointsVec(0,1) = 1; quadPointsVec(0,2) = 2; quadPointsVec(0,3) = 3;
  quadPointsVec(1,0) = 4; quadPointsVec(1,1) = 0; quadPointsVec(1,2) = 3; quadPointsVec(1,3) = 7;
  quadPointsVec(2,0) = 4; quadPointsVec(2,1) = 5; quadPointsVec(2,2) = 1; quadPointsVec(2,3) = 0;
  quadPointsVec(3,0) = 1; quadPointsVec(3,1) = 5; quadPointsVec(3,2) = 6; quadPointsVec(3,3) = 2;
  quadPointsVec(4,0) = 3; quadPointsVec(4,1) = 2; quadPointsVec(4,2) = 6; quadPointsVec(4,3) = 7;

  quadPointsVec(5,0) = 9; quadPointsVec(5,1) = 4; quadPointsVec(5,2) = 7; quadPointsVec(5,3) = 8;
  quadPointsVec(6,0) = 9; quadPointsVec(6,1) = 10; quadPointsVec(6,2) = 5; quadPointsVec(6,3) = 4;
  quadPointsVec(7,0) = 5; quadPointsVec(7,1) = 10; quadPointsVec(7,2) = 11; quadPointsVec(7,3) = 6;
  quadPointsVec(8,0) = 7; quadPointsVec(8,1) = 6; quadPointsVec(8,2) = 11; quadPointsVec(8,3) = 8;
  if(usingPML){
    quadPointsVec(9,0) = 13; quadPointsVec(9,1) = 9; quadPointsVec(9,2) = 8; quadPointsVec(9,3) = 12;
    quadPointsVec(10,0) = 14; quadPointsVec(10,1) = 15; quadPointsVec(10,2) = 10; quadPointsVec(10,3) = 9;
    quadPointsVec(11,0) = 10; quadPointsVec(11,1) = 16; quadPointsVec(11,2) = 17; quadPointsVec(11,3) = 11;
    quadPointsVec(12,0) = 8; quadPointsVec(12,1) = 11; quadPointsVec(12,2) = 18; quadPointsVec(12,3) = 19;

    quadPointsVec(13,0) = 21; quadPointsVec(13,1) = 14; quadPointsVec(13,2) =  9; quadPointsVec(13,3) = 13;
    quadPointsVec(14,0) = 15; quadPointsVec(14,1) = 22; quadPointsVec(14,2) = 16; quadPointsVec(14,3) = 10;
    quadPointsVec(15,0) = 11; quadPointsVec(15,1) = 17; quadPointsVec(15,2) = 23; quadPointsVec(15,3) = 18;
    quadPointsVec(16,0) = 12; quadPointsVec(16,1) =  8; quadPointsVec(16,2) = 19; quadPointsVec(16,3) = 20;
  }
  TPZManVector<int,maxNQuads> matIdsQuads(nQuads,-1);
  matIdsQuads[0] = matIdCore;
  matIdsQuads[1] = matIdCore;    matIdsQuads[2] = matIdCore;
  matIdsQuads[3] = matIdCore;    matIdsQuads[4] = matIdCore;

  matIdsQuads[5] = matIdCladding;    matIdsQuads[6] = matIdCladding;
  matIdsQuads[7] = matIdCladding;    matIdsQuads[8] = matIdCladding;

  if(usingPML){
    matIdsQuads[9] = matIdPMLxp;    matIdsQuads[10] = matIdPMLyp;
    matIdsQuads[11] = matIdPMLxm;    matIdsQuads[12] = matIdPMLym;
    matIdsQuads[13] = matIdPMLxpyp;    matIdsQuads[14] = matIdPMLxmyp;
    matIdsQuads[15] = matIdPMLxmym;    matIdsQuads[16] = matIdPMLxpym;
  }
  
  TPZManVector<int,maxNQuads> nDivQsi(nQuads,-1);
  TPZManVector<int,maxNQuads> nDivEta(nQuads,-1);
  //first hole
  nDivQsi[0]=nPtsCoreT;nDivEta[0]=nPtsCoreT;
  nDivQsi[1]=nPtsCoreR; nDivEta[1]=nPtsCoreT;
  nDivQsi[2]=nPtsCoreT; nDivEta[2]=nPtsCoreR;
  nDivQsi[3]=nPtsCoreR; nDivEta[3]=nPtsCoreT;
  nDivQsi[4]=nPtsCoreT; nDivEta[4]=nPtsCoreR;
  //cladding
  nDivQsi[5]=nPtsCladdingR; nDivEta[5]=nPtsCoreT;
  nDivQsi[6]=nPtsCoreT;     nDivEta[6]=nPtsCladdingR;
  nDivQsi[7]=nPtsCladdingR; nDivEta[7]=nPtsCoreT;
  nDivQsi[8]=nPtsCoreT;     nDivEta[8]=nPtsCladdingR;
  if(usingPML){
    nDivQsi[9]=nDivPml; nDivEta[9]=nPtsCoreT;
    nDivQsi[10]=nPtsCoreT;    nDivEta[10]=nDivPml;
    nDivQsi[11]=nDivPml; nDivEta[11]=nPtsCoreT;
    nDivQsi[12]=nPtsCoreT;     nDivEta[12]=nDivPml;

    nDivQsi[13]=nDivPml; nDivEta[13]=nDivPml;
    nDivQsi[14]=nDivPml;     nDivEta[14]=nDivPml;
    nDivQsi[15]=nDivPml; nDivEta[15]=nDivPml;
    nDivQsi[16]=nDivPml;     nDivEta[16]=nDivPml;
  }
  TPZManVector<bool,maxNQuads> side1NonLinearVec(nQuads,false);
  TPZManVector<bool,maxNQuads> side2NonLinearVec(nQuads,false);
  TPZManVector<bool,maxNQuads> side3NonLinearVec(nQuads,false);
  TPZManVector<bool,maxNQuads> side4NonLinearVec(nQuads,false);

  side4NonLinearVec[1] = true;
  side1NonLinearVec[2] = true;
  side2NonLinearVec[3] = true;
  side3NonLinearVec[4] = true;

  side2NonLinearVec[5] = true;
  side3NonLinearVec[6] = true;
  side4NonLinearVec[7] = true;
  side1NonLinearVec[8] = true;


  TPZManVector<TPZVec<REAL>,14> thetaVec(8,TPZVec<REAL>(2,0.));

  thetaVec[0][0] =-1 * M_PI/4;thetaVec[0][1] = 1 * M_PI/4;//quad1
  thetaVec[1][0] = 1 * M_PI/4;thetaVec[1][1] = 3 * M_PI/4;//quad2
  thetaVec[2][0] = 3 * M_PI/4;thetaVec[2][1] = 5 * M_PI/4;//quad3
  thetaVec[3][0] = 5 * M_PI/4;thetaVec[3][1] = 7 * M_PI/4;//quad4

  thetaVec[4][0] = 1 * M_PI/4;thetaVec[4][1] =-1 * M_PI/4;//quad5
  thetaVec[5][0] = 3 * M_PI/4;thetaVec[5][1] = 1 * M_PI/4;//quad6
  thetaVec[6][0] = 5 * M_PI/4;thetaVec[6][1] = 3 * M_PI/4;//quad7
  thetaVec[7][0] = 7 * M_PI/4;thetaVec[7][1] = 5 * M_PI/4;//quad8


  TPZManVector<TPZVec<REAL>,8> xcRef(8,xc);

  TPZManVector<int,5> matIdBoundVec(4,-1);
  matIdBoundVec[0] = matIdBC;//low
  matIdBoundVec[1] = matIdBC;//down
  matIdBoundVec[2] = matIdBC;//up
  matIdBoundVec[3] = matIdBC;//left

  TPZManVector<REAL,4> boundDistVec(4,-1);
  boundDistVec[0] = 0.;//low
  boundDistVec[1] = usingPML ? bound+dPML : bound;//right
  boundDistVec[2] = usingPML ? bound+dPML : bound;//up
  boundDistVec[3] = 0.;//left

  TPZManVector<REAL,14> rVec(14,rCore);

  auto gmesh = wgma::gmeshtools::CreateStructuredMesh(
    pointsVec, quadPointsVec, matIdsQuads, nDivQsi,nDivEta,
    side1NonLinearVec, side2NonLinearVec, side3NonLinearVec, side4NonLinearVec,
    thetaVec, xcRef, rVec, matIdBoundVec, boundDistVec);

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
    const REAL margin = 0.1 * rCore;
    TPZManVector<REAL,3> qsiPos(2,0);
    TPZManVector<REAL,3> xPos(3,0);
    TPZManVector<TPZGeoEl *,10>sons;
    int nel = gmesh->NElements();
    for(int iel =0 ; iel< nel; iel++){
      auto geo = gmesh->Element(iel);
      if(geo->Type() != EQuadrilateral) continue;
      
      geo->X(qsiPos,xPos);
      const REAL dist = sqrt((xPos[0]-xc[0])*(xPos[0]-xc[0])+(xPos[1]-xc[1])*(xPos[1]-xc[1]));
      if(std::abs(rCore - dist) < margin){
        if(geo->HasSubElement() == 0){
          geo->Divide(sons);
        }
      }
    }
  }
  gmesh->BuildConnectivity();
  //let us split the quadrilaterals into triangles
  wgma::gmeshtools::SplitQuadMeshIntoTriangles(gmesh);

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