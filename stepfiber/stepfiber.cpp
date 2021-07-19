/**
stepfiber.cpp

This target performs the modal analysis of a step-index optical fiber.

A HCurl-conforming approximation space is used for the transverse field components and
a H1-conforming approximation space for the axial component.
***/

//wgma includes
#include "cmeshtools.hpp"
#include "gmeshtools.hpp"
#include "pmltypes.hpp"
//pz includes
#include <MMeshType.h>                   //for MMeshType
#include <pzcmesh.h>                     //for TPZCompMesh
#include <pzgmesh.h>                     //for TPZGeoMesh
#include <tpzgeoelrefpattern.h>          //for TPZGeoElRefPattern
#include <TPZGeoMeshTools.h>             //for TPZGeoMeshTools::CreateGeoMeshOnGrid
#include <TPZRefPatternDataBase.h>       //for TPZRefPatternDataBase
#include <pzlog.h>                       //for TPZLogger
#include <TPZElectromagneticConstants.h> //for pzelectromag::cZero
#include <TPZEigenAnalysis.h>            //for TPZEigenAnalysis 
#include <TPZKrylovEigenSolver.h>        //for TPZKrylovEigenSolver
#include <TPZVTKGeoMesh.h>               //for exporting geomesh to vtk
#include <TPZSkylineNSymStructMatrix.h>  //non-symmetric skyline matrix storage
#include <TPZSpStructMatrix.h>           //non-symmetric sparse matrix storage
#include <pzbuildmultiphysicsmesh.h>     //for TPZBuildMultiphysicsMesh
#include <Electromagnetics/TPZWaveguideModalAnalysis.h> //for TPZMatWaveguideModalAnalysis
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
   @param[in] print whether to print the gmesh (.vtk and .txt formats)
   @param[in] prefix prefix used in the name of the gmesh files
   @param[out] matIdVec vector with the material identifiers: core/cladding/pmls/bcs
   @param[out] pmlTypeVec vector with the pml types (in the same order as `matIdVec`)
*/
TPZAutoPointer<TPZGeoMesh>
CreateStepFiberMesh(
  const REAL rCore, const REAL boundDist, const REAL dPML, const int nLayersPML,
  const int factor, bool refine, const REAL scale, const bool print,
  const std::string &prefix, TPZVec<int> &matIdVec,
  TPZVec<wgma::pmltype> &pmlTypeVec
  );


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
  constexpr REAL rCore{0.000008};
  //distance from the core from which the PML begins
  constexpr REAL boundDist{4*claddingReffIndex*2*M_PI};
  //pml width
  constexpr REAL dPML{7.5*claddingReffIndex*2*M_PI};
  //number of layers in the pml
  constexpr int nLayersPML{2};
  // operational wavelength
  constexpr STATE lambda{1.55e-6};
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
   the effective index neff).*/
  constexpr CSTATE target = -2.086289757;
  // Dimension of the krylov space to be used. Suggested to be at least nev * 10
  constexpr int krylovDim{500};


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
  TPZManVector<wgma::pmltype,8> pmlTypeVec;
  constexpr int factor{3};
  constexpr bool refine{true};
  //whether to print the geometric mesh in .vtk and .txt files
  constexpr bool print{true};
  //prefix for the gmesh files
  const std::string prefix = refine ? "ref" : "noref";
  // const REAL realRCore = simData.physicalOpts.stepFiberOpts.realRCore;
  // const REAL dPML = simData.physicalOpts.stepFiberOpts.dPML;
  // const REAL boundDist = simData.physicalOpts.stepFiberOpts.boundDist;
  // const REAL outerMaterialReffIndex = std::sqrt(std::real(simData.physicalOpts.erVec[1]));
  // const int nLayersPML = simData.physicalOpts.stepFiberOpts.nLayersPml;
  auto gmesh = CreateStepFiberMesh(rCore,boundDist,dPML,nLayersPML,
                                   factor,refine,scale,print,prefix,
                                   matIdVec,pmlTypeVec);
  return 0;
}


TPZAutoPointer<TPZGeoMesh>
CreateStepFiberMesh(
  const REAL realRCore, const REAL boundDist, const REAL dPML, const int nLayersPml,
  const int factor, bool refine, const REAL scale, const bool print,
  const std::string &prefix, TPZVec<int> &matIdVec,
  TPZVec<wgma::pmltype> &pmlTypeVec)
{

  using wgma::gmeshtools::EdgeData;
  using wgma::gmeshtools::QuadrilateralData;
  const int nPtsCoreR = factor * 2,//number of points in core (radial direction)
    nPtsCoreT = factor * 5,//number of points in core (tangential direction)
    nPtsCladdingR = factor * 4,
    nDivPml = factor * nLayersPml + 1;

  if(std::min<int>({nPtsCoreR,nPtsCoreT,nPtsCladdingR,nDivPml}) < 2 ) {
    std::cout<<"Mesh has not sufficient divisions."<<std::endl;
    std::cout<<"Minimum is 2."<<std::endl;
    DebugStop();
  }

  constexpr int nQuads = 17;
  constexpr int nPoints = 24;

  constexpr int matIdCore = 1, matIdCladding = 2, matIdBC= 18;
  constexpr int matIdPMLxp = 10,
    matIdPMLyp = 11,
    matIdPMLxm = 12,
    matIdPMLym = 13,
    matIdPMLxpym = 14,
    matIdPMLxpyp = 15,
    matIdPMLxmyp = 16,
    matIdPMLxmym = 17;

  const REAL rCore = scale * realRCore;
  const REAL bound = rCore + boundDist;
  TPZManVector<REAL,2> xc(2, 0.);
  xc[0] = 0.;
  xc[1] = 0.;



  ///all points
  TPZManVector<TPZVec<REAL>,33> pointsVec(nPoints,TPZVec<REAL>(2,0.));
  TPZManVector<EdgeData,33> edgesVec(0,EdgeData());
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
  //PML
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

  TPZFMatrix<int> quadPointsVec(nQuads,4);
  quadPointsVec(0,0) = 0; quadPointsVec(0,1) = 1; quadPointsVec(0,2) = 2; quadPointsVec(0,3) = 3;
  quadPointsVec(1,0) = 4; quadPointsVec(1,1) = 0; quadPointsVec(1,2) = 3; quadPointsVec(1,3) = 7;
  quadPointsVec(2,0) = 4; quadPointsVec(2,1) = 5; quadPointsVec(2,2) = 1; quadPointsVec(2,3) = 0;
  quadPointsVec(3,0) = 1; quadPointsVec(3,1) = 5; quadPointsVec(3,2) = 6; quadPointsVec(3,3) = 2;
  quadPointsVec(4,0) = 3; quadPointsVec(4,1) = 2; quadPointsVec(4,2) = 6; quadPointsVec(4,3) = 7;

  quadPointsVec(5,0) = 9; quadPointsVec(5,1) = 4; quadPointsVec(5,2) = 7; quadPointsVec(5,3) = 8;
  quadPointsVec(6,0) = 9; quadPointsVec(6,1) = 10; quadPointsVec(6,2) = 5; quadPointsVec(6,3) = 4;
  quadPointsVec(7,0) = 5; quadPointsVec(7,1) = 10; quadPointsVec(7,2) = 11; quadPointsVec(7,3) = 6;
  quadPointsVec(8,0) = 7; quadPointsVec(8,1) = 6; quadPointsVec(8,2) = 11; quadPointsVec(8,3) = 8;

  quadPointsVec(9,0) = 13; quadPointsVec(9,1) = 9; quadPointsVec(9,2) = 8; quadPointsVec(9,3) = 12;
  quadPointsVec(10,0) = 14; quadPointsVec(10,1) = 15; quadPointsVec(10,2) = 10; quadPointsVec(10,3) = 9;
  quadPointsVec(11,0) = 10; quadPointsVec(11,1) = 16; quadPointsVec(11,2) = 17; quadPointsVec(11,3) = 11;
  quadPointsVec(12,0) = 8; quadPointsVec(12,1) = 11; quadPointsVec(12,2) = 18; quadPointsVec(12,3) = 19;

  quadPointsVec(13,0) = 21; quadPointsVec(13,1) = 14; quadPointsVec(13,2) =  9; quadPointsVec(13,3) = 13;
  quadPointsVec(14,0) = 15; quadPointsVec(14,1) = 22; quadPointsVec(14,2) = 16; quadPointsVec(14,3) = 10;
  quadPointsVec(15,0) = 11; quadPointsVec(15,1) = 17; quadPointsVec(15,2) = 23; quadPointsVec(15,3) = 18;
  quadPointsVec(16,0) = 12; quadPointsVec(16,1) =  8; quadPointsVec(16,2) = 19; quadPointsVec(16,3) = 20;

  TPZManVector<int,nQuads> matIdsQuads(nQuads,-1);
  matIdsQuads[0] = matIdCore;
  matIdsQuads[1] = matIdCore;    matIdsQuads[2] = matIdCore;
  matIdsQuads[3] = matIdCore;    matIdsQuads[4] = matIdCore;

  matIdsQuads[5] = matIdCladding;    matIdsQuads[6] = matIdCladding;
  matIdsQuads[7] = matIdCladding;    matIdsQuads[8] = matIdCladding;
  matIdsQuads[9] = matIdPMLxp;    matIdsQuads[10] = matIdPMLyp;
  matIdsQuads[11] = matIdPMLxm;    matIdsQuads[12] = matIdPMLym;
  matIdsQuads[13] = matIdPMLxpyp;    matIdsQuads[14] = matIdPMLxmyp;
  matIdsQuads[15] = matIdPMLxmym;    matIdsQuads[16] = matIdPMLxpym;

  TPZManVector<int,nQuads> nDivQsi(nQuads,-1);
  TPZManVector<int,nQuads> nDivEta(nQuads,-1);
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

  nDivQsi[9]=nDivPml; nDivEta[9]=nPtsCoreT;
  nDivQsi[10]=nPtsCoreT;    nDivEta[10]=nDivPml;
  nDivQsi[11]=nDivPml; nDivEta[11]=nPtsCoreT;
  nDivQsi[12]=nPtsCoreT;     nDivEta[12]=nDivPml;

  nDivQsi[13]=nDivPml; nDivEta[13]=nDivPml;
  nDivQsi[14]=nDivPml;     nDivEta[14]=nDivPml;
  nDivQsi[15]=nDivPml; nDivEta[15]=nDivPml;
  nDivQsi[16]=nDivPml;     nDivEta[16]=nDivPml;

  TPZManVector<bool,nQuads> side1NonLinearVec(nQuads,false);
  TPZManVector<bool,nQuads> side2NonLinearVec(nQuads,false);
  TPZManVector<bool,nQuads> side3NonLinearVec(nQuads,false);
  TPZManVector<bool,nQuads> side4NonLinearVec(nQuads,false);

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


  TPZManVector<TPZVec<REAL>*,8> xcRef(8,&xc);

  TPZManVector<int,5> matIdBoundVec(4,-1);
  matIdBoundVec[0] = matIdBC;//low
  matIdBoundVec[1] = matIdBC;//down
  matIdBoundVec[2] = matIdBC;//up
  matIdBoundVec[3] = matIdBC;//left

  TPZManVector<REAL,4> boundDistVec(4,-1);
  boundDistVec[0] = 0.;//low
  boundDistVec[1] = bound+dPML;//down
  boundDistVec[2] = bound+dPML;//up
  boundDistVec[3] = 0.;//left

  TPZManVector<REAL,14> rVec(14,rCore);

  auto gmesh = wgma::gmeshtools::CreateStructuredMesh(
    pointsVec, edgesVec, quadPointsVec, matIdsQuads, nDivQsi,nDivEta,
    side1NonLinearVec, side2NonLinearVec, side3NonLinearVec, side4NonLinearVec,
    thetaVec, xcRef, matIdBoundVec, boundDistVec, rVec);
  matIdVec.Resize(11);
  matIdVec[0]=matIdCore;
  matIdVec[1]=matIdCladding;
  pmlTypeVec.Resize(8);
  pmlTypeVec[0]=wgma::pmltype::xp;
  matIdVec[2]=matIdPMLxp;
  pmlTypeVec[1]=wgma::pmltype::yp;
  matIdVec[3]=matIdPMLyp;
  pmlTypeVec[2]=wgma::pmltype::xm;
  matIdVec[4]=matIdPMLxm;
  pmlTypeVec[3]=wgma::pmltype::ym;
  matIdVec[5]=matIdPMLym;
  pmlTypeVec[4]=wgma::pmltype::xpym;
  matIdVec[6]=matIdPMLxpym;
  pmlTypeVec[5]=wgma::pmltype::xpyp;
  matIdVec[7]=matIdPMLxpyp;
  pmlTypeVec[6]=wgma::pmltype::xmyp;
  matIdVec[8]=matIdPMLxmyp;
  pmlTypeVec[7]=wgma::pmltype::xmym;
  matIdVec[9]=matIdPMLxmym;
  
  matIdVec[10]=matIdBC;

  gmesh->BuildConnectivity();

  //refpatquadtri.rpt
  long nel = gmesh->NElements();

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

  if(print){
    std::string meshFileName = prefix + "gmesh";
    const size_t strlen = meshFileName.length();
    meshFileName.append(".vtk");
    std::ofstream outVTK(meshFileName.c_str());
    meshFileName.replace(strlen, 4, ".txt");
    std::ofstream outTXT(meshFileName.c_str());

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
    gmesh->Print(outTXT);
    outTXT.close();
    outVTK.close();
  }

  return gmesh;
}
