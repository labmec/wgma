#include <cmeshtools.hpp>

#include <pzgmesh.h>
#include <pzcmesh.h>

template<class MATPML, class MATVOL>
void
wgma::cmeshtools::AddRectangularPMLRegion(const int matId,
                                    const int alphax, const int alphay,
                                    const wgma::pml::type type,
                                    const std::set<int> &volmats,
                                    TPZAutoPointer<TPZGeoMesh> gmesh,
                                    TPZAutoPointer<TPZCompMesh> cmesh)
{

  //let us find the (xmin,xmax) and (ymin,ymax) of the PML region
  REAL xMax = -1e20, xMin = 1e20, yMax = -1e20, yMin = 1e20;
  for (auto geo : gmesh->ElementVec()){
    if (geo && geo->MaterialId() == matId) {
      for (int iNode = 0; iNode < geo->NCornerNodes(); ++iNode) {
        TPZManVector<REAL, 3> co(3);
        geo->Node(iNode).GetCoordinates(co);
        const REAL &xP = co[0];
        const REAL &yP = co[1];
        if (xP > xMax) {
          xMax = xP;
        }
        if (xP < xMin) {
          xMin = xP;
        }
        if (yP > yMax) {
          yMax = yP;
        }
        if (yP < yMin) {
          yMin = yP;
        }
      }
    }
  }


  //now we compute xBegin, yBegin, attx, atty and d for the material ctor
  const bool attx = wgma::pml::attx(type);
  const bool atty = wgma::pml::atty(type);
  
  REAL xBegin{-1}, yBegin{-1}, dX{-01101991.}, dY{-01101991.};
  REAL boundPosX{-01101991.}, boundPosY{-01101991.};
  if(attx){
    const int xdir = wgma::pml::xinfo(type);
    dX = xMax - xMin;
    xBegin = xdir > 0 ? xMin : xMax;
    boundPosX = xBegin;
    boundPosY = (yMax + yMin)/2;
  }

  if(atty){
    const int ydir = wgma::pml::yinfo(type);
    dY = yMax - yMin;
    yBegin = ydir > 0 ? yMin : yMax;
    boundPosX = (xMax + xMin)/2;
    boundPosY = yBegin;
  }

  if(attx && atty){
    boundPosX = xBegin;
    boundPosY = yBegin;
  }
  //find the neighbouring material
  const auto neighMatId =
    FindPMLNeighbourMaterial(gmesh, matId, volmats, boundPosX, boundPosY);

  auto neighMat = dynamic_cast<MATVOL*>(
    cmesh->FindMaterial(neighMatId));
  if(!neighMat){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\n neighbouring material not found in mesh, aborting...."<<std::endl;
    DebugStop();
  }
  
  auto pmlMat = new MATPML(matId, *neighMat);
  if(attx) pmlMat->SetAttX(xBegin, alphax, dX);
  if(atty) pmlMat->SetAttY(yBegin, alphay, dY);
  cmesh->InsertMaterialObject(pmlMat);
  
}