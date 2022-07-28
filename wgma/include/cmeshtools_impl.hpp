#ifndef _CMESHTOOLS_IMPL_
#define _CMESHTOOLS_IMPL_
#include <cmeshtools.hpp>

#include <gmeshtools.hpp>

#include <pzgmesh.h>
#include <pzcmesh.h>
#include <Electromagnetics/TPZMatPML.h>
#include <TPZMatSingleSpace.h>
#include <TPZMatCombinedSpaces.h>
template<class MATVOL>
std::map<int,int>
wgma::cmeshtools::AddRectangularPMLRegion(const wgma::pml::data data,
                                          const std::set<int> &volmats,
                                          TPZAutoPointer<TPZGeoMesh> gmesh,
                                          TPZAutoPointer<TPZCompMesh> cmesh)
{

  REAL boundPosX{0}, boundPosY{0}, dX{0}, dY{0};


  wgma::gmeshtools::FindPMLWidth(gmesh, data.ids, data.t,
                                 boundPosX, dX, boundPosY, dY);
  
  

  std::map<int,int> all_neighs;
  for(auto id : data.ids){
    //check if neighbour has been set already, otherwise find it
    const int neigh_mat_id  = [&]{
      if(data.neigh.count(id) == 0){
        const auto neigh_mat_res =
          gmeshtools::FindPMLNeighbourMaterial(gmesh, id, volmats, boundPosX, boundPosY);
        if(neigh_mat_res.has_value() == false){
          PZError<<__PRETTY_FUNCTION__
                 <<"Could not find neighbouring material. Aborting...\n";
          DebugStop();
        }
        return neigh_mat_res.value();
      }else{
        return data.neigh.at(id);
      }
    }();
    
    all_neighs[id] = neigh_mat_id;
    
    auto *neighMat =
      dynamic_cast<MATVOL*>(cmesh->FindMaterial(neigh_mat_id));
    if(!neighMat){
      PZError<<__PRETTY_FUNCTION__;
      PZError<<"\n neighbouring material not found in mesh, aborting...."<<std::endl;
      DebugStop();
    }
  
    TPZMatPML<MATVOL> *pmlMat{nullptr};
    if constexpr (std::is_base_of_v<TPZMatCombinedSpaces, MATVOL>){
      pmlMat = new TPZCombinedSpacesPML<MATVOL>(id, *neighMat);
    }else{
      pmlMat = new TPZSingleSpacePML<MATVOL>(id, *neighMat);
    }

    const bool attx = wgma::pml::attx(data.t);
    const bool atty = wgma::pml::atty(data.t);
  
    if(attx) pmlMat->SetAttX(boundPosX, data.alphax, dX);
    if(atty) pmlMat->SetAttY(boundPosY, data.alphay, dY);
    cmesh->InsertMaterialObject(pmlMat);
  }
  

  
   

  return all_neighs;
}
#endif