#ifndef _CMESHTOOLS_IMPL_
#define _CMESHTOOLS_IMPL_
#include <cmeshtools.hpp>

#include <gmeshtools.hpp>

#include <pzgmesh.h>
#include <pzcmesh.h>
#include <Electromagnetics/TPZCartesianPML.h>
#include <Electromagnetics/TPZCylindricalPML.h>
#include <TPZMatSingleSpace.h>
#include <TPZMatCombinedSpaces.h>
#include <TPZSimpleTimer.h>

template<class MATVOL>
std::map<int,int>
wgma::cmeshtools::AddRectangularPMLRegion(const wgma::pml::cart::data data,
                                          const std::set<int> &volmats,
                                          TPZAutoPointer<TPZGeoMesh> gmesh,
                                          TPZAutoPointer<TPZCompMesh> cmesh)
{

  REAL boundPosX{0}, boundPosY{0}, boundPosZ{0}, dX{0}, dY{0}, dZ{0};




  wgma::gmeshtools::FindPMLWidth(gmesh, data.ids, data.t,
                                 boundPosX, dX,
                                 boundPosY, dY,
                                 boundPosZ, dZ);
  
  

  std::map<int,int> all_neighs;
  for(auto id : data.ids){
    //check if neighbour has been set already, otherwise find it
    const int neigh_mat_id  = [&]{
      if(data.neigh.count(id) == 0){
        const auto pmldim = data.dim;
        const auto neigh_mat_res =
          gmeshtools::FindPMLNeighbourMaterial(gmesh, pmldim, id, volmats);
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
  
    TPZCartesianPML<MATVOL> *pmlMat{nullptr};
    if constexpr (std::is_base_of_v<TPZMatCombinedSpaces, MATVOL>){
      pmlMat = new TPZCombinedSpacesCartesianPML<MATVOL>(id, *neighMat);
    }else{
      pmlMat = new TPZSingleSpaceCartesianPML<MATVOL>(id, *neighMat);
    }

    const bool attx = wgma::pml::cart::attx(data.t);
    const bool atty = wgma::pml::cart::atty(data.t);
    const bool attz = wgma::pml::cart::attz(data.t);
  
    if(attx) pmlMat->SetAttX(boundPosX, data.alphax, dX);
    if(atty) pmlMat->SetAttY(boundPosY, data.alphay, dY);
    if(attz) pmlMat->SetAttZ(boundPosZ, data.alphaz, dZ);

    // std::cout<<"pml";
    // for(auto i: data.ids){std::cout<<'\t'<<i;}
    // std::cout<<"\ntype "<<wgma::pml::to_string(data.t)
    //          <<"\nattx "<<attx<<" dx "<<dX<<" bx "<<boundPosX<<" ax "<<data.alphax
    //          <<"\natty "<<atty<<" dy "<<dY<<" by "<<boundPosY<<" ay "<<data.alphay
    //          <<"\nattz "<<attz<<" dz "<<dZ<<" bz "<<boundPosZ<<" az "<<data.alphaz
    //          <<std::endl;
      
    cmesh->InsertMaterialObject(pmlMat);
  }
  

  
   

  return all_neighs;
}


template<class MATVOL>
std::map<int,int>
wgma::cmeshtools::AddCylindricalPMLRegion(const wgma::pml::cyl::data data,
                                          const std::set<int> &volmats,
                                          TPZAutoPointer<TPZGeoMesh> gmesh,
                                          TPZAutoPointer<TPZCompMesh> cmesh)
{

  REAL rMin{0}, rMax{0}, boundPosZ{0}, dZ{0};


  wgma::gmeshtools::FindPMLWidth(gmesh, data.ids, data.t,
                                 rMin, rMax,
                                 boundPosZ, dZ);
  
  

  std::map<int,int> all_neighs;
  for(auto id : data.ids){
    //check if neighbour has been set already, otherwise find it
    const int neigh_mat_id  = [&]{
      if(data.neigh.count(id) == 0){
        const auto pmldim = data.dim;
        const auto neigh_mat_res =
          gmeshtools::FindPMLNeighbourMaterial(gmesh, pmldim, id, volmats);
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
  
    TPZCylindricalPML<MATVOL> *pmlMat{nullptr};
    if constexpr (std::is_base_of_v<TPZMatCombinedSpaces, MATVOL>){
      pmlMat = new TPZCombinedSpacesCylindricalPML<MATVOL>(id, *neighMat);
    }else{
      pmlMat = new TPZSingleSpaceCylindricalPML<MATVOL>(id, *neighMat);
    }

    const bool attr = wgma::pml::cyl::attr(data.t);
    const bool attz = wgma::pml::cyl::attz(data.t);
  
    if(attr) pmlMat->SetAttR(rMin, data.alphar, rMax-rMin);
    if(attz) pmlMat->SetAttZ(boundPosZ, data.alphaz, dZ);
    cmesh->InsertMaterialObject(pmlMat);
  }
  return all_neighs;
}
#endif