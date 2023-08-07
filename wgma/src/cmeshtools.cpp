#include "cmeshtools.hpp"
#include "pmltypes.hpp"
#include "cmeshtools_impl.hpp"
#include "gmeshtools.hpp"

#include <pzgmesh.h>
#include <pzcmesh.h>
#include <TPZBndCond.h>
#include <pzvec_extras.h>
#include <TPZSimpleTimer.h>



#include <regex>//for string search

using namespace wgma;


void
cmeshtools::SetupGmshMaterialData(
  const TPZVec<std::map<std::string,int>> &gmshmats,
  const std::map<std::string,std::pair<CSTATE,CSTATE>> &matmap,
  const std::map<std::string,wgma::bc::type> &bcmap,
  TPZVec<STATE> alphaPML,
  wgma::cmeshtools::PhysicalData &data,
  int dim)
{

  auto &matinfo = data.matinfovec;
  auto &pmlvec = data.pmlvec;
  auto &bcvec = data.bcvec;
  pmlvec.resize(0);

  if(alphaPML.size() < dim){
    PZError<<__PRETTY_FUNCTION__
           <<"not enough PML coefficients!\n"
           <<"dim = "<<dim<<" # pml coeff "<<alphaPML.size()
           <<"\nAborting...";
    DebugStop();
  }

  
  alphaPML.Resize(3,0);
  const auto &alphaPMLx = alphaPML[0];
  const auto &alphaPMLy = alphaPML[1];
  const auto &alphaPMLz = alphaPML[2];

  if(dim == -1){
    //find maximum dim
    for(int idim = 3; idim >=0; idim--){
      if(!gmshmats[idim].empty()){
        dim = idim;
        break;
      }
    }
  }


  //lambda for easier creation of PML mats
  auto SetupIfPml = [&pmlvec, &alphaPMLx, &alphaPMLy, &alphaPMLz](const std::string &name,
                                                                  int id,
                                                                  int dim) -> bool{

    constexpr auto pmlname{"pml"};
    //std::regex_constants::icase - ignores case
    const auto rx = std::regex{ pmlname ,std::regex_constants::icase };
    const bool ispml = std::regex_search(name, rx);
    if(ispml){
      auto cart_pml = wgma::pml::cart::IdentifyAndSetupPML(name, id, dim,
                                                           alphaPMLx,alphaPMLy,
                                                           alphaPMLz);
      if(cart_pml){
        pmlvec.push_back(cart_pml);
        return true;
      }

      auto cyl_pml = wgma::pml::cyl::IdentifyAndSetupPML(name, id, dim,
                                                         alphaPMLx,
                                                         alphaPMLz);
      if(cyl_pml){
        pmlvec.push_back(cyl_pml);
        return true;
      }

      std::cout<<"material "<<name
               <<" with id "<<id
               <<" could not be identified as a PML material!"<<std::endl;
      DebugStop();
    
      return true;
    }
    return false;
  };
  
  for(auto mat : gmshmats[dim]){
    const std::string name = mat.first;
    const auto id = mat.second;

    if(auto ispml = SetupIfPml(name,id,dim); !ispml){
      if(matmap.find(name) == matmap.end()){//material not found
        std::cout<<"info: mat "<<name<<" id "<<id<<" not found"
                 <<"\nSkipping..."<<std::endl;
      }else{
        const auto pos = matinfo.size();
        matinfo.push_back(std::make_tuple(
                            id,
                            matmap.at(name).first,
                            matmap.at(name).second)
                          );
      }
    }
  }
  //bc materials
  {
    for(auto bc : gmshmats[dim-1]){
      const std::string name = bc.first;
      const auto id = bc.second;

      if(auto ispml = SetupIfPml(name,id,dim-1); !ispml){
        /*sometimes we need 1d materials with a given id for other purposes,
          such as circumference arcs. soo, not finding it is not a problem*/
        if(bcmap.find(name) != bcmap.end()){
          const int ibc = bcvec.size();
          bcvec.resize(ibc+1);
          bcvec[ibc].id = id;
          bcvec[ibc].t = bcmap.at(name);
          bcvec[ibc].name = name;
        }
      }
    }
  }
    
}

void
cmeshtools::FindDirichletConnects(TPZAutoPointer<TPZCompMesh> cmesh,
                                  std::set<int64_t> &boundConnects,
                                  const std::set<int> & matIds)
{
  boundConnects.clear();
  std::set<int> all_mats = matIds;
  //if set is empty, we insert all materials in it
  if(all_mats.size() == 0){
    for(auto [id,ptr] : cmesh->MaterialVec()){
      all_mats.insert(id);
    }
  }

  for (int iel = 0; iel < cmesh->NElements(); iel++) {
    TPZCompEl *cel = cmesh->ElementVec()[iel];
    if (cel == nullptr) {
      continue;
    }
    if (cel->Reference() == nullptr) {//there is no associated geometric el
      continue;
    }
    TPZBndCond *mat = dynamic_cast<TPZBndCond *>(
        cmesh->MaterialVec()[cel->Reference()->MaterialId()]);
    if (mat && mat->Type() == 0 && all_mats.count(mat->Id())) {//check for dirichlet bcs
      std::set<int64_t> boundConnectsEl;
      cel->BuildConnectList(boundConnectsEl);
      for(auto val : boundConnectsEl){
        if (boundConnects.find(val) == boundConnects.end()) {
          boundConnects.insert(val);
        }
      }
    }
  }
}

void
cmeshtools::FilterBoundaryEquations(TPZAutoPointer<TPZCompMesh> cmesh,
                                    TPZVec<int64_t> &activeEquations,
                                    std::set<int64_t> &boundConnects)
{
  TPZSimpleTimer timer ("Filter dirichlet eqs");  

  FindDirichletConnects(cmesh, boundConnects);

  //certainly we have less equations than this, but we will avoid repeated resizes
  activeEquations.Resize(cmesh->NEquations());
  int neq = 0;
  for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
    if (boundConnects.find(iCon) == boundConnects.end()) {
      TPZConnect &con = cmesh->ConnectVec()[iCon];
      const auto condensed = con.IsCondensed();
      const auto hasdep = con.HasDependency();
      const auto seqnum = con.SequenceNumber();
      const auto pos = cmesh->Block().Position(seqnum);
      const auto blocksize = cmesh->Block().Size(seqnum);
      
      if(condensed || hasdep || seqnum < 0 || !blocksize) { continue; }
      const auto vs = neq;
      for (auto ieq = 0; ieq < blocksize; ieq++) {
        activeEquations[vs + ieq] = pos + ieq;
      }
      neq += blocksize;
    }
  }
  activeEquations.Resize(neq);
}

void cmeshtools::RemovePeriodicity(TPZAutoPointer<TPZCompMesh> cmesh)
{
  for(auto &con : cmesh->ConnectVec()){
    con.RemoveDepend();
  }
}