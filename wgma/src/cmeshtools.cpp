#include "cmeshtools.hpp"
#include "pmltypes.hpp"
#include "cmeshtools_impl.hpp"
#include "gmeshtools.hpp"

#include <pzgmesh.h>
#include <pzcmesh.h>
#include <TPZBndCond.h>
#include <TPZSimpleTimer.h>


#include <regex>//for string search

using namespace wgma;


void
cmeshtools::SetupGmshMaterialData(
  const TPZVec<std::map<std::string,int>> &gmshmats,
  const std::map<std::string,std::pair<CSTATE,CSTATE>> &matmap,
  const std::map<std::string,wgma::bc::type> &bcmap,
  const STATE alphaPMLx,
  const STATE alphaPMLy,
  wgma::cmeshtools::PhysicalData &data)
{
  auto &matinfo = data.matinfovec;
  auto &pmlvec = data.pmlvec;
  auto &bcvec = data.bcvec;
  pmlvec.resize(0);
    
  for(auto mat : gmshmats[2]){
    const std::string name = mat.first;
    const auto id = mat.second;
    constexpr auto pmlname{"pml"};

    //std::regex_constants::icase - ignores case
    const auto rx = std::regex{ pmlname ,std::regex_constants::icase };
    const bool ispml = std::regex_search(name, rx);
      
    if(ispml){
      
      const auto pos = pmlvec.size();
      bool found{false};
      const TPZVec<std::string> pmlnames =
        {"xpyp", "xmyp", "xmym", "xpym", "xp", "yp", "xm", "ym"};
      const TPZVec<wgma::pml::type> pmltypes =
        {wgma::pml::type::xpyp,wgma::pml::type::xmyp,
         wgma::pml::type::xmym,wgma::pml::type::xpym,
         wgma::pml::type::xp,wgma::pml::type::yp,
         wgma::pml::type::xm,wgma::pml::type::ym};
      
      for(int ipml = 0; ipml < pmlnames.size(); ipml++){
        const std::string pmlname = pmlnames[ipml];
        const auto rx = std::regex{ pmlname, std::regex_constants::icase };
        const bool test = std::regex_search(name, rx);
        if(test){
          pmlvec.resize(pos+1);
          pmlvec[pos].id = id;
          pmlvec[pos].alphax = alphaPMLx;
          pmlvec[pos].alphay = alphaPMLy;
          pmlvec[pos].t = pmltypes[ipml];
          pmlvec[pos].name = name;
          found = true;
          break;
        }
      }
      //pml was not identified
      if(!found){
        std::cout<<"error: mat "<<name<<" id "<<id<<" not found"
                 <<"\nAborting..."<<std::endl;
        DebugStop();
      }
    }else{
      if(matmap.find(name) == matmap.end()){//material not found
        std::cout<<"error: mat "<<name<<" id "<<id<<" not found"
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
    for(auto bc : gmshmats[1]){
      const auto name = bc.first;
      const auto id = bc.second;
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

void
cmeshtools::FilterBoundaryEquations(TPZAutoPointer<TPZCompMesh> cmesh,
                                    TPZVec<int64_t> &activeEquations,
                                    std::set<int64_t> &boundConnects)
{
  TPZSimpleTimer timer ("Filter dirichlet eqs");
  TPZManVector<int64_t, 1000> allConnects;
  boundConnects.clear();

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
    if (mat && mat->Type() == 0) {//check for dirichlet bcs
      std::set<int64_t> boundConnectsEl;
      std::set<int64_t> depBoundConnectsEl;
      std::set<int64_t> indepBoundConnectsEl;
      cel->BuildConnectList(indepBoundConnectsEl, depBoundConnectsEl);
      cel->BuildConnectList(boundConnectsEl);
      for(auto val : boundConnectsEl){
        if (boundConnects.find(val) == boundConnects.end()) {
          boundConnects.insert(val);
        }
      }
    }
  }
  for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
    if (boundConnects.find(iCon) == boundConnects.end()) {
      TPZConnect &con = cmesh->ConnectVec()[iCon];
      if (con.HasDependency())
        continue;
      const auto seqnum = con.SequenceNumber();
      if(seqnum < 0) { continue; }
      const auto pos = cmesh->Block().Position(seqnum);
      const auto blocksize = cmesh->Block().Size(seqnum);
      if (blocksize == 0){ continue; }

      const auto vs = activeEquations.size();
      activeEquations.Resize(vs + blocksize);
      for (auto ieq = 0; ieq < blocksize; ieq++) {
        activeEquations[vs + ieq] = pos + ieq;
      }
    }
  }
}

void cmeshtools::RemovePeriodicity(TPZAutoPointer<TPZCompMesh> cmesh)
{
  for(auto &con : cmesh->ConnectVec()){
    con.RemoveDepend();
  }
}