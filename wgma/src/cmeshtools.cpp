#include "cmeshtools.hpp"
#include "pmltypes.hpp"
#include "cmeshtools_impl.hpp"
#include "gmeshtools.hpp"

#include <pzgmesh.h>
#include <pzcmesh.h>
#include <TPZBndCond.h>
#include <pzvec_extras.h>
#include <TPZSimpleTimer.h>
#include <TPZParallelUtils.h>


#include <regex>//for string search

using namespace wgma;

void wgma::cmeshtools::PrintCompMesh(TPZAutoPointer<TPZCompMesh> cmesh,
                                     std::string filename)
{
  const std::string cmeshFileNameTxt{filename+".txt"};
  std::ofstream cmeshFileTxt(cmeshFileNameTxt);
  cmesh->Print(cmeshFileTxt);
}


void
cmeshtools::SetupGmshMaterialData(
  const TPZVec<std::map<std::string,int>> &gmshmats,
  const std::map<std::string,std::pair<CSTATE,CSTATE>> &matmap,
  const std::map<std::string,wgma::bc::type> &bcmap,
  TPZVec<CSTATE> alphaPML,
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

void cmeshtools::SetPeriodic(TPZAutoPointer<TPZCompMesh> &cmesh,
                             const std::map<int64_t,int64_t> &periodic_els)
{
  auto gmesh = cmesh->Reference();
  gmesh->ResetReference();
  cmesh->LoadReferences();
  //let us copy the connects
  for(auto [dep, indep] : periodic_els){
    //geometric elements
    auto *dep_gel = gmesh->Element(dep);
    const auto *indep_gel = gmesh->Element(indep);
    //computational element
    auto *indep_cel = indep_gel->Reference();
    auto *dep_cel = dep_gel->Reference();
    //not necessarily all elements belong to the comp mesh
    if(!indep_cel && !dep_cel){continue;}
    if(!indep_cel || !dep_cel){
      PZError<<__PRETTY_FUNCTION__
             <<"\nError: found only one of periodic els!\n";
      if(!indep_cel){
        PZError<<"Could not find comp el associated with geo el "<<indep<<std::endl;
      }
      if(!dep_cel){
        PZError<<"Could not find comp el associated with geo el "<<dep<<std::endl;
      }
    }
    //number of connects
    const auto n_dep_con = dep_cel->NConnects();
    const auto n_indep_con = indep_cel->NConnects();
    //just to be sure
    if(n_dep_con != n_indep_con){
      PZError<<__PRETTY_FUNCTION__
             <<"\nindep cel "<<indep_cel->Index()
             <<" has "<<n_indep_con<<" connects"
             <<"\ndep cel "<<dep_cel->Index()
             <<" has "<<n_dep_con<<" connects"<<std::endl;
      DebugStop();
    }

    //now we create dependencies between connects
    for(auto ic = 0; ic < n_indep_con; ic++){
      const auto indep_ci = indep_cel->ConnectIndex(ic);
      const auto dep_ci = dep_cel->ConnectIndex(ic);

      auto &dep_con = dep_cel->Connect(ic);
      const auto ndof = dep_con.NDof(cmesh);
      if(ndof==0) {continue;}
      constexpr int64_t ipos{0};
      constexpr int64_t jpos{0};
      
      TPZFMatrix<REAL> mat(ndof,ndof);
      mat.Identity();
      dep_con.AddDependency(dep_ci, indep_ci, mat, ipos,jpos,ndof,ndof);
    } 
  }
  cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();
}

void cmeshtools::RemovePeriodicity(TPZAutoPointer<TPZCompMesh> cmesh)
{
  for(auto &con : cmesh->ConnectVec()){
    con.RemoveDepend();
  }
}

void cmeshtools::ExtractSolFromMesh(TPZAutoPointer<TPZCompMesh> mesh_dest,
                                    TPZAutoPointer<TPZCompMesh> mesh_orig,
                                    TPZFMatrix<CSTATE> &sol_dest)
{
  sol_dest.Zero();
  mesh_orig->LoadReferences();
  const TPZFMatrix<CSTATE> &sol_orig = mesh_orig->Solution();
  const auto &block_dest = mesh_dest->Block();
  const auto &block_orig = mesh_orig->Block();
  const int nel = mesh_dest->NElements();
  pzutils::ParallelFor(0,nel,[mesh_dest,&block_dest,&block_orig,
                              &sol_dest,&sol_orig](int iel){
    auto el_dest = mesh_dest->Element(iel);
    if(el_dest->Dimension() < mesh_dest->Dimension()){return;}
    //geometric mesh has orig mesh as reference
    auto el_orig = el_dest->Reference()->Reference();
    if(!el_orig){DebugStop();}
    const auto ncon = el_orig->NConnects();
    if(ncon != el_dest->NConnects()){
      DebugStop();
    }
    for(auto icon = 0; icon < ncon; icon++){
      auto &con_orig = el_orig->Connect(icon);
      auto &con_dest = el_dest->Connect(icon);
        
      const int64_t dfseq_dest = con_dest.SequenceNumber();
      const int64_t pos_dest = block_dest.Position(dfseq_dest);

      const int64_t dfseq_orig = con_orig.SequenceNumber();
      const int64_t pos_orig = block_orig.Position(dfseq_orig);
        
      const int nshape = block_dest.Size(dfseq_dest);
      if(nshape!= block_orig.Size(dfseq_orig)){
        DebugStop();
      }
      for(int is = 0; is < nshape; is++){
        sol_dest.Put(pos_dest+is,0,sol_orig.Get(pos_orig+is,0));
      }
        
    }
  });
}

int64_t
cmeshtools::RestrictDofs(TPZAutoPointer<TPZCompMesh> scatt_mesh,
                         TPZAutoPointer<TPZCompMesh> modal_mesh,
                         const int nm,
                         const std::set<int64_t> &dirichlet_connects)
{
  //so we can access the dofs
  TPZFMatrix<CSTATE> &modal_sol = modal_mesh->Solution();

  const int nmodes_total = modal_sol.Cols();
  if(nm > nmodes_total){
    std::cout<<"nm "<<nm<<" bigger than n of computed modes: "
             <<nmodes_total<<std::endl;
    DebugStop();
  }
  
  std::set<int64_t> dependent_connects;
  //this "dummy" connect will impose the restriction over the modal domain

  //TPZCompMesh::AllocateNewConnect(int nshape, int nstate, int order)
  const auto indep_con_id = scatt_mesh->AllocateNewConnect(nm,1,1);
  auto &indep_con = scatt_mesh->ConnectVec()[indep_con_id];
  //the geometric mesh will point to the scattering mesh
  scatt_mesh->LoadReferences();

      
  const auto &modal_block = modal_mesh->Block();
  for(auto modal_el : modal_mesh->ElementVec()){
    //select only volumetric elements, skip boundary els
    if(modal_el->Dimension() < modal_mesh->Dimension()){continue;}
    auto scatt_el = modal_el->Reference()->Reference();
    if(!scatt_el){DebugStop();}
    const auto ncon = scatt_el->NConnects();
    if(ncon != modal_el->NConnects()){
      DebugStop();
    }
    for(auto icon = 0; icon < ncon; icon++){
      auto &scatt_con = scatt_el->Connect(icon);
      auto &modal_con = modal_el->Connect(icon);
        
      const int64_t modal_dfseq = modal_con.SequenceNumber();
      const int64_t modal_pos = modal_block.Position(modal_dfseq);
      const int nshape = modal_block.Size(modal_dfseq);
      const auto dep_con_id = scatt_el->ConnectIndex(icon);
      //for lower order els we might not have shape functions
      //associated with a given connect
      if(nshape == 0){continue;}
      
      if(dependent_connects.count(dep_con_id) == 0 &&
         dirichlet_connects.count(dep_con_id) == 0 ){
        dependent_connects.insert(dep_con_id);
        scatt_con.AddDependency(dep_con_id, indep_con_id, modal_sol, modal_pos, 0, nshape, nm);
      }
    }
  }

  scatt_mesh->InitializeBlock();
  return indep_con_id;
}