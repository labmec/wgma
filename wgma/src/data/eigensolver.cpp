#include <data/eigensolver.hpp>
#include <nlohmann/json.hpp>
#include "pzerror.h"

using nlohmann::json;
using namespace wgma::data;

Eigensolver::Eigensolver(const json &json){
  if(json.contains("eigensolver") == 0){
    std::cout<<"could not find key geom in json file! Aborting..."
             <<std::endl;
    DebugStop();
  }
  
  auto eigendata = json["eigensolver"];
  
  
  if(eigendata.contains("slepc")){
    slepc=eigendata["slepc"];
  }
  if(eigendata.contains("krylovdim")){
    krylovdim=eigendata["krylovdim"];
  }
  neigenpairs=eigendata["neigenpairs"];
  target=eigendata["target"];
  sort=eigendata["sorting"];
}