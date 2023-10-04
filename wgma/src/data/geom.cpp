#include <data/geom.hpp>
#include <nlohmann/json.hpp>
#include "pzerror.h"

using nlohmann::json;
using namespace wgma::data;

Geom::Geom(const json &json, const STATE lambda){
  if(json.contains("geom") == 0){
    std::cout<<"could not find key geom in json file! Aborting..."
             <<std::endl;
    DebugStop();
  }

  auto geomdata = json["geom"];
  
  if(geomdata.contains("meshfile")){
    meshfile=geomdata["meshfile"];
  }
  if(geomdata.contains("print_gmsh")){
    print=geomdata["print_gmsh"];
  }
  if(geomdata.contains("arc3d")){
    arc3d=geomdata["arc3d"];
  }
  if(geomdata.contains("arcfile")){
    arcfile=geomdata["arcfile"];
  }
  if(geomdata.contains("scale")){
    if(geomdata["scale"].is_string()){
      if(geomdata["scale"]=="k0"){
        scale = lambda/(2*M_PI);
      }else{
        std::cout<<"could not interpret scale value "<<geomdata["scale"]
                 <<"\nAborting..."<<std::endl;
        DebugStop();
      }
    }else{
      scale=geomdata["scale"];
    }
  }
}