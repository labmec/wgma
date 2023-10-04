#include <data/wgmadata.hpp>
//#include "tpzautopointer.h"
#include <json_util.hpp>
#include <pzerror.h>
#include <pzmanvector.h>

using namespace wgma::data;
using nlohmann::json;

WgmaData::WgmaData(const json &simdata,
                   const TPZVec<std::map<std::string, int>>&gmshmats)
{

  porder = simdata["porder"];
  std::map<std::string, std::pair<CSTATE, CSTATE>> modal_mats;
  std::map<std::string, wgma::bc::type> modal_bcs;
  
  TPZManVector<CSTATE,3> alphapml;
  //now we read material info from json
  for(auto &[key,val] : simdata["materials"].items()){
    //reffractive index
    if(val.contains("isotropic") && val["isotropic"] == false){
      std::cout<<"not yet implemented!"<<std::endl;
      DebugStop();
    }
    const CSTATE nmat =  val["n"];
    const CSTATE urmat = val["ur"];
    modal_mats[key] = {nmat*nmat,urmat};
  }
  //now we read boundary info from json
  for(auto &[key,val] : simdata["bcs"].items()){
    //allowed values: PEC,PMC,PERIODIC
    modal_bcs[key] = wgma::bc::from_string(val.get<std::string>());
  }
  const int modaldim = simdata["dim"];

  if(simdata.contains("pml")){
    auto pmldata = simdata["pml"];
    if(pmldata.contains("alphax")){
      //cartesian pml
      alphapml.push_back(pmldata["alphax"]);
      if(pmldata.contains("alphay")){
        alphapml.push_back(pmldata["alphay"]);
      }
      if(pmldata.contains("alphaz")){
        alphapml.push_back(pmldata["alphaz"]);
      }
    }else if(pmldata.contains("alphar")){
      //cylindrical pml (todo: what about spherical?)
      alphapml.push_back(pmldata["alphar"]);
      if(pmldata.contains("alphaz")){
        alphapml.push_back(pmldata["alphaz"]);
      }
    }else{
      DebugStop();
    }
    pmlpattern = pmldata.value("pattern","*");
  }else{
    alphapml={};
  }
  //now we fill our modal_data attribute
  wgma::cmeshtools::SetupGmshMaterialData(gmshmats, modal_mats, modal_bcs,
                                          alphapml, modal_data,
                                          modaldim);
}

WgmaData1D::WgmaData1D(const nlohmann::json &data,
                       const TPZVec<std::map<std::string, int>>&gmshmats)
  :WgmaData(data,gmshmats)
{
  mode = wgma::planarwg::string_to_mode(data.value("mode","Invalid"));
  if(mode == wgma::planarwg::mode::Invalid){
    PZError<<"could not read mode from json file. value :\n "
           <<data["mode"]<<"\nAborting..."<<std::endl;
    DebugStop();
  }
}