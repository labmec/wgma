#include <cartesian_pml.hpp>
#include <tpzautopointer.h>
#include <pzvec.h>
#include <regex>

namespace wgma::pml::cart{
  /*TODO: write unit tests for this.
   it is not hard, just test all possible names and a few false ones
   must remember: test also with cylindrical pml names*/
  TPZAutoPointer<wgma::pml::data>
  IdentifyAndSetupPML(const std::string &name,
                      const int id,
                      const int dim,
                      const CSTATE alphaX,
                      const CSTATE alphaY,
                      const CSTATE alphaZ)
  {
    TPZVec<std::string> pmlnames;
    //we need to go backwards! otherwise we wont match the exact name
    for (auto t : wgma::pml::cart::typeReverseIterator() ) { //notice the parentheses!
      pmlnames.push_back("_"+wgma::pml::cart::to_string(t));
    }
                
    for(int ipml = 0; ipml < pmlnames.size(); ipml++){
      const std::string pmlname = pmlnames[ipml];
      const auto rx = std::regex{ pmlname, std::regex_constants::icase };
      const bool test = std::regex_search(name, rx);
      if(test){
        wgma::pml::cart::type type;
        bool found=false;
        for (auto t : wgma::pml::cart::typeIterator() ) {
          if("_"+wgma::pml::cart::to_string(t) == pmlname){
            type = t;
            found=true;
            break;
          }
        }
        if(!found){
          std::cout<<"could not identify pml "<<name<<std::endl;
          DebugStop();
        }
        auto pmldata = new wgma::pml::cart::data;
        pmldata->ids = {id};
        pmldata->alphax = alphaX;
        pmldata->alphay = alphaY;
        pmldata->alphaz = alphaZ;
        pmldata->t = type;
        pmldata->names = {name};
        pmldata->dim = dim;
        return pmldata;
      }
    }
    return nullptr;
  }
};