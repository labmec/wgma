#include <cylindrical_pml.hpp>
#include <tpzautopointer.h>
#include <pzvec.h>
#include <regex>
namespace wgma::pml::cyl{
  /*TODO: write unit tests for this.
   it is not hard, just test all possible names and a few false ones
   must remember: test also with cartesian pml names*/
  TPZAutoPointer<wgma::pml::data>
  IdentifyAndSetupPML(const std::string &name,
                      const int id,
                      const int dim,
                      const CSTATE alphaR,
                      const CSTATE alphaZ)
  {

    TPZVec<std::string> pmlnames;
    //we need to iterate backwards!
    for (auto t : wgma::pml::cyl::typeReverseIterator() ) { //notice the parentheses!
      pmlnames.push_back("_"+wgma::pml::cyl::to_string(t));
    }
                
    for(int ipml = 0; ipml < pmlnames.size(); ipml++){
      const std::string pmlname = pmlnames[ipml];
      const auto rx = std::regex{ pmlname, std::regex_constants::icase };
      const bool test = std::regex_search(name, rx);
      if(test){
        wgma::pml::cyl::type type;
        bool found=false;
        for (auto t : wgma::pml::cyl::typeIterator() ) {
          if("_"+wgma::pml::cyl::to_string(t) == pmlname){
            type = t;
            found = true;
            break;
          }
        }
        if(!found){
          std::cout<<"could not identify pml "<<name<<std::endl;
          DebugStop();
        }
        auto pmldata = new wgma::pml::cyl::data;
        pmldata->ids = {id};
        pmldata->alphar = alphaR;
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