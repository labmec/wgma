#include <cylindrical_pml.hpp>
#include <tpzautopointer.h>
#include <pzvec.h>
#include <regex>
namespace wgma::pml::cyl{
  TPZAutoPointer<wgma::pml::data>
  IdentifyAndSetupPML(const std::string &name,
                      const int id,
                      const int dim,
                      const STATE alphaR,
                      const STATE alphaZ)
  {

    TPZVec<std::string> pmlnames;
    //we need to iterate backwards!
    for (auto t : wgma::pml::cyl::typeReverseIterator() ) { //notice the parentheses!
      pmlnames.push_back(wgma::pml::cyl::to_string(t));
    }
                
    for(int ipml = 0; ipml < pmlnames.size(); ipml++){
      const std::string pmlname = pmlnames[ipml];
      const auto rx = std::regex{ pmlname, std::regex_constants::icase };
      const bool test = std::regex_search(name, rx);
      if(test){
        wgma::pml::cyl::type type;
        for (auto t : wgma::pml::cyl::typeIterator() ) {
          if(wgma::pml::cyl::to_string(t) == pmlname){
            type = t;
            break;
          }
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