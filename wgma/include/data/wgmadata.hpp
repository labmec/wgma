#ifndef _WGMADATA_HPP_
#define _WGMADATA_HPP_

#include <map>
#include <string>

#include "cmeshtools.hpp"
#include "modetypes.hpp"
#include "nlohmann/json.hpp"

template <class T>
class TPZVec;
namespace wgma::data{

/**
   @brief Class containing all needed data for computing modal analysis of electromagnetic waveguides.
 */
  class WgmaData{
  public:
    //!Fills data from json instance and information from gmsh
    WgmaData(const nlohmann::json &, const TPZVec<std::map<std::string, int>>&);
    //!physical info of domains and bcs
    wgma::cmeshtools::PhysicalData modal_data;
    //!polynomial order of hcurl space
    int porder{-1};
    //!pattern that will be used to filter PMLs by name
    std::string pmlpattern{"invalid"};
  };


  class WgmaData1D : public WgmaData{
  public:
    //!Fills data from json instance and information from gmsh
    WgmaData1D(const nlohmann::json &, const TPZVec<std::map<std::string, int>>&);
    wgma::planarwg::mode mode{wgma::planarwg::mode::TE};
  };
};
#endif