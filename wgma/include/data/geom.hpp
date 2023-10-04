#ifndef _GEOMDATA_HPP_
#define _GEOMDATA_HPP_

#include <nlohmann/json.hpp>
#include <pzreal.h>
#include <string>

namespace wgma::data{
  /**
     @brief Class containing all needed data for reading a gmsh file
  */
  class Geom{
  public:
    /** @brief Fills data from json instance looking for "geom" key
        @note Scale can be set "k0" to correspond to 2pi/lambda
     */
    Geom(const nlohmann::json &, const STATE lambda=1.55);
    //!.msh file to be read
    std::string meshfile;
    //!.csv file in which curved entities information is stored
    std::string arcfile;
    //!verbosity lvl during reading of .msh file
    bool verbosity{false};
    //!whether to print NeoPZ geometric mesh
    bool print{false};
    //!whether to represent curved entities with non-linear mapping
    bool arc3d{false};
    //!additional scaling factor applied to geometry
    REAL scale{1};
  };
};
#endif