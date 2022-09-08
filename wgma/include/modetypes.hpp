#ifndef _MODE_TYPES_HPP_
#define _MODE_TYPES_HPP_

#include <string>

namespace wgma::planarwg{
  enum class mode{TE, TM};//< modes approximated in planar waveguides

  //! for easier output
  inline std::string mode_to_string(mode m){
    switch(m){
    case mode::TE:
      return "TE";
    case mode::TM:
      return "TM";
    }
  }
};

#endif