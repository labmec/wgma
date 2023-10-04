#ifndef _MODE_TYPES_HPP_
#define _MODE_TYPES_HPP_

#include <string>
#include <map>

namespace wgma::planarwg{
  enum class mode{Invalid,TE, TM};//< modes approximated in planar waveguides

  //! for easier output
  inline std::string mode_to_string(mode m){
    switch(m){
    case mode::Invalid:
      return "Invalid";
    case mode::TE:
      return "TE";
    case mode::TM:
      return "TM";
    }
  }
  //! for easier input
  inline mode string_to_mode(const std::string name){
    static const std::map<std::string,mode> stringvals{
      {"TE",mode::TE},
      {"TM",mode::TM}
    };

    auto itr = stringvals.find(name);
    if( itr != stringvals.end() ) {
      return itr->second;
    }
    return mode::Invalid; 
  }
};

#endif