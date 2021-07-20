#ifndef _BCTYPES_HPP_
#define _BCTYPES_HPP_
#include <ostream>

namespace wgma{
  namespace bc{
    //! Enum describing boundary condition types
    enum class type {
      PEC,//!< Perfect electric conductor (dirichlet)
      PMC//!< Perfect magnetic conductor (neumann)
    };

    int to_int(type t){
      switch(t){
      case type::PEC: return 0;
      case type::PMC: return 1;
      }
    }

    std::string to_string(type t){
      switch(t){
      case type::PEC: return "PEC";
      case type::PMC: return "PMC";
      }
    }

    std::ostream& operator<<( std::ostream& out, const type& t ){
      return out << to_string(t);
    }

    //! Data structure for easier creation of PML regions
    struct data{
      type t;
      int id;
      explicit data(int i, type tp) : id(i), t(tp) {}
      data() : t(type::PEC), id(-100) {}
    };
  };
};
#endif /* _BCTYPES_HPP_ */
