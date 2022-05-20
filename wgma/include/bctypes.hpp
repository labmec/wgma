#ifndef _BCTYPES_HPP_
#define _BCTYPES_HPP_
#include <ostream>

namespace wgma{
  namespace bc{
    //! Enum describing boundary condition types
    enum class type {
      PEC,     //!< Perfect electric conductor (dirichlet)
      PMC,     //!< Perfect magnetic conductor (neumann)
      PERIODIC //!< Periodic BC (must come in pairs)
    };

    inline int to_int(type t) {
      switch (t) {
      case type::PEC:
        return 0;
      case type::PMC:
        return 1;
      case type::PERIODIC:
        return 2;
      }
    }

    inline std::string to_string(type t){
      switch(t){
      case type::PEC: return "PEC";
      case type::PMC: return "PMC";
      case type::PERIODIC: return "PERIODIC";
      }
    }

    inline std::ostream& operator<<( std::ostream& out, const type& t ){
      return out << to_string(t);
    }

    //! Data structure for easier creation of Boundary Conditions
    struct data{
      std::string name{"none"};
      int id{-100};
      type t{type::PEC};
      int volid{-100};
      explicit data(int i, type tp, int v = 0) : id(i), t(tp), volid(v) {}
      data() {}
    };
    
  };
};
#endif /* _BCTYPES_HPP_ */
