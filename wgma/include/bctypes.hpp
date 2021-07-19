#ifndef _BCTYPES_HPP_
#define _BCTYPES_HPP_


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

    //! Data structure for easier creation of PML regions
    struct data{
      type t;
      int id;
    };
  };
};
#endif /* _BCTYPES_HPP_ */
