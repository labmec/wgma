#ifndef _BCTYPES_HPP_
#define _BCTYPES_HPP_


namespace wgma{
  //! Enum describing boundary condition types
  enum class bctype {
    PEC,//!< Perfect electric conductor (dirichlet)
    PMC//!< Perfect magnetic conductor (neumann)
  };
};
#endif /* _BCTYPES_HPP_ */
