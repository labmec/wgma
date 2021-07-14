#ifndef _PMLTYPES_HPP_
#define _PMLTYPES_HPP_

namespace wgma{
  //! Enum describing pml attenuation direction
  enum class pmltype{
    xp=0,//!< attenuates with increasing x
    yp,//!< attenuates with increasing y
    xm,//!< attenuates with decreasing x
    ym,//!< attenuates with decreasing y
    xpyp,//!< attenuates with increasing x and y
    xmyp,//!< attenuates with decreasing x and increasing y
    xmym,//!< attenuates with decreasing x and y
    xpym//!< attenuates with increasing x and decreasing y
  };
}
#endif /* _PMLTYPES_HPP_ */
