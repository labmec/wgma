#ifndef _PMLTYPES_HPP_
#define _PMLTYPES_HPP_

#include <pzreal.h>

namespace wgma{
  namespace pml{
    //! Enum describing pml attenuation direction
    enum class type{
      xp=0,//!< attenuates with increasing x
      yp,//!< attenuates with increasing y
      xm,//!< attenuates with decreasing x
      ym,//!< attenuates with decreasing y
      xpyp,//!< attenuates with increasing x and y
      xmyp,//!< attenuates with decreasing x and increasing y
      xmym,//!< attenuates with decreasing x and y
      xpym//!< attenuates with increasing x and decreasing y
    };

    std::string to_string(type t){
      switch(t){
      case type::xp: return "xp";
      case type::yp: return "yp";
      case type::xm: return "xm";
      case type::ym: return "ym";
      case type::xpyp: return "xpyp";
      case type::xmyp: return "xmyp";
      case type::xmym: return "xmym";
      case type::xpym: return "xpym";
      }
    }

    std::ostream& operator<<( std::ostream& out, const type& t ){
      return out << to_string(t);
    }

    //!Queries for pml attenuation in the x-direction.
    bool attx(type t){
      if( t == type::yp || t == type::ym) return false;
      else return true;
    }
    /** @brief Info for pml attenuation in the x-direction.
        @return A positive (negative)  value means that it attenuates with
    increasing (decreasing) x. 
    Zero means that it does not attenuate in the x direction.*/
    int xinfo(type t){
      if( t == type::yp || t == type::ym) return 0;
      else if(t == type::xp || t == type:: xpyp || t == type::xpym) return 1;
      else return -1;
    }
    //!Queries for pml attenuation in the y-direction.
    bool atty(type t){
      if( t == type::xp || t == type::xm) return false;
      else return true;
    }
    /** @brief Info for pml attenuation in the x-direction.
        @return A positive (negative)  value means that it attenuates with
    increasing (decreasing) x. 
    Zero means that it does not attenuate in the x direction.*/
    int yinfo(type t){
      if( t == type::xp || t == type::xm) return 0;
      else if(t == type::yp || t == type:: xpyp || t == type::xmyp) return 1;
      else return -1;
    }

    //! Data structure for easier creation of PML regions
    struct data{
      type t;
      STATE alpha;
      int id;
      explicit data(int i, type tp, STATE a) : id(i), t(tp), alpha(a) {}
      data() : id(-100), alpha(0), t(type::xm) {}
    };
  };
  
}
#endif /* _PMLTYPES_HPP_ */
