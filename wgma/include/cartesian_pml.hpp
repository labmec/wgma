#ifndef _CART_PML_HPP_
#define _CART_PML_HPP_

#include <pmltypes.hpp>
#include <util.hpp>

template<class T>
class TPZAutoPointer;

namespace wgma::pml::cart{

    TPZAutoPointer<wgma::pml::data> IdentifyAndSetupPML(const std::string &name,
                                                        const int id,
                                                        const int dim,
                                                        const CSTATE alphaX,
                                                        const CSTATE alphaY,
                                                        const CSTATE alphaZ);
    //! Enum describing pml attenuation direction
    enum class type{
        xp=0,//!< attenuates with increasing x
        yp,//!< attenuates with increasing y
        zp,//!< attenuates with increasing z
        xm,//!< attenuates with decreasing x
        ym,//!< attenuates with decreasing y
        zm,//!< attenuates with decreasing z,
        xpyp,//!< attenuates with increasing x and increasing y
        xmyp,//!< attenuates with decreasing x and increasing y
        xmym,//!< attenuates with decreasing x and decreasing y
        xpym,//!< attenuates with increasing x and decreasing y
        xpzp,//!< attenuates with increasing x and increasing z
        xmzp,//!< attenuates with decreasing x and increasing z
        ypzp,//!< attenuates with increasing y and increasing z
        ymzp,//!< attenuates with decreasing y and increasing z
        xpzm,//!< attenuates with increasing x and decreasing z
        xmzm,//!< attenuates with decreasing x and decreasing z
        ypzm,//!< attenuates with increasing y and decreasing z
        ymzm,//!< attenuates with decreasing y and decreasing z
        xpypzm,//!< attenuates with increasing x, increasing y and decreasing z
        xmypzm,//!< attenuates with decreasing x, increasing y and decreasing z
        xmymzm,//!< attenuates with decreasing x, decreasing y and decreasing z
        xpymzm,//!< attenuates with increasing x, decreasing y and decreasing z
        xpypzp,//!< attenuates with increasing x, increasing y and increasing z
        xmypzp,//!< attenuates with decreasing x, increasing y and increasing z
        xmymzp,//!< attenuates with decreasing x, decreasing y and increasing z
        xpymzp,//!< attenuates with increasing x, decreasing y and increasing z
    };

    typedef wgma::util::Iterator<type, type::xp, type::xpymzp> typeIterator;
    typedef wgma::util::ReverseIterator<type, type::xp, type::xpymzp> typeReverseIterator;
    
    inline std::string to_string(type t){
        switch(t){
        case type::xp: return "xp";
        case type::yp: return "yp";
        case type::zp: return "zp";
        case type::xm: return "xm";
        case type::ym: return "ym";
        case type::zm: return "zm";
        case type::xpyp: return "xpyp";
        case type::xmyp: return "xmyp";
        case type::xmym: return "xmym";
        case type::xpym: return "xpym";
        case type::xpzp: return "xpzp";
        case type::xmzp: return "xmzp";
        case type::ypzp: return "ypzp";
        case type::ymzp: return "ymzp";
        case type::xpzm: return "xpzm";
        case type::xmzm: return "xmzm";
        case type::ypzm: return "ypzm";
        case type::ymzm: return "ymzm";
        case type::xpypzm: return "xpypzm";
        case type::xmypzm: return "xmypzm";
        case type::xmymzm: return "xmymzm";
        case type::xpymzm: return "xpymzm";
        case type::xpypzp: return "xpypzp";
        case type::xmypzp: return "xmypzp";
        case type::xmymzp: return "xmymzp";
        case type::xpymzp: return "xpymzp";
        }
        return "error";
    }

    inline std::ostream& operator<<( std::ostream& out, const type& t ){
        return out << to_string(t);
    }

    //!Queries for pml attenuation in the x-direction.
    inline bool attx(type t){
        if( t == type::yp || t == type::ym || t == type::zp || t == type::zm) return false;
        else if (t == type::ymzm || t == type::ymzp ||
                 t == type::ypzm || t == type::ypzp) return false;
        else return true;
    }
    /** @brief Info for pml attenuation in the x-direction.
        @return A positive (negative)  value means that it attenuates with
        increasing (decreasing) x. 
        Zero means that it does not attenuate in the x direction.*/
    inline int xinfo(type t){
        if(!attx(t)){return 0;}
        else if (t == type::xp ||
                 t == type::xpyp || t == type::xpym ||
                 t == type::xpzp || t == type::xpzm ||
                 t == type::xpypzp || t == type::xpymzp ||
                 t == type::xpypzm || t == type::xpymzm){return 1;}
        else {return -1;}
    }
    //!Queries for pml attenuation in the y-direction.
    inline bool atty(type t){
        if( t == type::xp || t == type::xm || t == type::zp || t == type::zm ) return false;
        else if (t == type::xmzm || t == type::xmzp ||
                 t == type::xpzm || t == type::xpzp) return false;
        else return true;
    }
    /** @brief Info for pml attenuation in the x-direction.
        @return A positive (negative)  value means that it attenuates with
        increasing (decreasing) x. 
        Zero means that it does not attenuate in the x direction.*/
    inline int yinfo(type t){
        if(!atty(t)){return 0;}
        else if (t == type::yp ||
                 t == type::xpyp || t == type::xmyp ||
                 t == type::ypzm || t == type::ypzp ||
                 t == type::xmypzp || t == type::xpypzp ||
                 t == type::xmypzm || t == type::xpypzm){return 1;}
        else {return -1;}
    }

    //!Queries for pml attenuation in the z-direction.
    inline bool attz(type t){
        if( t == type::xp || t == type::xm || t == type::yp || t == type::ym) return false;
        else if (t == type::xmym || t == type::xmyp ||
                 t == type::xpym || t == type::xpyp) return false;
        else return true;
    }
    /** @brief Info for pml attenuation in the x-direction.
        @return A positive (negative)  value means that it attenuates with
        increasing (decreasing) x. 
        Zero means that it does not attenuate in the x direction.*/
    inline int zinfo(type t){
        if(!attz(t)){return 0;}
        else if (t == type::zp ||
                 t == type::xpzp || t == type::xmzp ||
                 t == type::ymzp || t == type::ypzp ||
                 t == type::xmypzp || t == type::xpypzp ||
                 t == type::xmymzp || t == type::xpymzp){return 1;}
        else {return -1;}
    }

    /** @brief Data structure for easier creation of PML regions.
        For automatic setting of PML (finding neighbours), leave attribute
        neigh empty and use only one id at a time.
        For pseudo-PML regions (i.e., periodic meshes), use attribute
        neigh for setting material properties for each PML region
    */
    struct data : public wgma::pml::data{
        type t{type::xm};
        CSTATE alphax{0};
        CSTATE alphay{0};
        CSTATE alphaz{0};
        explicit data(std::set<int> i, type tp, CSTATE ax, CSTATE ay, CSTATE az,
                      std::map<int,int> n = {{}}) :
            wgma::pml::data(i,n), t(tp), alphax(ax), alphay(ay), alphaz(az) {}
        data() = default;
    protected:
        void is_concrete_pml() override {}
    };
};
#endif /* _CART_PML_HPP_ */
