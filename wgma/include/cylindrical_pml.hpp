#ifndef _CYL_PML_HPP_
#define _CYL_PML_HPP_

#include <pmltypes.hpp>
#include <util.hpp>

template<class T>
class TPZAutoPointer;

namespace wgma::pml::cyl{

    TPZAutoPointer<wgma::pml::data>
    IdentifyAndSetupPML(const std::string &name,
                        const int id,
                        const int dim,
                        const STATE alphaR,
                        const STATE alphaZ);
    
    //! Enum describing pml attenuation direction
    enum class type{
        rp=0,//!< attenuates with increasing radius
        zp,//!< attenuates with increasing z
        zm,//!< attenuates with decreasing z,
        rpzp,//!< attenuates with increasing radius and increasing z
        rpzm//!< attenuates with increasing radius and decreasing z
    };

    typedef wgma::util::Iterator<type, type::rp, type::rpzm> typeIterator;
    
    inline std::string to_string(type t){
        switch(t){
        case type::rp: return "rp";
        case type::zp: return "zp";
        case type::zm: return "zm";
        case type::rpzp: return "rpzp";
        case type::rpzm: return "rpzm";
        }
    }
    
    inline std::ostream& operator<<( std::ostream& out, const type& t ){
        return out << to_string(t);
    }

    //!Queries for pml attenuation in the radial direction.
    inline bool attr(type t){
        if( t == type::rp || t == type::rpzp || t == type::rpzm)
        {return true;}
        else{ return false;}
    }
    /** @brief Info for pml attenuation in the  radial direction.*/
    inline int rinfo(type t){
        if(!attr(t)){return 0;}
        else {return 1;}
    }
    //!Queries for pml attenuation in the z-direction.
    inline bool attz(type t){
        if( t == type::zp || t == type::zm || t == type::rpzp || t == type::rpzm){
            return true;
        }else {return false;}
    }
    /** @brief Info for pml attenuation in the x-direction.
        @return A positive (negative)  value means that it attenuates with
        increasing (decreasing) x. 
        Zero means that it does not attenuate in the x direction.*/
    inline int zinfo(type t){
        if(!attz(t)){return 0;}
        else if (t == type::zp || t == type::rpzp){return 1;}
        else {return -1;}
    }

    /** @brief Data structure for easier creation of PML regions.
        For automatic setting of PML (finding neighbours), leave attribute
        neigh empty and use only one id at a time.
        For pseudo-PML regions (i.e., periodic meshes), use attribute
        neigh for setting material properties for each PML region
    */
    struct data :public wgma::pml::data{
        type t{type::rp};
        STATE alphar{0};
        STATE alphaz{0};
        explicit data(std::set<int> i, type tp, STATE ar, STATE az,
                      std::map<int,int> n = {{}}) :
            wgma::pml::data(i,n),t(tp), alphar(ar), alphaz(az) {}
        data() = default;
    protected:
        void is_concrete_pml() override {}
    };
};
#endif /* _CYL_PML_HPP_ */
