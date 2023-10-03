#ifndef _JSON_UTIL_HPP_
#define _JSON_UTIL_HPP_
#include <nlohmann/json.hpp>
#include <complex>




namespace std {
    
  using json = nlohmann::json;
  template< class T > void to_json(json &j, const std::complex< T > &p) {
    j = json {p.real(), p.imag()};
  }
    
  template< class T > void from_json(const json &j, std::complex< T > &p) {
    if(j.is_number()){
      p = j.get<T>();
    }
    else{
      if(j.size()!=2){
        //error
      }
      p.real(j.at(0));
      p.imag(j.at(1));
    }
  }
}
#endif