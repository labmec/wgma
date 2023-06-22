#ifndef _PMLTYPES_HPP_
#define _PMLTYPES_HPP_

#include <pzreal.h>
#include <set>
#include <map>

namespace wgma::pml{
  //!generic data for PML types
  struct data{
    int dim{-1};
    std::set<int> ids{-100};
    std::set<std::string> names{"no_name"};//<name from the region associated with PML
    std::map<int,int> neigh{};
    explicit data(std::set<int> i, std::map<int,int> n={{}}) :
      ids(i), neigh(n) {}
    data() = default;
    //we need to declare virtual dtor, so rule of 5
    virtual ~data()= default;
    data(const data&) = default;
    data& operator=(const data&) = default;
    data(data&&) = default;
    data& operator=(data&&) = default;
  protected:
    virtual void is_concrete_pml() = 0;
  };
}

#include <cartesian_pml.hpp>
#include <cylindrical_pml.hpp>

#endif /* _PMLTYPES_HPP_ */
