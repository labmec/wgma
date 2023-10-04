#ifndef _EIGENDATA_HPP_
#define _EIGENDATA_HPP_

#include <nlohmann/json.hpp>
#include <pzreal.h>
#include <TPZEigenSort.h>
#include <string>


enum class TPZEigenSort;

namespace wgma::data{
  
  /**
     @brief Class that stores basic eigensolver info read from json file

     json format:

     "eigensolver":{
       "slepc":true,//bool
       "krylovdim":300,//int
       "neigenpairs":150,//int
       "target":2.4025,
       "sorting":"TargetRealPart"
       //options are: (see TPZEigenSort.h)
       //"AbsAscending" Ascending magnitude
       //"AbsDescending" Descending magnitude
       //"RealAscending" Ascending real part
       //"RealDescending" Descending real part
       //"ImagAscending" Ascending imaginary part
       //"ImagDescending" Descending imaginary part
       //"TargetRealPart" Real part closest to target
       //"TargetImagPart" Imaginary part closest to target
       //"TargetMagnitude"Magnitude closest to target
     }
  */
  class Eigensolver{
  public:
    /** @brief Fills data from json instance looking for "solver" key
     */
    Eigensolver(const nlohmann::json &);
    //!whether to use slepc solver (only for generalised evp)
    bool slepc{false};
    //!(maximum) dimension of krylov subspace
    int krylovdim{1};
    //!number of sought eigenpairs
    int neigenpairs{1};
    //!value around which eigenvalues are sought
    CSTATE target{-1};
    //!sorting rule for computed eigenvalues
    TPZEigenSort sort{TPZEigenSort::Invalid};
  };
};
#endif