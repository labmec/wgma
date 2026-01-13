#include "util.hpp"

#include <rapidcsv.h>

#include <complex>
#include <filesystem>

namespace wgma::util{

  std::string ExtractPath(const std::string filepath){
 
    const auto pos = filepath.find_last_of("/\\");
 
    if(pos != filepath.npos){//found dir delimiter
      return filepath.substr(0,pos);
    }else{
      return "";
    }
  }
  
  void CreatePath(const std::string filepath){
 
    auto path = ExtractPath(filepath);
 
    if(path.size() != 0){
      CreatePath(path);
    }
    std::filesystem::create_directory(filepath);
  }


  std::complex<double> GetRefIndexFromCSV(const std::string &name, double wl){
    //we have tested already, file exists
    rapidcsv::Document doc(name, rapidcsv::LabelParams(-1, -1),rapidcsv::SeparatorParams(' '));
    auto wlvec = doc.GetColumn<double>(0);
    auto nvec = doc.GetColumn<double>(1);
    auto kvec = doc.GetColumn<double>(2);
    const auto nwl = wlvec.size();
    auto wlpos = std::lower_bound(wlvec.begin(),wlvec.end(),wl) - wlvec.begin();
    if(wlpos == 0 || wlpos == nwl){
      //we go back one position
      if(wlpos == nwl) wlpos--;
      //tolerance
      if(std::abs(wlvec[wlpos]-wl) < 1e-3){
        const std::complex<double> n {doc.GetCell<double>(wlpos,1),doc.GetCell<double> (wlpos,2)};
        return n;
      }else{
        std::cerr<<__PRETTY_FUNCTION__
                 <<"\nInvalid range for wl "<<wl
                 <<"\nMinimum wavelength in file "<<name<<" is "<<wlvec[wlpos]<<std::endl;
        std::abort();
      }
    }
    const auto wl1 = wlvec[wlpos-1];
    const std::complex<double> n1 {nvec[wlpos-1],-kvec[wlpos-1]};
    const auto wl2 = wlvec[wlpos];
    const std::complex<double> n2 {nvec[wlpos],-kvec[wlpos]};
    const std::complex<double> n = n1 + (wl-wl1)*(n2-n1)/(wl2-wl1);
    return n;
  }
};