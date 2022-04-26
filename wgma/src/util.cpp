#include "util.hpp"

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
};