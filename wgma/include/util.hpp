#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include <string>

namespace wgma::util{

  /**
     @brief Extracts path from file name
     @note string will be split both looking for '/' and '\'.
  */
  std::string ExtractPath(const std::string filepath);
  /**
     @brief Recursively creates all non-existing directories of given path.
     @note string will be split both looking for '/' and '\'.
  */
  void CreatePath(const std::string path);
};

#endif /* _UTIL_HPP_ */


