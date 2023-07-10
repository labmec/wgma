#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#include <string>
#include <type_traits>
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

  //! Allows for iterating over enums taken from https://stackoverflow.com/a/31836401
  template < typename C, C beginVal, C endVal>
  class Iterator {
    typedef typename std::underlying_type<C>::type val_t;
    int val;
  public:
    Iterator(const C & f) : val(static_cast<val_t>(f)) {}
    Iterator() : val(static_cast<val_t>(beginVal)) {}
    Iterator operator++() {
      ++val;
      return *this;
    }
    C operator*() { return static_cast<C>(val); }
    Iterator begin() { return *this; } //default ctor is good
    Iterator end() {
      static const Iterator endIter=++Iterator(endVal); // cache it
      return endIter;
    }
    bool operator!=(const Iterator& i) { return val != i.val; }
  };
  //! Reverse iterator for enums
  template < typename C, C beginVal, C endVal>
  class ReverseIterator {
    typedef typename std::underlying_type<C>::type val_t;
    int val;
  public:
    ReverseIterator(const C & f) : val(static_cast<val_t>(f)) {}
    ReverseIterator() : val(static_cast<val_t>(endVal)) {}
    ReverseIterator operator++() {
      --val;
      return *this;
    }
    C operator*() { return static_cast<C>(val); }
    ReverseIterator begin() { return *this; } //default ctor is good
    ReverseIterator end() {
      static const ReverseIterator endIter=++ReverseIterator(beginVal); // cache it
      return endIter;
    }
    bool operator!=(const ReverseIterator& i) { return val != i.val; }
  };
};

#endif /* _UTIL_HPP_ */


