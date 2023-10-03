cmake_minimum_required(VERSION 3.13)

project (WaveguideAnalysis LANGUAGES CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# this option is useful for static analyzers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#this is the folder where additional findMYLIB.cmake modules must be put
#in order to use find_package(MYLIB)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# Finding the pz package if not in neopz buildtree
if(NOT CMAKE_IS_PZ_BUILDTREE)
  find_package(NeoPZ REQUIRED)
endif()

check_pz_opt(PZ_USING_MKL _PZ_HAS_MKL)
#enabling SLEPc library (enables PETSc, BLAS, LAPACK, ARPACK, SUPERLU and METIS)
option(USING_SLEPC "Whether the SLEPc library will be linked in" OFF)

if(NOT _PZ_HAS_MKL AND NOT USING_SLEPC)
  message(FATAL_ERROR "This project requires NeoPZ to "
    "have been configured using MKL or USING_SLEPC should be on")
endif()


add_subdirectory(wgma)
if(USING_SLEPC)
  include(cmake/enable_slepc.cmake)
  enable_slepc(wgma)
endif()

#now we check for json library
find_package(nlohmann_json 3.11.2 QUIET)
if(NOT TARGET nlohmann_json::nlohmann_json)
  #could not find package, downloading it
  include(FetchContent)

  # Optional: set this to ON if your target publicly links to nlohmann_json and needs to install() 
  # set(JSON_Install ON)

  if(NOT json_SOURCE_DIR)
    message(STATUS "Downloading nlohmann_json")
  endif()
  FetchContent_Declare(json
    GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
    GIT_PROGRESS TRUE
    GIT_SHALLOW TRUE
    GIT_TAG v3.11.2)

  FetchContent_MakeAvailable(json)

endif()
target_link_libraries(wgma PUBLIC nlohmann_json::nlohmann_json)

#let us add all the examples
add_subdirectory(examples)