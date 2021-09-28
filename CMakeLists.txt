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
add_subdirectory(wr90)
add_subdirectory(stepfiber)