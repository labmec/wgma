cmake_minimum_required(VERSION 3.13)

project (WaveguideAnalysis LANGUAGES CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# this option is useful for static analyzers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Finding the pz package if not in neopz buildtree
if(NOT CMAKE_IS_PZ_BUILDTREE)
  find_package(NeoPZ REQUIRED)
endif()

check_pz_opt(PZ_USING_LAPACK _PZ_HAS_LAPACK)

if(NOT _PZ_HAS_LAPACK)
  message(FATAL_ERROR "This project requires NeoPZ to "
    "have been configured using LAPACK")
endif()


add_subdirectory(wgma)
add_subdirectory(wr90)
add_subdirectory(stepfiber)