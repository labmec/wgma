
# - Try to find SLEPC
# Once done this will define
#
#  SLEPC_FOUND            - system has SLEPc
#  SLEPC_INCLUDE_DIRS     - include directories for SLEPc
#  SLEPC_LIBRARY_DIRS     - library directories for SLEPc
#  SLEPC_LIBARIES         - libraries for SLEPc
#  SLEPC_STATIC_LIBARIES  - ibraries for SLEPc (static linking, undefined if not required)
#  SLEPC_VERSION          - version of SLEPc
#  SLEPC_VERSION_MAJOR    - First number in SLEPC_VERSION
#  SLEPC_VERSION_MINOR    - Second number in SLEPC_VERSION
#  SLEPC_VERSION_SUBMINOR - Third number in SLEPC_VERSION


#=============================================================================
# Copyright (C) 2010-2016 Garth N. Wells, Anders Logg and Johannes Ring
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

unset(SLEPC_FOUND CACHE)
unset(SLEPC_INCLUDE_DIRS CACHE)
unset(SLEPC_LIBRARY_DIRS CACHE)
unset(SLEPC_LIBARIES CACHE)
unset(SLEPC_STATIC_LIBARIES CACHE)
unset(SLEPC_VERSION CACHE)
unset(SLEPC_VERSION_MAJOR CACHE)
unset(SLEPC_VERSION_MINOR CACHE)
unset(SLEPC_VERSION_SUBMINOR CACHE)

# Load CMake pkg-config module
find_package(PkgConfig REQUIRED)

# Find SLEPc pkg-config file
set(ENV{PKG_CONFIG_PATH} "${SLEPC_DIR}/${PETSC_ARCH}/lib/pkgconfig:${SLEPC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
set(ENV{PKG_CONFIG_PATH} "${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:${PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
set(ENV{PKG_CONFIG_PATH} "${PETSC_DIR}/${PETSC_ARCH}:${PETSC_DIR}:$ENV{PKG_CONFIG_PATH}")
pkg_search_module(SLEPC slepc)

# Extract major, minor, etc from version string
if (SLEPC_VERSION)
  message(STATUS "Found SLEPc version ${SLEPC_VERSION}")
  string(REPLACE "." ";" VERSION_LIST ${SLEPC_VERSION})
  list(GET VERSION_LIST 0 SLEPC_VERSION_MAJOR)
  list(GET VERSION_LIST 1 SLEPC_VERSION_MINOR)
  list(GET VERSION_LIST 2 SLEPC_VERSION_SUBMINOR)
endif()

# Configure SLEPc IMPORT (this involves creating an 'imported' target
# and attaching 'properties')
if (SLEPC_FOUND AND NOT TARGET SLEPC::slepc)
  add_library(SLEPC::slepc INTERFACE IMPORTED)

  # Add include paths
  set_property(TARGET SLEPC::slepc PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${SLEPC_INCLUDE_DIRS})

  # Add libraries
  unset(_libs)
  if(SLEPC_LIBRARIES)
    message(STATUS "SLEPC libs:")
  endif()
  foreach (lib ${SLEPC_LIBRARIES})
    unset(LIB_${lib} CACHE)
    #find_library(LIB_${lib} NAMES ${lib} PATHS ${SLEPC_LIBRARY_DIRS} NO_DEFAULT_PATH)
    find_library(LIB_${lib} ${lib} HINTS ${SLEPC_LIBRARY_DIRS})
    list(APPEND _libs ${LIB_${lib}})
    message(STATUS "LIB_${lib} ${LIB_${lib}}")
  endforeach()
  set_property(TARGET SLEPC::slepc PROPERTY INTERFACE_LINK_LIBRARIES "${_libs}")
endif()

if (SLEPC_FOUND AND NOT TARGET SLEPC::slepc_static)
  add_library(SLEPC::slepc_static INTERFACE IMPORTED)

  # Add libraries (static)
  unset(_libs)
  foreach (lib ${SLEPC_STATIC_LIBRARIES})
    find_library(LIB_${lib} ${lib} HINTS ${SLEPC_STATIC_LIBRARY_DIRS} NO_DEFAULT_PATH)
    list(APPEND _libs ${LIB_${lib}})
  endforeach()
  set_property(TARGET SLEPC::slepc_static PROPERTY INTERFACE_LINK_LIBRARIES "${_libs}")

endif()


# Compile and run test
if (SLEPC_FOUND)
  # Create SLEPc test program
  set(SLEPC_TEST_LIB_CPP
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/slepc_test_lib.cpp")
  file(WRITE ${SLEPC_TEST_LIB_CPP} "
#include \"petsc.h\"
#include \"slepceps.h\"
int main()
{
  PetscErrorCode ierr;
  int argc = 0;
  char** argv = NULL;
  ierr = SlepcInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
  EPS eps;
  ierr = EPSCreate(PETSC_COMM_SELF, &eps); CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);
  ierr = SlepcFinalize(); CHKERRQ(ierr);
  return 0;
}
")
  # Try to run test program (shared linking)
  try_run(
    SLEPC_TEST_LIB_EXITCODE
    SLEPC_TEST_LIB_COMPILED
    ${CMAKE_CURRENT_BINARY_DIR}
    ${SLEPC_TEST_LIB_CPP}
    CMAKE_FLAGS
    "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
    LINK_LIBRARIES PETSC::petsc SLEPC::slepc  ${CMAKE_REQUIRED_LIBRARIES}
    COMPILE_OUTPUT_VARIABLE SLEPC_TEST_LIB_COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE SLEPC_TEST_LIB_OUTPUT
    )
  if (SLEPC_TEST_LIB_COMPILED AND SLEPC_TEST_LIB_EXITCODE EQUAL 0)
    message(STATUS "Test SLEPC_TEST_RUNS with shared library linking - Success")
    set(SLEPC_TEST_RUNS TRUE)

    # Static libraries not required, so unset
    set_property(TARGET SLEPC::slepc_static PROPERTY INTERFACE_LINK_LIBRARIES )
  else()
    message(STATUS "Test SLEPC_TEST_RUNS with shared library linking - Failed")
    message(STATUS "Tried to compile with\n${SLEPC_TEST_LIB_COMPILE_OUTPUT}")
    # Try to run test program (static linking)
    try_run(
      SLEPC_TEST_STATIC_LIBS_EXITCODE
      SLEPC_TEST_STATIC_LIBS_COMPILED
      ${CMAKE_CURRENT_BINARY_DIR}
      ${SLEPC_TEST_LIB_CPP}
      CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
      LINK_LIBRARIES PETSC::petsc SLEPC::slepc_static ${CMAKE_REQUIRED_LIBRARIES}
      COMPILE_OUTPUT_VARIABLE SLEPC_TEST_STATIC_LIBS_COMPILE_OUTPUT
      RUN_OUTPUT_VARIABLE SLEPC_TEST_STATIC_LIBS_OUTPUT
      )
    if (SLEPC_TEST_STATIC_LIBS_COMPILED AND SLEPC_TEST_STATIC_LIBS_EXITCODE EQUAL 0)
      message(STATUS "Test SLEPC_TEST__RUNS with static linking - Success")
      set(SLEPC_TEST_RUNS TRUE)
    else()
      message(STATUS "Test SLEPC_TETS_RUNS with static linking - Failed")
      set(SLEPC_TEST_RUNS FALSE)
    endif()
  endif()
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
if (SLEPC_FOUND)
  find_package_handle_standard_args(SLEPc
    REQUIRED_VARS SLEPC_FOUND SLEPC_TEST_RUNS
    VERSION_VAR SLEPC_VERSION
    FAIL_MESSAGE "SLEPc could not be configured.")
else()
  find_package_handle_standard_args(SLEPc
    REQUIRED_VARS SLEPC_FOUND
    FAIL_MESSAGE "SLEPc could not be found. Be sure to set SLEPC_DIR.")
endif()
