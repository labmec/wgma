# - Try to find PETSc
# Once done this will define
#
#  PETSC_FOUND             - system has PETSc
#  PETSC_INCLUDE_DIRS      - include directories for PETSc
#  PETSC_LIBRARY_DIRS      - library directories for PETSc
#  PETSC_LIBRARIES         - libraries for PETSc
#  PETSC_STATIC_LIBRARIES  - libraries for PETSc (static linking, undefined if not required)
#  PETSC_VERSION           - version for PETSc
#  PETSC_VERSION_MAJOR     - First number in PETSC_VERSION
#  PETSC_VERSION_MINOR     - Second number in PETSC_VERSION
#  PETSC_VERSION_SUBMINOR  - Third number in PETSC_VERSION
#  PETSC_SCALAR_SIZE       - sizeof(PetscScalar)
#
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
unset(PETSC_FOUND CACHE)
unset(PETSC_INCLUDE_DIRS CACHE)
unset(PETSC_LIBRARY_DIRS CACHE)
unset(PETSC_LIBRARIES CACHE)
unset(PETSC_STATIC_LIBRARIES CACHE)
unset(PETSC_VERSION CACHE)
unset(PETSC_VERSION_MAJOR CACHE)
unset(PETSC_VERSION_MINOR CACHE)
unset(PETSC_VERSION_SUBMINOR CACHE)
unset(PETSC_SCALAR_SIZE CACHE)
unset(HAVE_PETSC_SCALAR_SIZE CACHE)
unset(PETSC_INT_SIZE CACHE)
unset(HAVE_PETSC_INT_SIZE CACHE)

# Outline:
# 1. Get flags from PETSc-generated pkg-config file
# 2. Test compile and run program using shared library linking
# 3. If shared library linking fails, test with static library linking

# Load pkg-config module (provided by CMake)
find_package(PkgConfig REQUIRED)
find_package(MPI REQUIRED CXX)
# Find PETSc pkg-config file.
set(ENV{PKG_CONFIG_PATH} "${PETSC_DIR}/lib/pkgconfig:${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
pkg_search_module(PETSC PETSc)

# Extract major, minor, etc from version string
if (PETSC_VERSION)
  message(STATUS "Found PETSc version ${PETSC_VERSION}")
  string(REPLACE "." ";" VERSION_LIST ${PETSC_VERSION})
  list(GET VERSION_LIST 0 PETSC_VERSION_MAJOR)
  list(GET VERSION_LIST 1 PETSC_VERSION_MINOR)
  list(GET VERSION_LIST 2 PETSC_VERSION_SUBMINOR)
endif()

# Configure PETSc IMPORT (this involves creating an 'imported' target
# and attaching 'properties')
if (PETSC_FOUND AND NOT TARGET PETSC::petsc)
  add_library(PETSC::petsc INTERFACE IMPORTED)
  # Add include paths
  set_property(TARGET PETSC::petsc PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${PETSC_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH})
  # Add libraries
  unset(_libs)
  foreach (lib ${PETSC_LIBRARIES})
    find_library(LIB_${lib} NAMES ${lib} PATHS ${PETSC_LIBRARY_DIRS} NO_DEFAULT_PATH)
    list(APPEND _libs ${LIB_${lib}})
  endforeach()
  set_property(TARGET PETSC::petsc PROPERTY INTERFACE_LINK_LIBRARIES "${_libs}" ${MPI_CXX_LIBRARIES})
endif()

## Configure PETSc 'static' IMPORT (this involves creating an
## 'imported' target and attaching 'properties')
if (PETSC_FOUND AND NOT TARGET PETSC::petsc_static)
  add_library(PETSC::petsc_static INTERFACE IMPORTED)

  # Add libraries (static)
  unset(_libs)
  foreach (lib ${PETSC_STATIC_LIBRARIES})
    find_library(LIB_${lib} ${lib} HINTS ${PETSC_STATIC_LIBRARY_DIRS})
    list(APPEND _libs ${LIB_${lib}})
  endforeach()
  set_property(TARGET PETSC::petsc_static PROPERTY INTERFACE_LINK_LIBRARIES "${_libs}")
endif()

# Attempt to build and run PETSc test program
if (PETSC_FOUND)

  # Create PETSc test program
  set(PETSC_TEST_LIB_CPP
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/petsc_test_lib.cpp")
  file(WRITE ${PETSC_TEST_LIB_CPP} "
#include \"petscts.h\"
#include \"petsc.h\"
int main()
{
  PetscErrorCode ierr;
  TS ts;
  int argc = 0;
  char** argv = NULL;
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
")

  # Try to run test program (shared linking)
  try_run(
    PETSC_TEST_LIB_EXITCODE
    PETSC_TEST_LIB_COMPILED
    ${CMAKE_CURRENT_BINARY_DIR}
    ${PETSC_TEST_LIB_CPP}
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
    LINK_LIBRARIES PETSC::petsc ${CMAKE_REQUIRED_LIBRARIES}
    COMPILE_OUTPUT_VARIABLE PETSC_TEST_LIB_COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE PETSC_TEST_LIB_OUTPUT)
  # Check program output
  if (PETSC_TEST_LIB_COMPILED AND PETSC_TEST_LIB_EXITCODE EQUAL 0)

    message(STATUS "Test PETSC_TEST_RUNS with shared library linking - Success")
    set(PETSC_TEST_RUNS TRUE)

    # Static libraries not required, so unset
    set_property(TARGET PETSC::petsc_static PROPERTY INTERFACE_LINK_LIBRARIES)

  else()

    message(STATUS "Test PETSC_TEST_RUNS with shared library linking - Failed")

    # Try to run test program (static linking)
    try_run(
      PETSC_TEST_LIB_EXITCODE
      PETSC_TEST_LIB_COMPILED
      ${CMAKE_CURRENT_BINARY_DIR}
      ${PETSC_TEST_LIB_CPP}
      CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
      LINK_LIBRARIES PETSC::petsc_static ${CMAKE_REQUIRED_LIBRARIES}
      COMPILE_OUTPUT_VARIABLE PETSC_TEST_LIB_COMPILE_OUTPUT
      RUN_OUTPUT_VARIABLE PETSC_TEST_LIB_OUTPUT)

    if (PETSC_TEST_LIB_COMPILED AND PETSC_TEST_LIB_EXITCODE EQUAL 0)
      message(STATUS "Test PETSC_TEST_RUNS static linking - Success")
      set(PETSC_TEST_RUNS TRUE)
    else()
      message(STATUS "Test PETSC_TEST_RUNS static linking - Failed")
      set(PETSC_TEST_RUNS FALSE)
    endif()

  endif()
endif()

# Check sizeof(PetscInt)
if (PETSC_INCLUDE_DIRS)
  include(CheckTypeSize)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${PETSC_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH})
  set(CMAKE_EXTRA_INCLUDE_FILES petscsys.h)
  check_type_size("PetscInt" PETSC_INT_SIZE LANGUAGE CXX)
  unset(CMAKE_EXTRA_INCLUDE_FILES)
  unset(CMAKE_REQUIRED_INCLUDES)
endif()

# Check compatibility of PetscScalar and STATE
if (PETSC_INCLUDE_DIRS)
  include(CheckTypeSize)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${PETSC_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH})
  set(CMAKE_EXTRA_INCLUDE_FILES petscsys.h)
  check_type_size("PetscScalar" PETSC_SCALAR_SIZE LANGUAGE CXX)
  set(PETSC_TEST_SCALAR_CPP
    "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/petsc_test_scalar.cpp")
  file(WRITE ${PETSC_TEST_SCALAR_CPP} "
  #include \"petscts.h\"
  #include \"petsc.h\"
  int main()
  {
    #if defined(PETSC_USE_COMPLEX)
      return 0;
      #else
      return 1;
      #endif
    return 0;
  }
  ")
  try_run(
    PETSC_TEST_SCALAR_EXITCODE
    PETSC_TEST_SCALAR_COMPILED
    ${CMAKE_CURRENT_BINARY_DIR}
    ${PETSC_TEST_SCALAR_CPP}
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}"
    LINK_LIBRARIES PETSC::petsc ${CMAKE_REQUIRED_LIBRARIES}
    COMPILE_OUTPUT_VARIABLE PETSC_TEST_SCALAR_COMPILE_OUTPUT
    RUN_OUTPUT_VARIABLE PETSC_TEST_SCALAR_OUTPUT)

  if (PETSC_TEST_SCALAR_COMPILED AND PETSC_TEST_SCALAR_EXITCODE EQUAL 0)
    ##PETSC_USE_COMPLEX IS DEFINED
    get_property(PZ_STATE_TYPE GLOBAL PROPERTY "PZ_STATE_TYPE")
    if(PZ_STATE_TYPE STREQUAL "float" AND PETSC_SCALAR_SIZE EQUAL 8)
    elseif(PZ_STATE_TYPE STREQUAL "double" AND PETSC_SCALAR_SIZE EQUAL 16)
    else()
      message(WARNING "PZ STATE is: ${PZ_STATE_TYPE}")
      message(WARNING "PetscScalar size is ${PETSC_SCALAR_SIZE}")
      message(FATAL_ERROR "PetscScalar and STATE type are incompatible")
    endif()
  elseif(PETSC_TEST_SCALAR_COMPILED AND PETSC_TEST_SCALAR_EXITCODE EQUAL 1)
    message(FATAL_ERROR "PETSc should be configured with complex scalars")
  endif()
endif()

unset(CMAKE_EXTRA_INCLUDE_FILES)
unset(CMAKE_REQUIRED_INCLUDES)

# Standard package handling
include(FindPackageHandleStandardArgs)
if (PETSC_FOUND)
  find_package_handle_standard_args(PETSc
    REQUIRED_VARS PETSC_FOUND PETSC_TEST_RUNS VERSION_VAR PETSC_VERSION
    FAIL_MESSAGE "PETSc could not be configured.")
else()
  find_package_handle_standard_args(PETSc
    REQUIRED_VARS PETSC_FOUND
    FAIL_MESSAGE "PETSc could not be found. Be sure to set PETSC_DIR.")
endif()
