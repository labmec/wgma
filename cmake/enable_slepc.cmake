function(enable_slepc target)
  include(cmake/enable_petsc.cmake)
  enable_petsc(${target})
  if(NOT DEFINED SLEPC_DIR)
	  if(DEFINED ENV{SLEPC_DIR})
		  set(SLEPC_DIR $ENV{SLEPC_DIR})
	  else()
		  #########################################################
		  #				The	value of SLEPC_DIR is estabilished by       #
		  # 				your installation directory.				          #
		  #########################################################
		  message(FATAL_ERROR "SLEPC_DIR environment variable is not defined.")
	  endif()
  endif()
  find_package(SLEPc REQUIRED)
  get_target_property(SLEPC_INCLUDE_DIRECTORIES SLEPC::slepc INTERFACE_INCLUDE_DIRECTORIES)
  target_include_directories(${target} PRIVATE ${SLEPC_INCLUDE_DIRECTORIES})
  target_compile_definitions(${target} PUBLIC WGMA_USING_SLEPC)
  target_link_libraries(${target} PRIVATE SLEPC::slepc)
endfunction()