function(enable_petsc target)
  if(NOT DEFINED PETSC_DIR)
		if(DEFINED ENV{PETSC_DIR})
			set(PETSC_DIR $ENV{PETSC_DIR})
		else()
			#########################################################
		  #				The	value of PETSC_DIR is estabilished by       #
		  # 				your installation directory.				          #
		  #########################################################
		  message(FATAL_ERROR "SLEPC_DIR environment variable is not defined.")
		endif()
	endif()

	if(NOT DEFINED PETSC_ARCH)
		if(DEFINED ENV{PETSC_ARCH})
			set(PETSC_ARCH $ENV{PETSC_ARCH})
		endif()
	endif()

	find_package(PETSc REQUIRED)
	get_target_property(PETSC_INCLUDE_DIRECTORIES PETSC::petsc INTERFACE_INCLUDE_DIRECTORIES)
  target_include_directories(${target} PRIVATE ${PETSC_INCLUDE_DIRECTORIES})
  target_link_libraries(${target} PRIVATE PETSC::petsc MPI::MPI_CXX)
  set(PETSC_ARCH ${PETSC_ARCH} PARENT_SCOPE)
endfunction()