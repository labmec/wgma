add_library(wgma SHARED "")
target_compile_features(wgma PUBLIC cxx_std_17)


target_include_directories(wgma PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )

set(headers
  include/modetypes.hpp
  include/bctypes.hpp
  include/pmltypes.hpp
  include/cartesian_pml.hpp
  include/cylindrical_pml.hpp
  include/gmeshtools.hpp
  include/cmeshtools.hpp
  include/cmeshtools_impl.hpp
  include/wganalysis.hpp
  include/scattering.hpp
  include/slepcepshandler.hpp
  include/twisted_wgma.hpp
  include/precond.hpp
  include/util.hpp
  include/post/eldata.hpp
  include/post/integrator.hpp
  include/post/solutionnorm.hpp
  include/post/orthosol.hpp
  )

set(sources
  src/cmeshtools.cpp
  src/gmeshtools.cpp
  src/cartesian_pml.cpp
  src/cylindrical_pml.cpp
  src/wganalysis.cpp
  src/scattering.cpp
  src/slepcepshandler.cpp
  src/twisted_wgma.cpp
  src/precond.cpp
  src/util.cpp
  src/post/integrator.cpp
  src/post/solutionnorm.cpp
  src/post/orthosol.cpp
  )

target_sources(wgma PUBLIC ${headers})
target_sources(wgma PRIVATE ${sources})

if(CMAKE_IS_PZ_BUILDTREE)
  target_link_libraries(wgma PUBLIC pz)
else()
  target_link_libraries(wgma PUBLIC NeoPZ::pz)
  target_include_directories(wgma PUBLIC ${PZ_INCLUDE_DIRS})
endif()