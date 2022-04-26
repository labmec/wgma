add_library(wgma SHARED "")
target_compile_features(wgma PUBLIC cxx_std_17)


target_include_directories(wgma PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )

set(headers
  include/modetypes.hpp
  include/bctypes.hpp
  include/pmltypes.hpp
  include/gmeshtools.hpp
  include/cmeshtools.hpp
  include/cmeshtools_impl.hpp
  include/wganalysis.hpp
  include/scattering.hpp
  include/slepcepshandler.hpp
  include/util.hpp
  )

set(sources
  src/cmeshtools.cpp
  src/gmeshtools.cpp
  src/wganalysis.cpp
  src/scattering.cpp
  src/slepcepshandler.cpp
  src/util.cpp
  )

target_sources(wgma PUBLIC ${headers})
target_sources(wgma PRIVATE ${sources})

if(CMAKE_IS_PZ_BUILDTREE)
  target_link_libraries(wgma PUBLIC pz)
else()
  target_link_libraries(wgma PUBLIC NeoPZ::pz)
  target_include_directories(wgma PUBLIC ${PZ_INCLUDE_DIRS})
endif()