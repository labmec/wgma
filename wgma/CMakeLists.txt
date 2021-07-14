add_library(wgma SHARED "")
target_compile_features(wgma PUBLIC cxx_std_17)


target_include_directories(wgma PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )

set(headers
  include/bctypes.hpp
  include/gmeshtools.hpp
  include/cmeshtools.hpp
  )

set(sources
  src/cmeshtools.cpp
  src/gmeshtools.cpp
  )


target_sources(wgma PUBLIC ${headers})
target_sources(wgma PRIVATE ${sources})

if(CMAKE_IS_PZ_BUILDTREE)
  target_link_libraries(wgma PUBLIC pz)
else()
  target_link_libraries(wgma PUBLIC NeoPZ::pz)
  target_include_directories(wgma PUBLIC ${PZ_INCLUDE_DIRS})
endif()