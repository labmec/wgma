add_pz_target(
  NAME wr90
  SOURCES wr90.cpp)
target_link_libraries(wr90 PUBLIC wgma)

add_pz_target(
  NAME stepfiber
  SOURCES stepfiber.cpp)
target_link_libraries(stepfiber PUBLIC wgma)

add_pz_target(
  NAME ribwg
  SOURCES ribwg.cpp
  FILES meshes/ribWG.msh)
target_link_libraries(ribwg PUBLIC wgma)

add_pz_target(
  NAME ecf
  SOURCES ecf.cpp
  FILES meshes/ecf.msh meshes/ecf_closed.msh)
target_link_libraries(ecf PUBLIC wgma)

add_pz_target(
  NAME planar_wg
  SOURCES planar_wg.cpp 
  FILES meshes/planar_wg.msh meshes/wg_disc.msh)
target_link_libraries(planar_wg PUBLIC wgma)

add_pz_target(
  NAME pcwg
  SOURCES pcwg.cpp
  FILES meshes/pcwg.msh meshes/pcwg_circdata.csv
  FILES meshes/pcwg_test.msh meshes/pcwg_test_circdata.csv
  )
target_link_libraries(pcwg PUBLIC wgma)


add_pz_target(
  NAME planar_wg_periodic
  SOURCES planar_wg_periodic.cpp
  FILES meshes/periodic_test.msh
  )
target_link_libraries(planar_wg_periodic PUBLIC wgma)
