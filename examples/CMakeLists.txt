add_pz_target(
  NAME wr90
  SOURCES wr90.cpp)
target_link_libraries(wr90 PUBLIC wgma)

add_pz_target(
  NAME stepfiber
  SOURCES stepfiber.cpp
  FILES meshes/stepfiber.msh meshes/stepfiber_circdata.csv
  )
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


list(APPEND slab_disc_files ../py/scripts/slab_disc_validation.py ../py/scripts/slab_disc_nmodes.py)
foreach(file ${slab_disc_files})
  get_filename_component(filename ${file} NAME)
  configure_file(${file} scripts/${filename})
endforeach(file)
add_pz_target(
  NAME slab_disc
  SOURCES slab_disc.cpp
  FILES meshes/slab_disc.msh meshes/slab_disc_validation.msh
  )
target_link_libraries(slab_disc PUBLIC wgma)

add_pz_target(
  NAME sf3d
  SOURCES sf3d.cpp
  FILES meshes/sf3d.msh meshes/sf3d_cyldata.csv
  )
  target_link_libraries(sf3d PUBLIC wgma)

add_pz_target(
  NAME sf3d_cyl
  SOURCES sf3d_cyl.cpp
  )
target_link_libraries(sf3d_cyl PUBLIC wgma)  

add_pz_target(
  NAME pcf
  SOURCES pcf.cpp
  FILES meshes/pcf.msh meshes/pcf_circdata.csv
  )
target_link_libraries(pcf PUBLIC wgma)  


add_pz_target(
  NAME aniso_pcf
  SOURCES aniso_pcf.cpp
  FILES meshes/pcf.msh meshes/pcf_circdata.csv
)
target_link_libraries(aniso_pcf PUBLIC wgma)

add_pz_target(
  NAME acoustic_pcf
  SOURCES acoustic_pcf.cpp
  FILES meshes/pcf.msh meshes/pcf_circdata.csv
)
target_link_libraries(acoustic_pcf PUBLIC wgma)

add_pz_target(
  NAME acoustic_pcf2
  SOURCES acoustic_pcf2.cpp
  FILES meshes/pcf.msh meshes/pcf_circdata.csv
)
target_link_libraries(acoustic_pcf2 PUBLIC wgma)



file(GLOB meta_surf_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "../py/scripts/meta_surf*json")
list(APPEND meta_surf_files ../py/scripts/meta_surf_1d_freq_sweep.py)
list(APPEND meta_surf_files ../py/scripts/meta_surf_1d_indices.py)
set(meta_surf_data ../py/scripts/data/AZ.dat ../py/scripts/data/Cu1.dat)

foreach(file ${meta_surf_files})
  get_filename_component(filename ${file} NAME)
  configure_file(${file} scripts/${filename})
endforeach(file)
foreach(file ${meta_surf_data})
  get_filename_component(filename ${file} NAME)
  configure_file(${file} scripts/data/${filename})
endforeach(file)

add_pz_target(
  NAME meta_surf_1d
  SOURCES meta_surf_1d.cpp
  FILES meshes/meta_surf_1d.msh
)
target_link_libraries(meta_surf_1d PUBLIC wgma)