add_pz_target(
  NAME stepfiber
  SOURCES stepfiber.cpp
  REQUIRED PZ_USING_LAPACK)

target_link_libraries(stepfiber PUBLIC wgma)