add_pz_target(
  NAME wr90
  SOURCES wr90.cpp
  REQUIRED PZ_USING_LAPACK)

target_link_libraries(wr90 PUBLIC wgma)