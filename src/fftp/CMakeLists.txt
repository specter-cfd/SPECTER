add_library(fftp_lib STATIC
  fftw_wrappers.h
  fprecision_mod.F90
  fftp_mod.F90
  fcgram_mod.f90
  fftp.F90
)

# Include directories
target_include_directories(fftp_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/fftp
  ${CMAKE_BINARY_DIR}/src/fftp
)

# Make sure dependencies are satisfied
add_dependencies(fftp_lib timer_lib)

# Link libraries
target_link_libraries(fftp_lib PUBLIC project_mpi project_fft timer_lib)
