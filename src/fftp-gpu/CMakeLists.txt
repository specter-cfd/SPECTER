# FFTP-GPU: parallel FFT with cuFFT/hipFFT on the device and FFTW on
# the host (requires P_GPU)
add_library(fftp_lib STATIC
  fftw_wrappers.h
  gfft_defs.h
  fprecision_mod.F90
  fftp_mod.F90
  gfft_bindings.F90
  fcgram_mod.f90
  fftp.F90
)

# Include directories
target_include_directories(fftp_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/fftp-gpu
  ${CMAKE_BINARY_DIR}/src/fftp-gpu
  ${CMAKE_BINARY_DIR}/src/memmgt
)

# Make sure dependencies are satisfied
add_dependencies(fftp_lib timer_lib gdevice_lib)

# Link MPI, the FFT libraries and the device management
target_link_libraries(fftp_lib PUBLIC project_mpi project_fft gdevice_lib timer_lib)
