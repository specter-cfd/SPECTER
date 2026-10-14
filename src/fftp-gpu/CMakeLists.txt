# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# fftp-gpu: FFT paralela con FC-Gram, con cuFFT/hipFFT en el dispositivo
# y FFTW en el host (requiere P_GPU)
add_library(fftp_lib STATIC
  fftw_wrappers.h
  gfft_defs.h
  fprecision_mod.F90
  fftp_mod.F90
  gfft_bindings.F90
  fcgram_mod.f90
  fftp.F90
)

# La precisión, IKIND, CSIZE y NSTRIP se definen para todo el proyecto
# en el CMakeLists.txt de la raíz

# Directorios de inclusión (gdevice.mod queda en build/src/memmgt)
target_include_directories(fftp_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/fftp-gpu
  ${CMAKE_BINARY_DIR}/src/fftp-gpu
  ${CMAKE_BINARY_DIR}/src/memmgt
)

# Dependencias
add_dependencies(fftp_lib timer_lib gdevice_lib)

# Bibliotecas
target_link_libraries(fftp_lib PUBLIC project_mpi project_fft gdevice_lib timer_lib)

