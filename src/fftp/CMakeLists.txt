# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# fftp: FFT paralela con FC-Gram sobre FFTW
add_library(fftp_lib STATIC
  fftw_wrappers.h
  fprecision_mod.F90
  fftp_mod.F90
  fcgram_mod.f90
  fftp.F90
)

# La precisión, IKIND, CSIZE y NSTRIP se definen para todo el proyecto
# en el CMakeLists.txt de la raíz

# Directorios de inclusión
target_include_directories(fftp_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/fftp
  ${CMAKE_BINARY_DIR}/src/fftp
)

# Dependencias
add_dependencies(fftp_lib timer_lib)

# Bibliotecas
target_link_libraries(fftp_lib PUBLIC project_mpi project_fft timer_lib)
