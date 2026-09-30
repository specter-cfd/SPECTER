# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# mpiio: escritura y lectura de los campos con MPI-IO
add_library(mpiio_lib STATIC
  binary_mod.f90
  binary_io.f90
)

# Directorios de inclusión
target_include_directories(mpiio_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/mpiio
  ${CMAKE_BINARY_DIR}/src/mpiio
)

# Dependencias
add_dependencies(mpiio_lib fftp_lib timer_lib utils_lib)

# Bibliotecas
target_link_libraries(mpiio_lib PUBLIC project_mpi fftp_lib timer_lib utils_lib)
