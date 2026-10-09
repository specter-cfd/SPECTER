# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# pdes: la lectura de &status (status_mod.f90), la clase base de los
# solvers (equation_base.f90) y el solver HD (hd_solver.f90)
add_library(pdes_lib STATIC
  status_mod.f90
  equation_base.f90
  hd_solver.f90
)

# Directorios de inclusión. memmgt_lib no publica la carpeta de sus
# módulos, así que se agrega acá, como en GHOST.
target_include_directories(pdes_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/pdes
  ${CMAKE_BINARY_DIR}/src/pdes
  ${CMAKE_BINARY_DIR}/src/memmgt
)

# Dependencias
add_dependencies(pdes_lib fftp_lib pseudo_lib mpiio_lib memmgt_lib boundary_lib)

# Bibliotecas
target_link_libraries(pdes_lib PUBLIC project_mpi fftp_lib pseudo_lib mpiio_lib memmgt_lib boundary_lib)
