# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# pdes: por ahora solo la lectura de &status (status_mod.f90)
add_library(pdes_lib STATIC
  status_mod.f90
)

# Directorios de inclusión
target_include_directories(pdes_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/pdes
  ${CMAKE_BINARY_DIR}/src/pdes
)

# Dependencias
add_dependencies(pdes_lib fftp_lib)

# Bibliotecas
target_link_libraries(pdes_lib PUBLIC project_mpi fftp_lib)

