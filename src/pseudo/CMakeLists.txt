# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# pseudo: grilla y rutinas pseudoespectrales (solo lo que usa HD)
add_library(pseudo_lib STATIC
  pseudospec_mod.f90
  module_dns.f90
  pseudospec_hd.f90
)

# Directorios de inclusión
target_include_directories(pseudo_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/pseudo
)

# Dependencias
add_dependencies(pseudo_lib fftp_lib)

# Bibliotecas
target_link_libraries(pseudo_lib PUBLIC project_mpi fftp_lib)
