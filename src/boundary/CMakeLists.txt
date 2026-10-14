# Propio de SPECTER (no existe en GHOST): ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# boundary: condiciones de borde de la velocidad. boundary_mod.f90
# incluye adentro vboundary.f90, que por eso no se compila por separado.
add_library(boundary_lib STATIC boundary_mod.f90)

# Directorios de inclusión
target_include_directories(boundary_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/boundary
  ${CMAKE_BINARY_DIR}/src/boundary
)

# Dependencias
add_dependencies(boundary_lib fftp_lib pseudo_lib)

# Bibliotecas
target_link_libraries(boundary_lib PUBLIC project_mpi fftp_lib pseudo_lib)
