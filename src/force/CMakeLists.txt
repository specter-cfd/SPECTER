# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# force: forzados (por ahora, solo los de la velocidad)
add_library(force_lib STATIC
  force_base.f90
  force_velocity.f90
  force_factory.f90
)

# Directorios de inclusión
target_include_directories(force_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/force
  ${CMAKE_BINARY_DIR}/src/force
)

# Dependencias
add_dependencies(force_lib pdes_lib utils_lib)

# Bibliotecas (pdes_lib trae memmgt, pseudo, mpiio, fftp y los directorios de sus módulos)
target_link_libraries(force_lib PUBLIC pdes_lib utils_lib)
