# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# ic: condiciones iniciales (por ahora, solo las de la velocidad)
add_library(ic_lib STATIC
  ic_base.f90
  ic_velocity.f90
  ic_factory.f90
)

# Directorios de inclusión
target_include_directories(ic_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/ic
  ${CMAKE_BINARY_DIR}/src/ic
)

# Dependencias
add_dependencies(ic_lib pdes_lib utils_lib)

# Bibliotecas (pdes_lib trae memmgt, pseudo, mpiio, fftp y los directorios de sus módulos)
target_link_libraries(ic_lib PUBLIC pdes_lib utils_lib)
