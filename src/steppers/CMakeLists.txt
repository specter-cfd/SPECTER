# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# steppers: integradores temporales (por ahora solo el de Canuto)
add_library(stepper_lib STATIC
  gstepper_base.f90
  canuto_stepper.f90
)

# Directorios de inclusión
target_include_directories(stepper_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/steppers
  ${CMAKE_BINARY_DIR}/src/steppers
)

# Dependencias
add_dependencies(stepper_lib pdes_lib)

# Bibliotecas (pdes_lib trae memmgt, pseudo, fftp y los directorios de sus módulos)
target_link_libraries(stepper_lib PUBLIC pdes_lib)
