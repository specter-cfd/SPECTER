# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# pseudo: grilla y rutinas pseudoespectrales (solo lo que usa HD)
add_library(pseudo_lib STATIC
  pseudospec_mod.fpp
  module_dns.f90
  pseudospec_hd.f90
)

# pseudospec_mod.fpp pasa por el preprocesador, en formato libre
set_source_files_properties(pseudospec_mod.fpp
  PROPERTIES Fortran_PREPROCESS ON Fortran_FORMAT FREE)

# Grilla y FC-Gram, que el Makefile viejo pasaba a cpp en "make edit"
target_compile_definitions(pseudo_lib PRIVATE
  NX_=${NX} NY_=${NY} NZ_=${NZ}
  CX_=${CX} CY_=${CY} CZ_=${CZ}
  OX_=${OX} OY_=${OY} OZ_=${OZ})

# Directorios de inclusión
target_include_directories(pseudo_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/pseudo
)

# Dependencias
add_dependencies(pseudo_lib fftp_lib)

# Bibliotecas
target_link_libraries(pseudo_lib PUBLIC project_mpi fftp_lib)
