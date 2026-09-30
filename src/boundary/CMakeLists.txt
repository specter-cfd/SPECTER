# boundary: condiciones de borde. boundary_mod.fpp incluye adentro
# vboundary.f90 (y sboundary/bboundary en otros solvers), que por eso
# no se compilan por separado.
add_library(boundary_lib STATIC boundary_mod.fpp)

# boundary_mod.fpp pasa por el preprocesador, en formato libre
set_source_files_properties(boundary_mod.fpp
  PROPERTIES Fortran_PREPROCESS ON Fortran_FORMAT FREE)

# El solver, que el Makefile pasaba a cpp en "make edit"
target_compile_definitions(boundary_lib PRIVATE ${SOLVER}_SOL)

# Directorios de inclusión
target_include_directories(boundary_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/boundary
  ${CMAKE_BINARY_DIR}/src/boundary
)

# Dependencias
add_dependencies(boundary_lib fftp_lib pseudo_lib)

# Bibliotecas
target_link_libraries(boundary_lib PUBLIC project_mpi fftp_lib pseudo_lib)
