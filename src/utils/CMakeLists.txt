# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# gtimer: sin dependencias dentro de SPECTER
add_library(timer_lib STATIC gtimer_mod.f90)
# gutils: usa fprecision, commtypes, mpivars y threads (de fftp) y gtimer
add_library(utils_lib STATIC gutils_mod.f90)

# Directorios de inclusión
target_include_directories(timer_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/utils
  ${CMAKE_BINARY_DIR}/src/utils
)
target_include_directories(utils_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/utils
  ${CMAKE_BINARY_DIR}/src/utils
  ${CMAKE_BINARY_DIR}/src/${FFTP}
)

# Dependencias
add_dependencies(utils_lib timer_lib fftp_lib)

# Bibliotecas
target_link_libraries(timer_lib PUBLIC project_mpi)
target_link_libraries(utils_lib PUBLIC project_mpi timer_lib fftp_lib)
