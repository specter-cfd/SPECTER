# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# pseudobase: grilla, números de onda, caja y plan del FC-Gram, que
# necesita la gestión de memoria antes de poder armar los núcleos
# pseudoespectrales
add_library(pseudobase_lib STATIC
  pseudospec_mod.f90
)

target_include_directories(pseudobase_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/memmgt
)

add_dependencies(pseudobase_lib fftp_lib gmem_lib)
target_link_libraries(pseudobase_lib PUBLIC project_mpi fftp_lib gmem_lib)

# pseudo: los núcleos pseudoespectrales de HD; toman sus temporarios
# del pool del workspace (memmgt_lib)
add_library(pseudo_lib STATIC
  pseudospec_hd.f90
)

target_include_directories(pseudo_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/memmgt
)

# Dependencias
add_dependencies(pseudo_lib fftp_lib pseudobase_lib memmgt_lib)

# Bibliotecas
target_link_libraries(pseudo_lib PUBLIC project_mpi fftp_lib pseudobase_lib memmgt_lib)
