# pseudobase: the basic modules (grid, wavenumbers, box, fft plans,
# random numbers), needed by the memory management before the
# pseudospectral kernels can be built
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

# pseudo: the pseudospectral kernels; they take their temporaries
# from the workspace pool (memmgt_lib)
add_library(pseudo_lib STATIC
  pseudospec_hd.f90
)

target_include_directories(pseudo_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/pseudo
  ${CMAKE_BINARY_DIR}/src/memmgt
)

# Make sure dependencies are satisfied
add_dependencies(pseudo_lib fftp_lib pseudobase_lib memmgt_lib)

# Link libraries
target_link_libraries(pseudo_lib PUBLIC project_mpi fftp_lib pseudobase_lib memmgt_lib)
