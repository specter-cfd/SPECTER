# boundary: boundary conditions of the velocity (boundary_mod.f90
# includes vboundary.f90)
add_library(boundary_lib STATIC boundary_mod.f90)

target_include_directories(boundary_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/boundary
  ${CMAKE_BINARY_DIR}/src/boundary
  ${CMAKE_BINARY_DIR}/src/memmgt
)

# Make sure dependencies are satisfied
add_dependencies(boundary_lib fftp_lib pseudo_lib memmgt_lib)

# Link libraries
target_link_libraries(boundary_lib PUBLIC project_mpi fftp_lib pseudo_lib memmgt_lib)
