add_library(pdes_lib STATIC
  status_mod.f90
  equation_base.f90
  hd_solver.f90
  equation_factory.f90
)

target_include_directories(pdes_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/pdes
  ${CMAKE_BINARY_DIR}/src/pdes
  ${CMAKE_BINARY_DIR}/src/memmgt
)

# Make sure dependencies are satisfied
add_dependencies(pdes_lib fftp_lib pseudo_lib mpiio_lib memmgt_lib boundary_lib)

# Link libraries
target_link_libraries(pdes_lib PUBLIC project_mpi fftp_lib pseudo_lib mpiio_lib memmgt_lib boundary_lib)
