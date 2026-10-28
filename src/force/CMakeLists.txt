add_library(force_lib STATIC
  force_base.f90
  force_velocity.f90
  force_factory.f90
)

target_include_directories(force_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/force
  ${CMAKE_BINARY_DIR}/src/force
)

# Make sure dependencies are satisfied
add_dependencies(force_lib pdes_lib utils_lib)

# Link libraries
target_link_libraries(force_lib PUBLIC pdes_lib utils_lib)
