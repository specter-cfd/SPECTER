add_library(ic_lib STATIC
  ic_base.f90
  ic_velocity.f90
  ic_factory.f90
)

target_include_directories(ic_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/ic
  ${CMAKE_BINARY_DIR}/src/ic
)

# Make sure dependencies are satisfied
add_dependencies(ic_lib pdes_lib utils_lib)

# Link libraries
target_link_libraries(ic_lib PUBLIC pdes_lib utils_lib)
