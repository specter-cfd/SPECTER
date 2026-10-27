add_library(stepper_lib STATIC
  gstepper_base.f90
  canuto_stepper.f90
  stepper_factory.f90
)

target_include_directories(stepper_lib PUBLIC
  ${CMAKE_SOURCE_DIR}/src/steppers
  ${CMAKE_BINARY_DIR}/src/steppers
)

# Make sure dependencies are satisfied
add_dependencies(stepper_lib pdes_lib)

# Link libraries
target_link_libraries(stepper_lib PUBLIC pdes_lib)
