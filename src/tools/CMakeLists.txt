# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# tools: BOOTS, el regrillador (sin la biblioteca para Python de GHOST)
add_executable(BOOTS "${SRC_DIR}/tools/boots3D.f90")

# En SPECTER pseudo_lib no llega por mpiio_lib, como en GHOST, así que
# se enlaza aparte; utils_lib llega por mpiio_lib. Each archive must
# appear once: with amdflang a duplicate makes the device link fail
# with duplicate symbols.
target_link_libraries(BOOTS
    PRIVATE
        project_includes
        project_gpu
        mpiio_lib
        pseudo_lib
)

set_target_properties(BOOTS PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
)
