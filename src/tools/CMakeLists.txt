# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt en la raíz
# tools: BOOTS, el regrillador (sin la biblioteca para Python de GHOST)
add_executable(BOOTS "${SRC_DIR}/tools/boots3D.f90")

# En SPECTER pseudo_lib no llega por utils_lib ni por mpiio_lib, como
# en GHOST, así que se enlaza aparte (y con ella fftp_lib)
target_link_libraries(BOOTS
    PRIVATE
        project_includes
        project_gpu
        utils_lib
        mpiio_lib
        pseudo_lib
)

set_target_properties(BOOTS PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
)
