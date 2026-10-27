# ---------------------- SPECTER tools --------------------------------
# Boots3D:
add_executable(BOOTS "${SRC_DIR}/tools/boots3D.f90")

# utils_lib comes in through mpiio_lib
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
