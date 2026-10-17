# -------------------------------------------------------------------
# SPECTER UserConfig CMake file
# Cambios respecto del de GHOST: ver CMAKE_CHANGES_FROM_GHOST.txt
#
# Cada valor se puede cambiar al configurar, sin editar este archivo,
# por ejemplo desde build/: "cmake .. -DPRECISION=SINGLE"
# Si las bibliotecas no están en el path, usar -DCMAKE_PREFIX_PATH
# y -DFFTW3_ROOT, por ejemplo
# "cmake .. -DCMAKE_PREFIX_PATH=/opt/openmpi-5.0.8 -DFFTW3_ROOT=/opt/fftw-3.3.10"
# -------------------------------------------------------------------

# ----------------------- Opciones / variantes (GHOST) ---------------
# Las de GHOST. Tienen efecto PRECISION y las de GPU (P_GPU, GPU_VENDOR
# y GPU_ARCH; probado con amdflang en las AMD MI210, gfx90a, y con el
# flang de LLVM 23 en las NVIDIA Tesla P4, sm_61). P_HYBRID
# está bloqueada: queda en OFF salvo con P_GPU, que la enciende.
set(P_HYBRID   OFF      CACHE BOOL   "Enable hybrid OpenMP (bloqueada en SPECTER)" FORCE)
option(P_GPU    "Enable OpenMP offload to GPUs (implies P_HYBRID)" OFF)
set(GPU_VENDOR "NVIDIA" CACHE STRING "GPU vendor for P_GPU (NVIDIA or AMD)")
set(GPU_ARCH   ""       CACHE STRING "GPU architecture for P_GPU (e.g. cc80, sm_80, gfx90a)")
set(PRECISION  "DOUBLE" CACHE STRING "Precision (SINGLE or DOUBLE)")
set(FFTP       "fftp"   CACHE STRING "FFT library (fftp)")
set(IOLIB      "mpiio"  CACHE STRING "IO library (mpiio)")

# ----------------------- Advanced set up params ---------------------
# IKIND: 4 en máquinas de 32 bits, 8 en máquinas de 64 bits
# CSIZE: 8 si la cache L1 es <= 64 kb, 16 si es de 128 kb
# NSTRIP: strip mining (en general 1)
set(IKIND  8  CACHE STRING "Integer kind for pointers")
set(CSIZE  32 CACHE STRING "Cache size")
set(NSTRIP 1  CACHE STRING "Strip mining")

# ----------------------- Solver (SPECTER) ---------------------------
# Solo HD está soportado
set(SOLVER "HD" CACHE STRING "Solver (HD)")
