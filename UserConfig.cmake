# -------------------------------------------------------------------
# SPECTER UserConfig CMake file
#
# If libraries are not in the path, use -DCMAKE_PREFIX_PATH
# and -DFFTW3_ROOT to set paths, e.g., from the build directory do
# "cmake .. -DCMAKE_PREFIX_PATH=/opt/openmpi-5.0.8 -DFFTW3_ROOT=/opt/fftw-3.3.10"
# -------------------------------------------------------------------

# ----------------------- Options / variants ------------------------
# P_HYBRID is locked to OFF (P_GPU turns it on)
set(P_HYBRID   OFF      CACHE BOOL   "Enable hybrid OpenMP (locked in SPECTER)" FORCE)
option(P_GPU    "Enable OpenMP offload to GPUs (implies P_HYBRID)" OFF)
set(GPU_VENDOR "NVIDIA" CACHE STRING "GPU vendor for P_GPU (NVIDIA or AMD)")
set(GPU_ARCH   ""       CACHE STRING "GPU architecture for P_GPU (e.g. cc80, sm_80, gfx90a)")
set(PRECISION  "DOUBLE" CACHE STRING "Precision (SINGLE or DOUBLE)")
set(FFTP       "fftp"   CACHE STRING "FFT library (fftp)")
set(IOLIB      "mpiio"  CACHE STRING "IO library (mpiio)")

# ----------------------- Advanced set up params --------------------
# IKIND    4 in 32 bits machines, 8 in 64 bits machines
# CSIZE    8 if L1 cache is <= 64 kb, 16 if L1 cache is 128 kb
# NSTRIP   controls strip mining (often set to 1)
set(IKIND  8  CACHE STRING "Integer kind for pointers")
set(CSIZE  32 CACHE STRING "Cache size")
set(NSTRIP 1  CACHE STRING "Strip mining")

# ----------------------- Solver ------------------------------------
set(SOLVER "HD" CACHE STRING "Solver (HD)")
