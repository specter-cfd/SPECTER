#!/bin/bash
# Builds SPECTER with OpenMP offload for NVIDIA GPUs.
# Usage: build_nvidia.sh [llvm|nvhpc] [ARCH] [SINGLE|DOUBLE]
#   llvm : LLVM flang with NVPTX support and a CUDA-aware MPI library.
#   nvhpc: nvfortran with its own MPI; needs compute capability >=7.0
# P_GPU implies P_HYBRID and selects the fftp-gpu backend (cuFFT).
# The binaries go to SPECTER/bin/nvidia; the build tree to
# build/nvidia-<llvm|nvhpc>.
set -e
TC=${1:-llvm}; PREC=${3:-DOUBLE}
SRC=$(cd $(dirname $0)/.. && pwd)
BLD=$SRC/build/nvidia-$TC
mkdir -p $BLD && cd $BLD
if [ $TC = llvm ]; then
  ARCH=${2:-sm_61}
  module purge; module load cmake/4.1.2 cuda/12.4 hwloc/2.12.0
  module load llvm-nvptx/23.1.0; module load openmpi5/5.0.7
  CUDA_HOME=/opt/ohpc/pub/utils/cuda/12.4
  CC=clang FC=mpifort cmake $SRC \
    -DFFTW3_ROOT=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11 \
    -DP_GPU=ON -DGPU_VENDOR=NVIDIA -DGPU_ARCH=$ARCH -DCUDA_PATH=$CUDA_HOME \
    -DPRECISION=$PREC -DBIN_DIR=$SRC/bin/nvidia
else
  ARCH=${2:-cc70}
  # NVHPC compilers older than 24.11 have a bug in use_device_addr
  # and don't offload properly, use LLVM compilers instead if you
  # need support for older GPUs or CUDA versions.
  module purge; module load nvhpc-stack/24.5; module load nvhpc/24.5 cmake/4.1.2
  NV=/opt/ohpc/pub/compiler/nvhpc/24.5/Linux_x86_64/24.5
  CC=nvc FC=$NV/comm_libs/12.4/hpcx/hpcx-2.19/ompi/bin/mpifort cmake $SRC \
    -DFFTW3_ROOT=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11 \
    -DP_GPU=ON -DGPU_VENDOR=NVIDIA -DGPU_ARCH=$ARCH -DCUDA_PATH=$NV/cuda/12.4 \
    -DPRECISION=$PREC -DBIN_DIR=$SRC/bin/nvidia
fi
make -j 16
