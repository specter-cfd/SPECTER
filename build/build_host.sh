#!/bin/bash
# Builds SPECTER for the host with gfortran (pure MPI: hybrid MPI+OpenMP
# without GPUs is still blocked in UserConfig.cmake).
# Usage: build_host.sh [SINGLE|DOUBLE]
# The binaries go to SPECTER/bin (BIN_DIR); the build tree to
# build/host. Run from any directory.
set -e
PREC=${1:-DOUBLE}
SRC=$(cd $(dirname $0)/.. && pwd)
BLD=$SRC/build/host
module purge; module load gnu15 openmpi5 cmake/4.1.2 fftw/3.3.11
mkdir -p $BLD && cd $BLD
CC=gcc FC=gfortran cmake $SRC \
  -DFFTW3_ROOT=/opt/ohpc/pub/libs/gnu15/fftw/3.3.11 \
  -DP_GPU=OFF -DPRECISION=$PREC -DBIN_DIR=$SRC/bin
make -j 16
