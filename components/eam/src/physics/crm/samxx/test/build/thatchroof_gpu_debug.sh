#!/bin/bash

unset YAKL_ARCH
unset NCRMS

export NCHOME=${NETCDF_PATH}
export NFHOME=/usr
export NCRMS=42
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export FFLAGS="-O3 -ffree-line-length-none"
export YAKL_HOME="`pwd`/../../../../../../../../externals/YAKL"
export YAKL_ARCH="CUDA"
export YAKL_CUDA_FLAGS="-arch sm_35 -O0 -g -G -DYAKL_DEBUG -DTHRUST_IGNORE_CUB_VERSION_CHECK"


