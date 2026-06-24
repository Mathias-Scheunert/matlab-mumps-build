#!/bin/bash
set -e

# Default values
: ${MUMPS_VERSION:=5.6.2}
: ${PREFIX:="${PWD}"}

# Get path to matlab / mex
MATLAB_ROOT=$(dirname $(dirname $(readlink -f $(which matlab))))
export MATLAB_ROOT="$MATLAB_ROOT"
MATLAB_MEX="$MATLAB_ROOT/bin/mex"

# Build mex-files
cd src/MUMPS_${MUMPS_VERSION}/MATLAB
cp Makefile Makefile.bak
cp ../../../Makefile_matlab ./Makefile
# Replace inconsitant / outdated c file
cp mumpsmex.c mumpsmex.c.bak
cp ../../../mumpsmex_64.c mumpsmex.c
#make MEX="${MATLAB_MEX} -largeArrayDims"
make

# Install
mkdir -p "${PREFIX}/lib/matlab"
cp *.m *.mex* *.mat "${PREFIX}/lib/matlab"

# Test MATLAB
echo " "
echo "running tests ..."
echo " "
cd /tmp
#export LD_LIBRARY_PATH=$MATLAB_ROOT/bin/glnxa64:$LD_LIBRARY_PATH
MATLABPATH="${PREFIX}/lib/matlab:${MATLABPATH}" matlab -nojvm -batch 'simple_example; zsimple_example; schur_example; diagainv_example; multiplerhs_example; sparserhs_example; polyfit(1:10, sin(1:10), 2)'

# checks:
# ldd ${PREFIX}/lib/matlab/dmumpsmex.mexa64
# -> list of dynamic linked libs *so must not contain metis, scotch, pord
# no "...=> not found" should be visible
