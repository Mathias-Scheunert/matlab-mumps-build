#!/bin/bash
set -e

# Default values
: ${MUMPS_VERSION:=5.6.2}
: ${PREFIX:="${PWD}"}

# Download and unpack
mkdir -p src && cd src
if [ ! -f "MUMPS_${MUMPS_VERSION}.tar.gz" ]; then
    wget -c -O MUMPS_${MUMPS_VERSION}.tar.gz https://coin-or-tools.github.io/ThirdParty-Mumps/MUMPS_${MUMPS_VERSION}.tar.gz
fi
if [ ! -d "MUMPS_${MUMPS_VERSION}" ]; then
    tar -xzf MUMPS_${MUMPS_VERSION}.tar.gz
fi
cd MUMPS_${MUMPS_VERSION}

# Build MUMPS static lib
cp ../../Makefile_mumps.inc Makefile.inc
# -k || true ignores all errors
make d openmp=1 -k || true
make z openmp=1 -k || true

# Move static libs abd header files to desired folders
cp lib/*.a "${PREFIX}/lib/"
cp include/*.h "${PREFIX}/include/"

# Test
cd examples
echo " "
echo "running tests ..."
echo " "
cp ../../../check_test .
./dsimpletest < input_simpletest_real | ./check_test
./zsimpletest < input_simpletest_cmplx | ./check_test
cd ..

