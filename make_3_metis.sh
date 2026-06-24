#!/bin/bash
set -e

# Default values
: ${METIS_VERSION:=5.1.0}
: ${PREFIX:="${PWD}"}

# Download and unpack
mkdir -p src && cd src
wget -c -O metis-${METIS_VERSION}.tar.gz https://sourceforge.net/projects/openfoam-extend/files/foam-extend-3.0/ThirdParty/metis-${METIS_VERSION}.tar.gz/download
tar -xzf metis-${METIS_VERSION}.tar.gz
cd metis-${METIS_VERSION}

# Ensure that correct CMake-version is used
sed -i '1s/.*/cmake_minimum_required(VERSION 3.5)/' CMakeLists.txt

# Modify metis header file from 32 to 64 bit
sed -i 's/#define IDXTYPEWIDTH 32/#define IDXTYPEWIDTH 64/' include/metis.h

# Build METIS static lib
# r64 float   -> =1 for double (incompatible with metis-mumps interface!)
# i64 integer -> =1 for 64 bit 
make config prefix=${PREFIX} openmp=1 i64=1 r64=0
make

# Install i.e. copy content to desired paths
make install
