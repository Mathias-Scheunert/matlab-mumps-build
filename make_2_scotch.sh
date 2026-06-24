#!/bin/bash
set -e

# Default values
: ${SCOTCH_VERSION:=v7.0.3}
: ${PREFIX:="${PWD}"}

# Download and unpack
mkdir -p src && cd src
git clone https://gitlab.inria.fr/scotch/scotch.git
cd scotch
git checkout tags/${SCOTCH_VERSION}
cd src

# Build SCOTCH static lib (DEBIAN system!)
# -> the linux packages bison and flex are required! (UBUNTU: sudo apt-get install bison flex)
cp ./Make.inc/Makefile.inc.x86-64_pc_linux2 Makefile.inc

# Add additional flags to ensure 64bit Integer and avoid error instead of warning in case of pointer issues (may occur with >gcc14)
sed -i '/^CFLAGS[[:space:]]*=/ s/$/ -DINTSIZE64 -std=gnu99 -Wno-error=incompatible-pointer-types/' Makefile.inc

# Add further flags to tell Fortran-Compiler that standard Integer have 8 Byte
make esmumps \
     FFLAGS="-O3 -fdefault-integer-8"

# Install i.e. copy content to desired paths
make install prefix=${PREFIX}
