Forked from: https://github.com/xmjiao/mumps4m-openmp and https://github.com/blechta/mumps-matlab-recipes \
Original Authors: Xiangmin Jiao, Jan Blechta \
Code in this repository was generated and refined with the assistance of Google Gemini and OpenAI ChatGPT.

Builds 32-bit MUMPS and statically links dependent libraries \
-> to get the 64-bit variant please checkout branch **mumps-64**

Recipe tested with Ubuntu 24.04 LTS, MATLAB 2024b
- dependencies:\
openblas 0.3.26 \
mumps    5.6.2 \
metis    5.1.0 \
scotch   7.0.3 \
mumps    5.6.2

- package requirements:
bison flex zlib-dev gfortran

- compiler requirements: \
gcc-9 gfortran-9 \
https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/system-requirements-release-2021b-supported-compilers.pdf \
may use update-alternatives under Linux: \
https://linuxconfig.org/how-to-switch-between-multiple-gcc-and-g-compiler-versions-on-ubuntu-20-04-lts-focal-fossa \

- prepend before calling MATLAB: \
```export MATLABPATH="<Path-to-mex-files>":<sth-like-/usr/local/MATLAB>"```

Simply execute ```make_...``` in order (and hope that it will also work as simple ;-)

## Literature
P. R. Amestoy, I. S. Duff, J. Koster and J.-Y. L’Excellent,
A fully asynchronous multifrontal solver using distributed dynamic
scheduling, SIAM Journal on Matrix Analysis and Applications,
Vol 23, No 1, pp 15-41 (2001).

P. R. Amestoy, A. Buttari, J.-Y. L’Excellent and T. Mary,
Performance and scalability of the block low-rank multifrontal
factorization on multicore architectures,
ACM Transactions on Mathematical Software,
Vol 45, Issue 1, pp 2:1-2:26 (2019)

## License

Copyright (c) 2018, Jan Blechta (Original code) \
Copyright (c) 2026, Mathias Scheunert / TU Bergakademie Freiberg (Modifications and additions)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://gnu.org>.

---

Portions of this software are based on code originally created by Jan Blechta
under the BSD 2-Clause License:

Copyright (c) 2018, Jan Blechta. All rights reserved.
