Forked from: https://github.com/xmjiao/mumps4m-openmp and https://github.com/blechta/mumps-matlab-recipes \
Original Authors: Xiangmin Jiao, Jan Blechta \
Thanks to ChatGPT and Google-AI!

Builds 64-bit MUMPS and statically links dependent libraries

Recipe tested with Ubuntu 26.04 LTS, MATLAB 2026a
- dependencies:\
openblas 0.3.26 \
mumps    5.6.2 \
metis    5.1.0 \
scotch   7.0.3

- package requirements:
bison flex zlib-dev gfortran

- compiler requirements: \
If you get some error msg, please check: https://de.mathworks.com/support/requirements/supported-compilers-linux.html

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


