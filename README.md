This is ZFP fork of https://github.com/LLNL/zfp with support of NEC SX-Aurora TSUBASA. Compression and decompression have been optimized for 1D,2D,3D and 4D arrays of type FP32 and FP64. Only fix bit rate mode is supported.

Compilation of the library requires [CMake](https://cmake.org/). To be able to compile tests, [CMake-VE](https://github.com/SX-Aurora/cmake/tree/nec-support/) is required.

    $ git clone https://github.com/LLNL/zfp.git
    $ export PATH=/opt/nec/ve/bin:/usr/local/ve/cmake-ve-3.20.2/bin:${PATH}
    $ cd zfp
    $ mkdir build
    $ cd build
    $ export CC=ncc
    $ export CXX=nc++

    # To compile with the tests
    $ cmake -DCMAKE_CXX_STANDARD=14 -DBUILD_EXAMPLES=ON -DZFP_WITH_VE=ON ../zfp

    # or to compile without the tests:
    $ cmake -DCMAKE_CXX_STANDARD=14 -DBUILD_EXAMPLES=ON -DBUILD_TESTING=OFF -DZFP_WITH_VE=ON ../zfp

    $ make -j


The library is compiled to run natively on NEC SX-Aurora TSUBASA, with support of following execution policies:
- zfp_exec_serial: serial execution, no vectorization.
- zfp_exec_omp: OpenMP multi-threaded execution, no vectorization.
- zfp_exec_ve: OpenMP multi-threaded execution with vectorization. If compiled with -DZFP_WITH_OPENMP=OFF, then run with vectorization only.



ZFP
===
[![Github Actions Build Status](https://github.com/LLNL/zfp/actions/workflows/main.yml/badge.svg?branch=develop)](https://github.com/LLNL/zfp/actions/workflows/main.yml)
[![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/qb3ld7j11segy52k/branch/develop?svg=true)](https://ci.appveyor.com/project/lindstro/zfp)
[![Documentation Status](https://readthedocs.org/projects/zfp/badge/?version=release0.5.5)](https://zfp.readthedocs.io/en/release0.5.5/?badge=release0.5.5)
[![Code Coverage](https://codecov.io/gh/LLNL/zfp/branch/develop/graph/badge.svg)](https://codecov.io/gh/LLNL/zfp)

zfp is a compressed format for representing multidimensional floating-point
and integer arrays.  zfp provides compressed-array classes that support high
throughput read and write random access to individual array elements.  zfp
also supports serial and parallel (OpenMP and CUDA) compression of whole
arrays, e.g., for applications that read and write large data sets to and
from disk.

zfp uses lossy but optionally error-bounded compression to achieve high
compression ratios.  Bit-for-bit lossless compression is also possible
through one of zfp's compression modes.  zfp works best for 2D, 3D, and 4D
arrays that exhibit spatial correlation, such as continuous fields from
physics simulations, natural images, regularly sampled terrain surfaces, etc.
zfp compression of 1D arrays is possible but generally discouraged.

zfp is freely available as open source and is distributed under a BSD license.
zfp is primarily written in C and C++ but also includes Python and Fortran
bindings.  zfp conforms to various language standards, including C89, C99,
C11, C++98, C++11, and C++14, and is supported on Linux, macOS, and Windows.


Quick Start
-----------

To download zfp, type:

    git clone https://github.com/LLNL/zfp.git

zfp may be built using either [CMake](https://cmake.org/) or
[GNU make](https://www.gnu.org/software/make/).  To use CMake, type:

    cd zfp
    mkdir build
    cd build
    cmake ..
    cmake --build . --config Release
    ctest

This builds the zfp library in the `build/lib` directory and the zfp
command-line executable in the `build/bin` directory.  It then runs
the regression tests.

zfp may also be built using GNU make:

    cd zfp
    make
    make test

Note: GNU builds are less flexible and do not support all available features,
e.g., CUDA support.

For further configuration and build instructions, please consult the
[documentation](https://zfp.readthedocs.io/en/latest/installation.html).


Documentation
-------------

Full HTML [documentation](http://zfp.readthedocs.io/) is available online.
A [PDF](http://readthedocs.org/projects/zfp/downloads/pdf/latest/) version
is also available.

Further information on the zfp software is included in these files:

- Change log: see [CHANGELOG.md](./CHANGELOG.md).
- Support and additional resources: see [SUPPORT.md](./SUPPORT.md).
- Code contributions: see [CONTRIBUTING.md](./CONTRIBUTING.md).


Authors
-------

zfp was originally developed by [Peter Lindstrom](https://people.llnl.gov/pl)
at [Lawrence Livermore National Laboratory](https://www.llnl.gov/).  Please
see the [Contributors Page](https://github.com/LLNL/zfp/graphs/contributors)
for a full list of contributors.

### Citing zfp

If you use zfp for scholarly research, please cite this paper:

* Peter Lindstrom.
  [Fixed-Rate Compressed Floating-Point Arrays](https://www.researchgate.net/publication/264417607_Fixed-Rate_Compressed_Floating-Point_Arrays).
  IEEE Transactions on Visualization and Computer Graphics, 20(12):2674-2683, December 2014.
  [doi:10.1109/TVCG.2014.2346458](http://doi.org/10.1109/TVCG.2014.2346458).

The algorithm implemented in the current version of zfp is described in the
[documentation](https://zfp.readthedocs.io/en/latest/algorithm.html) and in
the following paper:

* James Diffenderfer, Alyson Fox, Jeffrey Hittinger, Geoffrey Sanders, Peter Lindstrom.
  [Error Analysis of ZFP Compression for Floating-Point Data](https://www.researchgate.net/publication/324908266_Error_Analysis_of_ZFP_Compression_for_Floating-Point_Data).
  SIAM Journal on Scientific Computing, 41(3):A1867-A1898, June 2019.
  [doi:10.1137/18M1168832](http://doi.org/10.1137/18M1168832).


License
-------

zfp is distributed under the terms of the BSD 3-Clause license.  See
[LICENSE](./LICENSE) and [NOTICE](./NOTICE) for details.

SPDX-License-Identifier: BSD-3-Clause

LLNL-CODE-663824
