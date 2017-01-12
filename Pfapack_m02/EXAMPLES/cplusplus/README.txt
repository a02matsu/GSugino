--------
Overview
--------

The Fortran routines from PFAPACK can usually be called easily from
C or C++, too. The two examples here show how this can be achieved from
C++.

The C++ routines call the Fortran77 interface of PFAPACK. The header
file pfapack.h contains the necessary subroutine declarations.  When
calling Fortran from C, it is necessary to realize that Fortran passes
all of its parameters by reference, i.e. all arguments are actually
pointers to variables.

Interoperability between C and Fortran77 is not standardized, but
usually works (compiling the library and the program with compilers
from the same family, for example the GNU compiler family with g++ and
gfortran, often helps). Details depend on the specific architecture:
For example, usually the routine names in Fortran are appended by an
underscore (this is assumed in pfapack.h), but this can differ.  The
examples here can only serve as a guideline on how to call PFAPACK
from C.

-------------------------------
Overview of the example program
-------------------------------

example_dense.cc:
Compute the Pfaffian of a dense matrix

example_band.cc:
Compute the Pfaffian of a banded matrix
