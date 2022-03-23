# rminimax

Source code for designing flexible polynomial and rational approximations to mathematical functions.

## Requirements

This tool requires the following tools and/or libraries be installed and discoverable on your system:

* [CMake](https://cmake.org) - open-source, cross-platform family of tools designed to build, test and package software (v3.1 or later)
* [GNUplot](http://www.gnuplot.info/) - portable command-line driven graphing utility
* [Eigen](https://eigen.tuxfamily.org/) - C++ template library for linear algebra operations
* [GMP](https://gmplib.org/) - C library for arbitrary precision arithmetic on signed integers, rational numbers, and floating-point values
* [MPFR](https://www.mpfr.org/) - C library for multiple precision floating-point computations with correct rounding
* [MPREAL](https://github.com/advanpix/mpreal) - thin C++ wrapper for MPFR
* [fplll](https://github.com/fplll/fplll) - C++ library implementing several lattice algorithms
* [QSopt Exact](https://github.com/jonls/qsopt-ex) - an exact linear programming solver, written in C
* [flint](https://flintlib.org/) - C library for doing number theory
* [Arb](https://fredrikj.net/arb/) - C library for arbitrary-precision floating-point ball arithmetic

## Compilation
This is a CMake-based project, hence building the tool can be as simple as 
```
mkdir build
cd build
cmake ..
make
```

The compilation process generates the `rminimax` library that can be used to call various functions 
for generating flexible polynomial/rational approximations. It also generates the `ratapprox` tool, which uses 
this library and lets the user specify what kind of polynomial/rational approximations to generate from 
a convenient command line setting.

## Using the tool
### Input
The `ratapprox` tool is command line-based. The full list of options and parameters available can be 
recovered via help:
```
./ratapprox --help
```
which should return something similar to:
```
usage: ratapprox [OPTIONS]

--------------------------------------------------------------------------------
General options:
Option:                       Meaning:
--help                        Prints this help.
--log                         Prints log information during execution.
--function=[string]           The function to approximate.
--weight=[string]             Weight function. By default this is
                              the reciprocal of the function.
--num=[even|odd]              Specifies if the numerator should only
                              contain even/odd powered monomials.
--num=[string,string,...]     Specifies if the numerator should use a
                              custom basis where the functions are
                              given as parsable strings.
--den=[even|odd]              Specifies if the denominator should
                              contain only even/odd monomials.
--den=[string,string,...]     Specifies if the denominator should use
                              a custom basis where the functions are
                              given as parsable strings.
--denF=[int|string,...]       Specifies the list of floating point
                              formats for the den. coefficients.
                              They can either be numeric values or
                              strings, following Sollya notation:
                                  * HP - halfprecision
                                  * SG - singleprecision
                                  * D  - doubleprecision
                                  * DE - double extended
                                  * DD - double double
                                  * TD - triple double
                                  * QD - quad precision
                              In case it contains less elements than
                              the number of den. basis functions,
                              the last format is used for all
                              remaining coefficients. If it contains
                              more elements, then the extra formats
                              are discarded. By default, elements
                              will be optimized as double (D) values.
--numF=[int|string,...]       Similar to denF, only applied to the
                              numerator.
--type=[int,int]              Type of the rational approximation.
                              This option needs to be specified if
                              the monomial basis is used for the
                              numerator and denominator, or even/odd
                              alternatives.
--factorize                   Factorizes both the expressions in the
                              numerator and denominator into irreducible
                              degree one and degree two factors.
                              WARNING: results will be irrelevant if the
                              numerator and denominator are not polynomials
                               with a full set of basis functions.
--factorF=[string]            Specifies the floating point format for
                              the factorization coefficients. Can take
                              the same values as in numF and denF.
--dom=[double,double]         Approximation interval.
--output=[string]             Path to the output file.
                              By default string=./coeffs.sollya
```

### Output

The tool generates a Sollya file containing two lists with the coefficients of the final 
approximation that is generated. The first list contains the coefficients of the numerator, 
whereas the second includes the coefficients of the denominator.