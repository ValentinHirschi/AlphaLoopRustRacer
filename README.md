# AlphaLoopRustRacer
A repo for exploring alternative implementations of Local Unitarity integrands

# Install

Place this repository in the `alpha_loop` directory and generate the three test processes using the `<TEST_CASE>.aL` cards present in the repository.

If you want to be able to run arbitrary precision and quadruple precision, make sure to install Fortran `mpfr` with the tarball supplied as follows:

`tar -xzf mpfun20-mpfr-v24.tar.gz`
`cd mpfun20-mpfr-v24/fortran-var2`
`bash gnu-complib2.scr`
(if the last command complains, then edit the on-liner compilation in `gnu-complib2.scr` by adding your paths to `MPFR` and `GMP` if they are not standard)

# Run

Simply then run with:

`./rust_racer.py --n_inputs 2000 -a a_ddx_N3LO_SG_QG108`

If you want quadruple precisiona and arbitrary precision then the following four environment variables must be specified:

```
# e.g. : /Users/vjhirsch/HEP_programs/gmp-6.2.1_build
GMP_PATH
# e.g. : /Users/vjhirsch/HEP_programs/mpfr-4.1.0_build
MPFR_PATH
# e.g. : /Users/vjhirsch/HEP_programs/mpfr_fortran/mpfun20-mpfr-v24/fortran-var2
FORTRAN_MPFR_PATH
# e.g. : /Users/vjhirsch/MG5/3.0.2.py3/PLUGIN/tderivative_alphaloop/libraries/form/mppp-0.26
MPP_PATH
```

The current test cases available are:

`a_ddx_NLO_SG_QG0`
`a_ddx_NNLO_SG_QG3`
`a_ddx_N3LO_SG_QG108`
`a_ddx_N3LO_SG_QG58`
