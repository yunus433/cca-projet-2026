# Fast Computation of Multivariate Subresultants
CCA Project #13 - Yunus GURLEK, Frédéric XIA

Supervisor: Weijia WANG

This repo is the research project and implementation of fast computation of resultants, subresultants and subresultant polynomials in multivariate polynomials.

The repo includes:

- The research paper explaining general findings and implementations, /report,
- A presentation of the research project, /presentation,
- Implementation and testing of Euclidean GCD for polynomials in Z, /code/src/euclid_gcd,
- Implementation and testing of Euclidean Division for polynomials in Q, /code/src/euclidean_division,
- Testing and performance analysis of FLINT interpolation in Z, /code/src/interpolation,
- Implementation, testing, and benchmarking of our evoluation & interpolation algorithm for Z, Q and Zp for resultants and subresultants, /code/src/mpoly_resultant_interpolation and /code/src/mpoly_subresultant_interpolation,
- Implementation and testing of univariate resultant with pseudo remainder algorithm in Z, /code/src/pseudo_remainder,
- Implementation and testing of univariate resultants with naive (determinant of Sylvester matrix) algorithm in Z, /code/src/resultant,
- Implementation and testing of univariate subresultants with naive, pseudo remainder and euclid (for Q) methods in Z, /code/src/subresultant,
- Implementation and testing of univariate subresultant polynomials with naive and pseudo remainder methods in Z, /code/src/subresultant_polynomial,
- Implementations and test of sylvester matrices used in naive methods in Z, /code/src/sylvester_matrix,
- Three utility functions and tests: Normalize polynomial in Z, random uniform polynomial generation in Z and fmpz_poly_set_str_pretty (missing in FLINT normally, used for testcases), /code/src/utils/*,

Implementation that are not mentioned here are either not complete or not very important for this project.
