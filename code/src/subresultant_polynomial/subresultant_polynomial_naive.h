#ifndef SUBRESULTANT_POLYNOMIAL_NAIVE
#define SUBRESULTANT_POLYNOMIAL_NAIVE

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_mat.h>

#include "../sylvester_matrix/polynomial_sylvester_matrix.c"

int fmpz_subresultant_polynomials_naive(
  fmpz_poly_t * subresultant_polynomials,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
);

#endif