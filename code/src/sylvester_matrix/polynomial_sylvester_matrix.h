#ifndef POLYNOMIAL_SYLVESTER_MATRIX
#define POLYNOMIAL_SYLVESTER_MATRIX

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_mat.h>

int fmpz_poly_polynomial_sylvester_matrix(
  fmpz_poly_mat_t S,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
);

#endif