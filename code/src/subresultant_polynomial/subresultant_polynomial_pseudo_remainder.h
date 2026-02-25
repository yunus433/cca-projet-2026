#ifndef SUBRESULTANT_POLYNOMIAL_PSEUDO_REMAINDER
#define SUBRESULTANT_POLYNOMIAL_PSEUDO_REMAINDER


#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_mat.h>

int fmpz_subresultant_polynomial_pseudo_remainder(
  fmpz_poly_t * subresultant_polynomials,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
);

#endif