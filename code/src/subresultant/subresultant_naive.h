#ifndef SUBRESULTANT_NAIVE_H
#define SUBRESULTANT_NAIVE_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

#include "../sylvester_matrix/sylvester_matrix.h"

int fmpz_poly_subresultant_naive(
  fmpz_t * subresultants,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
);

#endif