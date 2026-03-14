#ifndef RESULTANT_NAIVE_H
#define RESULTANT_NAIVE_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

#include "../sylvester_matrix/sylvester_matrix.h"

int fmpz_poly_resultant_naive(
  fmpz_t resultant,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
);

#endif