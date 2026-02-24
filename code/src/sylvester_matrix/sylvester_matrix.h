#ifndef SYLVESTER_MATRIX_H
#define SYLVESTER_MATRIX_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

int fmpz_poly_sylvester_matrix(
  fmpz_mat_t S,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
);

#endif