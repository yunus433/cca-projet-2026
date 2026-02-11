#ifndef EUCLIDEAN_H
#define EUCLIDEAN_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

int fmpz_poly_euclidean_division(
  fmpz_poly_t q,
  fmpz_poly_t r,
  const fmpz_poly_t a,
  const fmpz_poly_t b
);

#endif