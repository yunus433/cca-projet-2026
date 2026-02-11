#ifndef SUBRES_H
#define SUB_RES_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

int fmpz_poly_subresultants_naive(
  fmpz_t * subresultants,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
);

#endif