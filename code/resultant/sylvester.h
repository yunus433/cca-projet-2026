#ifndef SYLV_H
#define SYLV_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

int fmpz_poly_sylvester_matrix(
  fmpz_mat_t S, // Not a pointer, already accessed by referenced in Sylvester types. No dynamic memory allocation either, so no need to free later on
  const fmpz_poly_t P,
  const fmpz_poly_t Q
);

#endif