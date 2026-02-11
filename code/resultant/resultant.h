#ifndef RESULTANT_H
#define RESULTANT_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

int fmpz_poly_resultant_naive(
  fmpz_t resultant,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
) ;
#endif