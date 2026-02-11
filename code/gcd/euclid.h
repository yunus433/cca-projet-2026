#ifndef EUCLID_H
#define EUCLID_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

int fmpz_poly_gcd_euclid(
  fmpz_poly_t r,
  const fmpz_poly_t p,
  const fmpz_poly_t q
) ;
#endif