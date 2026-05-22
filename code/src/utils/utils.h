#ifndef UTILS_H
#define UTILS_H

#include <string.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

void fmpz_poly_randtest_uniform(
  fmpz_poly_struct *f,
  flint_rand_struct *state,
  slong len,
  flint_bitcnt_t bits
);

int fmpz_poly_set_str_pretty(
  fmpz_poly_struct *poly,
  const char *str,
  const char *x
);

#endif