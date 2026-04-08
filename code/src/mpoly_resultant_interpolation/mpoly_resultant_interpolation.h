#ifndef MPOLY_RESULTANT_INTERPOLATION
#define MPOLY_RESULTANT_INTERPOLATION

#include <stdio.h>
#include <time.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mpoly.h> 

int fmpz_mpoly_resultant_interpolation(
  fmpz_mpoly_t R,
  const fmpz_mpoly_t P,
  const fmpz_mpoly_t Q,
  slong var,
  const fmpz_mpoly_ctx_struct *ctx,
  const flint_bitcnt_t bits
);

#endif