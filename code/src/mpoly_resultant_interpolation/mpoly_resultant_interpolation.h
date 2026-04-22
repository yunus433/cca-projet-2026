#ifndef MPOLY_RESULTANT_INTERPOLATION
#define MPOLY_RESULTANT_INTERPOLATION

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mpoly.h> 
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_mpoly.h> 
#include <flint/nmod_poly.h>
#include <flint/nmod_mpoly.h> 

typedef enum {
  BITSIZE_RANDOM,
  BITSIZE_RANDOM_POS,
  SMALL_ORDERED,
  SMALL_ORDERED_POS
} interpolation_mode_t;

int fmpz_mpoly_resultant_interpolation(
  fmpz_mpoly_t R,
  const fmpz_mpoly_t P,
  const fmpz_mpoly_t Q,
  slong var,
  const fmpz_mpoly_ctx_struct *ctx
);

int fmpz_mpoly_resultant_interpolation_mode(
  fmpz_mpoly_t R,
  const fmpz_mpoly_t P,
  const fmpz_mpoly_t Q,
  slong var,
  const fmpz_mpoly_ctx_struct *ctx,
  const flint_bitcnt_t bits,
  interpolation_mode_t mode 
);

int fmpq_mpoly_resultant_interpolation(
  fmpq_mpoly_t R,
  const fmpq_mpoly_t P,
  const fmpq_mpoly_t Q,
  slong var,
  const fmpq_mpoly_ctx_struct *ctx
);

int fmpq_mpoly_resultant_interpolation_mode(
  fmpq_mpoly_t R,
  const fmpq_mpoly_t P,
  const fmpq_mpoly_t Q,
  slong var,
  const fmpq_mpoly_ctx_struct *ctx,
  const flint_bitcnt_t bits,
  interpolation_mode_t mode 
);

int nmod_mpoly_resultant_interpolation(
  nmod_mpoly_t R,
  const nmod_mpoly_t P,
  const nmod_mpoly_t Q,
  slong var,
  const nmod_mpoly_ctx_struct *ctx
);

int nmod_mpoly_resultant_interpolation_mode(
  nmod_mpoly_t R,
  const nmod_mpoly_t P,
  const nmod_mpoly_t Q,
  slong var,
  const nmod_mpoly_ctx_struct *ctx,
  const flint_bitcnt_t bits,
  interpolation_mode_t mode 
);

#endif
