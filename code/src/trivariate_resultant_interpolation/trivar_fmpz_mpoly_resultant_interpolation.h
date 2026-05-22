#ifndef TRIVAR_MPOLY_RESULTANT_INTERPOLATION
#define TRIVAR_MPOLY_RESULTANT_INTERPOLATION

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

int trivar_fmpz_mpoly_resultant_interpolation(
  fmpz_mpoly_t R,
  const fmpz_mpoly_t P,
  const fmpz_mpoly_t Q,
  slong var,
  const fmpz_mpoly_ctx_struct *ctx
);

#endif
