#ifndef MPOLY_SUBRESULTANT_INTERPOLATION
#define MPOLY_SUBRESULTANT_INTERPOLATION

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

#include "../subresultant/subresultant_pseudo_remainder.h"

int fmpz_mpoly_subresultant_interpolation(
  fmpz_mpoly_t * subresultant_polynomials,
  const fmpz_mpoly_t P,
  const fmpz_mpoly_t Q,
  slong var,
  const fmpz_mpoly_ctx_struct *ctx
);

#endif