#include "mpoly_resultant_interpolation.h"

int fmpq_mpoly_resultant_interpolation(
  fmpq_mpoly_t R,
  const fmpq_mpoly_t P,
  const fmpq_mpoly_t Q,
  slong var,
  const fmpq_mpoly_ctx_struct *ctx
) {
  int success;
  fmpz_t degP, degQ;
  fmpz_mpoly_ctx_struct *zctx;

  fmpz_init(degP);
  fmpz_init(degQ);
  zctx = (fmpz_mpoly_ctx_struct *) ctx->zctx;

  success = (fmpz_mpoly_resultant_interpolation(
    R->zpoly, P->zpoly, Q->zpoly, var, zctx) == 0);

  if (!success) {
    fmpq_zero(R->content);
    fmpz_clear(degP);
    fmpz_clear(degQ);
    return -1;
  }

  fmpz_mpoly_degree_fmpz(degP, P->zpoly, var, zctx);
  fmpz_mpoly_degree_fmpz(degQ, Q->zpoly, var, zctx);

  if (fmpz_sgn(degP) >= 0 && fmpz_sgn(degQ) >= 0) {
    fmpq_t t;
    fmpq_init(t);

    success = fmpq_pow_fmpz(t, P->content, degQ)
           && fmpq_pow_fmpz(R->content, Q->content, degP);

    if (success)
      fmpq_mul(R->content, R->content, t);
    else
      fmpq_zero(R->content);

    fmpq_clear(t);
  } else {
    fmpq_zero(R->content);
  }

  fmpq_mpoly_reduce(R, ctx);

  fmpz_clear(degP);
  fmpz_clear(degQ);

  return success ? 0 : -1;
}

int fmpq_mpoly_resultant_interpolation_mode(
  fmpq_mpoly_t R,
  const fmpq_mpoly_t P,
  const fmpq_mpoly_t Q,
  slong var,
  const fmpq_mpoly_ctx_struct *ctx,
  const flint_bitcnt_t bits, // This is only used for mode BITSIZE_RANDOM
  interpolation_mode_t mode
) {
  int success;
  fmpz_t degP, degQ;
  fmpz_mpoly_ctx_struct *zctx;

  fmpz_init(degP);
  fmpz_init(degQ);
  zctx = (fmpz_mpoly_ctx_struct *) ctx->zctx;

  success = (fmpz_mpoly_resultant_interpolation_mode(
    R->zpoly, P->zpoly, Q->zpoly, var, zctx, bits, mode) == 0);

  if (!success) {
    fmpq_zero(R->content);
    fmpz_clear(degP);
    fmpz_clear(degQ);
    return -1;
  }

  fmpz_mpoly_degree_fmpz(degP, P->zpoly, var, zctx);
  fmpz_mpoly_degree_fmpz(degQ, Q->zpoly, var, zctx);

  if (fmpz_sgn(degP) >= 0 && fmpz_sgn(degQ) >= 0) {
    fmpq_t t;
    fmpq_init(t);

    success = fmpq_pow_fmpz(t, P->content, degQ)
           && fmpq_pow_fmpz(R->content, Q->content, degP);

    if (success)
      fmpq_mul(R->content, R->content, t);
    else
      fmpq_zero(R->content);

    fmpq_clear(t);
  } else {
    fmpq_zero(R->content);
  }

  fmpq_mpoly_reduce(R, ctx);

  fmpz_clear(degP);
  fmpz_clear(degQ);

  return success ? 0 : -1;
}
