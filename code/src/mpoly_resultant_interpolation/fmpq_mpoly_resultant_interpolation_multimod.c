#include "mpoly_resultant_interpolation.h"

int
fmpq_mpoly_resultant_interpolation_multimod(
    fmpq_mpoly_t R,
    const fmpq_mpoly_t A,
    const fmpq_mpoly_t B,
    slong var,
    const fmpq_mpoly_ctx_struct *ctx)
{
    int success;
    fmpz_t degA, degB;
    fmpz_mpoly_ctx_struct * zctx;

    fmpz_init(degA);
    fmpz_init(degB);
    zctx = (fmpz_mpoly_ctx_struct *) ctx->zctx;

    success = (fmpz_mpoly_resultant_interpolation_multimod(
        R->zpoly, A->zpoly, B->zpoly, var, zctx) == 0);

    if (!success)
    {
        fmpq_zero(R->content);
        fmpz_clear(degA);
        fmpz_clear(degB);
        return -1;
    }

    fmpz_mpoly_degree_fmpz(degA, A->zpoly, var, zctx);
    fmpz_mpoly_degree_fmpz(degB, B->zpoly, var, zctx);

    if (fmpz_sgn(degA) >= 0 && fmpz_sgn(degB) >= 0)
    {
        fmpq_t t;
        fmpq_init(t);

        success = fmpq_pow_fmpz(t, A->content, degB)
               && fmpq_pow_fmpz(R->content, B->content, degA);

        if (success)
            fmpq_mul(R->content, R->content, t);
        else
            fmpq_zero(R->content);

        fmpq_clear(t);
    }
    else
    {
        fmpq_zero(R->content);
    }

    fmpq_mpoly_reduce(R, ctx);

    fmpz_clear(degA);
    fmpz_clear(degB);

    return success ? 0 : -1;
}

int
fmpq_mpoly_resultant_interpolation_multimod_mode(
    fmpq_mpoly_t R,
    const fmpq_mpoly_t A,
    const fmpq_mpoly_t B,
    slong var,
    const fmpq_mpoly_ctx_struct *ctx,
    const flint_bitcnt_t bits, /* This is only used for mode BITSIZE_RANDOM */
    interpolation_mode_t mode)
{
    int success;
    fmpz_t degA, degB;
    fmpz_mpoly_ctx_struct * zctx;

    fmpz_init(degA);
    fmpz_init(degB);
    zctx = (fmpz_mpoly_ctx_struct *) ctx->zctx;

    success = (fmpz_mpoly_resultant_interpolation_multimod_mode(
        R->zpoly, A->zpoly, B->zpoly, var, zctx, bits, mode) == 0);

    if (!success)
    {
        fmpq_zero(R->content);
        fmpz_clear(degA);
        fmpz_clear(degB);
        return -1;
    }

    fmpz_mpoly_degree_fmpz(degA, A->zpoly, var, zctx);
    fmpz_mpoly_degree_fmpz(degB, B->zpoly, var, zctx);

    if (fmpz_sgn(degA) >= 0 && fmpz_sgn(degB) >= 0)
    {
        fmpq_t t;
        fmpq_init(t);

        success = fmpq_pow_fmpz(t, A->content, degB)
               && fmpq_pow_fmpz(R->content, B->content, degA);

        if (success)
            fmpq_mul(R->content, R->content, t);
        else
            fmpq_zero(R->content);

        fmpq_clear(t);
    }
    else
    {
        fmpq_zero(R->content);
    }

    fmpq_mpoly_reduce(R, ctx);

    fmpz_clear(degA);
    fmpz_clear(degB);

    return success ? 0 : -1;
}
