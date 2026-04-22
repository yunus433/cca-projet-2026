#include "mpoly_resultant_interpolation.h"

#include <flint/mpoly.h>
#include <flint/ulong_extras.h>

/*
   Set A <- B interpreted in symmetric range (-p/2, p/2].
   Assumes p fits in ulong.
*/
static void
_fmpz_mpoly_set_smod_from_nmod_mpoly(
    fmpz_mpoly_t A,
    const nmod_mpoly_t B,
    ulong p,
    const fmpz_mpoly_ctx_struct *ctx)
{
    slong i;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

    fmpz_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);

    for (i = 0; i < B->length; i++)
    {
        mpoly_monomial_set(A->exps + N*i, B->exps + N*i, N);
        fmpz_set_ui_smod(A->coeffs + i, B->coeffs[i], p);
    }

    _fmpz_mpoly_set_length(A, B->length, ctx);
}

/* Reduce integer sparse polynomial modulo p into nmod_mpoly. */
static void
_nmod_mpoly_set_fmpz_mpoly_modp(
    nmod_mpoly_t Ap,
    const fmpz_mpoly_t A,
    const nmod_mpoly_ctx_t pctx,
    const fmpz_mpoly_ctx_struct *zctx)
{
    slong i, k;
    slong N = mpoly_words_per_exp(A->bits, zctx->minfo);

    nmod_mpoly_fit_length_reset_bits(Ap, A->length, A->bits, pctx);

    k = 0;
    for (i = 0; i < A->length; i++)
    {
        ulong c = fmpz_get_nmod(A->coeffs + i, pctx->mod);
        if (c == 0)
            continue;

        mpoly_monomial_set(Ap->exps + N*k, A->exps + N*i, N);
        Ap->coeffs[k] = c;
        k++;
    }

    _nmod_mpoly_set_length(Ap, k, pctx);
}

static int
_nmod_mpoly_resultant_mod_prime_ctx(
    nmod_mpoly_t Rp,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    slong var,
    slong degA_expected,
    slong degB_expected,
    const nmod_mpoly_ctx_t pctx,
    const fmpz_mpoly_ctx_struct *zctx)
{
    int success;
    nmod_mpoly_t Ap, Bp;

    nmod_mpoly_init(Ap, pctx);
    nmod_mpoly_init(Bp, pctx);

    _nmod_mpoly_set_fmpz_mpoly_modp(Ap, A, pctx, zctx);
    _nmod_mpoly_set_fmpz_mpoly_modp(Bp, B, pctx, zctx);

    /* Prime precheck: preserve elimination-variable degrees modulo p. */
    if (nmod_mpoly_degree_si(Ap, var, pctx) != degA_expected ||
        nmod_mpoly_degree_si(Bp, var, pctx) != degB_expected)
    {
        success = 0;
        goto cleanup;
    }

    success = (nmod_mpoly_resultant_interpolation(Rp, Ap, Bp, var, pctx) == 0);

cleanup:
    nmod_mpoly_clear(Ap, pctx);
    nmod_mpoly_clear(Bp, pctx);

    return success;
}

static int
_nmod_mpoly_resultant_mod_prime_ctx_mode(
    nmod_mpoly_t Rp,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    slong var,
    slong degA_expected,
    slong degB_expected,
    const nmod_mpoly_ctx_t pctx,
    const fmpz_mpoly_ctx_struct *zctx,
    const flint_bitcnt_t bits,
    interpolation_mode_t mode)
{
    int success;
    nmod_mpoly_t Ap, Bp;

    nmod_mpoly_init(Ap, pctx);
    nmod_mpoly_init(Bp, pctx);

    _nmod_mpoly_set_fmpz_mpoly_modp(Ap, A, pctx, zctx);
    _nmod_mpoly_set_fmpz_mpoly_modp(Bp, B, pctx, zctx);

    /* Prime precheck: preserve elimination-variable degrees modulo p. */
    if (nmod_mpoly_degree_si(Ap, var, pctx) != degA_expected ||
        nmod_mpoly_degree_si(Bp, var, pctx) != degB_expected)
    {
        success = 0;
        goto cleanup;
    }

    success = (nmod_mpoly_resultant_interpolation_mode(Rp, Ap, Bp, var, pctx, bits, mode) == 0);

cleanup:
    nmod_mpoly_clear(Ap, pctx);
    nmod_mpoly_clear(Bp, pctx);

    return success;
}

/*
   H stores symmetric coefficients modulo M.
   Update to symmetric coefficients modulo M*p, matching Rp modulo p.
*/
static int
_fmpz_mpoly_crt_merge_with_nmod(
    fmpz_mpoly_t H,
    const fmpz_t M,
    const nmod_mpoly_t Rp,
    ulong p,
    const fmpz_mpoly_ctx_struct *ctx)
{
    int changed = 0;
    slong i, j, k;
    slong N;
    ulong * cmpmask;
    fmpz_mpoly_t T;
    fmpz_t crt, zero;

    FLINT_ASSERT(H->bits == Rp->bits);

    N = mpoly_words_per_exp(H->bits, ctx->minfo);
    cmpmask = FLINT_ARRAY_ALLOC(N, ulong);
    mpoly_get_cmpmask(cmpmask, N, H->bits, ctx->minfo);

    fmpz_mpoly_init(T, ctx);
    fmpz_mpoly_fit_length_reset_bits(T, H->length + Rp->length, H->bits, ctx);

    fmpz_init(crt);
    fmpz_init(zero);
    fmpz_zero(zero);

    i = 0;
    j = 0;
    k = 0;

    while (i < H->length || j < Rp->length)
    {
        int cmp;
        const ulong * exp;
        const fmpz * old_coeff;
        ulong rp_coeff;

        if (i >= H->length)
        {
            cmp = -1;
        }
        else if (j >= Rp->length)
        {
            cmp = 1;
        }
        else
        {
            cmp = mpoly_monomial_cmp(H->exps + N*i, Rp->exps + N*j, N, cmpmask);
        }

        if (cmp > 0)
        {
            exp = H->exps + N*i;
            old_coeff = H->coeffs + i;
            rp_coeff = 0;
            i++;
        }
        else if (cmp < 0)
        {
            exp = Rp->exps + N*j;
            old_coeff = NULL;
            rp_coeff = Rp->coeffs[j];
            j++;
        }
        else
        {
            exp = H->exps + N*i;
            old_coeff = H->coeffs + i;
            rp_coeff = Rp->coeffs[j];
            i++;
            j++;
        }

        fmpz_CRT_ui(crt, old_coeff ? old_coeff : zero, M, rp_coeff, p, 1);

        if (old_coeff)
            changed |= !fmpz_equal(crt, old_coeff);
        else
            changed |= !fmpz_is_zero(crt);

        if (!fmpz_is_zero(crt))
        {
            mpoly_monomial_set(T->exps + N*k, exp, N);
            fmpz_set(T->coeffs + k, crt);
            k++;
        }
    }

    _fmpz_mpoly_set_length(T, k, ctx);
    fmpz_mpoly_swap(H, T, ctx);

    fmpz_clear(crt);
    fmpz_clear(zero);
    fmpz_mpoly_clear(T, ctx);
    flint_free(cmpmask);

    return changed;
}

/* return: 1=match, 0=mismatch, -1=bad prime (skip) */
static int
_fmpz_mpoly_verify_candidate_mod(
    const fmpz_mpoly_t H,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    slong var,
    slong degA_expected,
    slong degB_expected,
    ulong p,
    const fmpz_mpoly_ctx_struct *zctx)
{
    int ok;
    int ret = -1;
    nmod_mpoly_ctx_t pctx;
    nmod_mpoly_t Rp, Hp;

    nmod_mpoly_ctx_init(pctx, zctx->minfo->nvars, zctx->minfo->ord, p);

    nmod_mpoly_init(Rp, pctx);
    nmod_mpoly_init(Hp, pctx);

    ok = _nmod_mpoly_resultant_mod_prime_ctx(
        Rp, A, B, var, degA_expected, degB_expected, pctx, zctx);
    if (!ok)
        goto cleanup;

    _nmod_mpoly_set_fmpz_mpoly_modp(Hp, H, pctx, zctx);

    if (Rp->bits > Hp->bits)
    {
        if (!nmod_mpoly_repack_bits_inplace(Hp, Rp->bits, pctx))
            goto cleanup;
    }
    else if (Rp->bits < Hp->bits)
    {
        if (!nmod_mpoly_repack_bits_inplace(Rp, Hp->bits, pctx))
            goto cleanup;
    }

    ret = nmod_mpoly_equal(Hp, Rp, pctx) ? 1 : 0;

cleanup:
    nmod_mpoly_clear(Rp, pctx);
    nmod_mpoly_clear(Hp, pctx);
    nmod_mpoly_ctx_clear(pctx);

    return ret;
}

int
fmpz_mpoly_resultant_interpolation_multimod(
    fmpz_mpoly_t R,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    slong var,
    const fmpz_mpoly_ctx_struct *ctx)
{
    int success = 0;
    slong attempts, good_primes, stable;
    slong degA_expected, degB_expected;
    ulong p;
    fmpz_t M;
    fmpz_mpoly_t Arep, Brep, H;
    flint_bitcnt_t bits;

    if (ctx->minfo->nvars != 2 || var < 0 || var >= 2)
        return -1;

    if (fmpz_mpoly_is_zero(A, ctx) || fmpz_mpoly_is_zero(B, ctx))
    {
        fmpz_mpoly_zero(R, ctx);
        return 0;
    }

    bits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_init(Arep, ctx);
    fmpz_mpoly_init(Brep, ctx);
    fmpz_mpoly_init(H, ctx);
    fmpz_init(M);

    if (!fmpz_mpoly_repack_bits(Arep, A, bits, ctx) ||
        !fmpz_mpoly_repack_bits(Brep, B, bits, ctx))
    {
        goto cleanup;
    }

    degA_expected = fmpz_mpoly_degree_si(Arep, var, ctx);
    degB_expected = fmpz_mpoly_degree_si(Brep, var, ctx);

    /* Start from primes just above 2^30 (about 31-bit primes). */
    p = UWORD(1) << 30;
    attempts = 0;
    good_primes = 0;
    stable = 0;

    while (attempts < 2500 && good_primes < 800)
    {
        int ok;
        nmod_mpoly_ctx_t pctx;
        nmod_mpoly_t Rp;

        p = n_nextprime(p, 1);
        attempts++;

        nmod_mpoly_ctx_init(pctx, ctx->minfo->nvars, ctx->minfo->ord, p);
        nmod_mpoly_init(Rp, pctx);

        ok = _nmod_mpoly_resultant_mod_prime_ctx(
            Rp, Arep, Brep, var, degA_expected, degB_expected, pctx, ctx);
        if (!ok)
        {
            nmod_mpoly_clear(Rp, pctx);
            nmod_mpoly_ctx_clear(pctx);
            continue;
        }

        if (good_primes == 0)
        {
            _fmpz_mpoly_set_smod_from_nmod_mpoly(H, Rp, p, ctx);
            fmpz_set_ui(M, p);
            good_primes = 1;
            stable = 0;
        }
        else
        {
            int changed;

            if (Rp->bits > H->bits)
            {
                if (!fmpz_mpoly_repack_bits_inplace(H, Rp->bits, ctx))
                {
                    nmod_mpoly_clear(Rp, pctx);
                    nmod_mpoly_ctx_clear(pctx);
                    goto cleanup;
                }
            }
            else if (Rp->bits < H->bits)
            {
                if (!nmod_mpoly_repack_bits_inplace(Rp, H->bits, pctx))
                {
                    nmod_mpoly_clear(Rp, pctx);
                    nmod_mpoly_ctx_clear(pctx);
                    goto cleanup;
                }
            }

            changed = _fmpz_mpoly_crt_merge_with_nmod(H, M, Rp, p, ctx);
            fmpz_mul_ui(M, M, p);
            good_primes++;
            stable = changed ? 0 : (stable + 1);
        }

        nmod_mpoly_clear(Rp, pctx);
        nmod_mpoly_ctx_clear(pctx);

        if (good_primes >= 3 && (stable >= 2 || (good_primes % 4) == 0))
        {
            slong checks = 0;
            slong tries = 0;
            ulong q = p;

            while (tries < 14 && checks < 2)
            {
                int v;
                q = n_nextprime(q, 1);
                tries++;

                v = _fmpz_mpoly_verify_candidate_mod(
                    H, Arep, Brep, var, degA_expected, degB_expected, q, ctx);
                if (v < 0)
                    continue;
                if (v == 0)
                {
                    checks = -1;
                    break;
                }
                checks++;
            }

            if (checks >= 2)
            {
                fmpz_mpoly_set(R, H, ctx);
                success = 1;
                break;
            }
        }
    }

cleanup:
    fmpz_mpoly_clear(Arep, ctx);
    fmpz_mpoly_clear(Brep, ctx);
    fmpz_mpoly_clear(H, ctx);
    fmpz_clear(M);

    return success ? 0 : -1;
}

int
fmpz_mpoly_resultant_interpolation_multimod_mode(
    fmpz_mpoly_t R,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    slong var,
    const fmpz_mpoly_ctx_struct *ctx,
    const flint_bitcnt_t bits_mode, /* This is only used for mode BITSIZE_RANDOM */
    interpolation_mode_t mode)
{
    int success = 0;
    slong attempts, good_primes, stable;
    slong degA_expected, degB_expected;
    ulong p;
    fmpz_t M;
    fmpz_mpoly_t Arep, Brep, H;
    flint_bitcnt_t bits;

    if (ctx->minfo->nvars != 2 || var < 0 || var >= 2)
        return -1;

    if (fmpz_mpoly_is_zero(A, ctx) || fmpz_mpoly_is_zero(B, ctx))
    {
        fmpz_mpoly_zero(R, ctx);
        return 0;
    }

    bits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_init(Arep, ctx);
    fmpz_mpoly_init(Brep, ctx);
    fmpz_mpoly_init(H, ctx);
    fmpz_init(M);

    if (!fmpz_mpoly_repack_bits(Arep, A, bits, ctx) ||
        !fmpz_mpoly_repack_bits(Brep, B, bits, ctx))
    {
        goto cleanup;
    }

    degA_expected = fmpz_mpoly_degree_si(Arep, var, ctx);
    degB_expected = fmpz_mpoly_degree_si(Brep, var, ctx);

    /* Start from primes just above 2^30 (about 31-bit primes). */
    p = UWORD(1) << 30;
    attempts = 0;
    good_primes = 0;
    stable = 0;

    while (attempts < 2500 && good_primes < 800)
    {
        int ok;
        nmod_mpoly_ctx_t pctx;
        nmod_mpoly_t Rp;

        p = n_nextprime(p, 1);
        attempts++;

        nmod_mpoly_ctx_init(pctx, ctx->minfo->nvars, ctx->minfo->ord, p);
        nmod_mpoly_init(Rp, pctx);

        ok = _nmod_mpoly_resultant_mod_prime_ctx_mode(
            Rp, Arep, Brep, var, degA_expected, degB_expected, pctx, ctx, bits_mode, mode);
        if (!ok)
        {
            nmod_mpoly_clear(Rp, pctx);
            nmod_mpoly_ctx_clear(pctx);
            continue;
        }

        if (good_primes == 0)
        {
            _fmpz_mpoly_set_smod_from_nmod_mpoly(H, Rp, p, ctx);
            fmpz_set_ui(M, p);
            good_primes = 1;
            stable = 0;
        }
        else
        {
            int changed;

            if (Rp->bits > H->bits)
            {
                if (!fmpz_mpoly_repack_bits_inplace(H, Rp->bits, ctx))
                {
                    nmod_mpoly_clear(Rp, pctx);
                    nmod_mpoly_ctx_clear(pctx);
                    goto cleanup;
                }
            }
            else if (Rp->bits < H->bits)
            {
                if (!nmod_mpoly_repack_bits_inplace(Rp, H->bits, pctx))
                {
                    nmod_mpoly_clear(Rp, pctx);
                    nmod_mpoly_ctx_clear(pctx);
                    goto cleanup;
                }
            }

            changed = _fmpz_mpoly_crt_merge_with_nmod(H, M, Rp, p, ctx);
            fmpz_mul_ui(M, M, p);
            good_primes++;
            stable = changed ? 0 : (stable + 1);
        }

        nmod_mpoly_clear(Rp, pctx);
        nmod_mpoly_ctx_clear(pctx);

        if (good_primes >= 3 && (stable >= 2 || (good_primes % 4) == 0))
        {
            slong checks = 0;
            slong tries = 0;
            ulong q = p;

            while (tries < 14 && checks < 2)
            {
                int v;
                q = n_nextprime(q, 1);
                tries++;

                v = _fmpz_mpoly_verify_candidate_mod(
                    H, Arep, Brep, var, degA_expected, degB_expected, q, ctx);
                if (v < 0)
                    continue;
                if (v == 0)
                {
                    checks = -1;
                    break;
                }
                checks++;
            }

            if (checks >= 2)
            {
                fmpz_mpoly_set(R, H, ctx);
                success = 1;
                break;
            }
        }
    }

cleanup:
    fmpz_mpoly_clear(Arep, ctx);
    fmpz_mpoly_clear(Brep, ctx);
    fmpz_mpoly_clear(H, ctx);
    fmpz_clear(M);

    return success ? 0 : -1;
}
