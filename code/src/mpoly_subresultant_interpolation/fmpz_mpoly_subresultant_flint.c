#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_factor.h>
#include <flint/mpoly.h>
#include <flint/mpoly_univar.h>

/*
    Suite de sous-résultants bivariés.

    subresultants[i] :
        i-ème sous-résultant

    len :
        nombre de sous-résultants

    var :
        variable principale
*/

int mpoly_univar_subresultant_ducos(
    mpoly_univar_t * seq,
    slong * len,
    mpoly_univar_t B,
    mpoly_univar_t A,
    mpoly_void_ring_t R)
{
    *len = 0;
    slong i, j, k, aJ, ae;
    fmpz_t n, d, e, J, z1, alpha;
    int iexists, jexists, kexists;
    void * u, * v, * w, * s;
    mpoly_univar_t C, D, H, T;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(fmpz_cmp(B->exps + 0, A->exps + 0) >= 0);
    FLINT_ASSERT(fmpz_sgn(A->exps + 0) >= 0);

    if (fmpz_is_zero(A->exps + 0))
    {
        mpoly_univar_fit_length(seq[*len], 1, R);

        seq[*len]->length = 1;

        fmpz_zero(seq[*len]->exps + 0);

        POW_FMPZ(
            COEFF(seq[*len], 0),
            COEFF(A, 0),
            B->exps + 0);

        (*len)++;

        return 1;
    }

    fmpz_init(n);
    fmpz_init(d);
    fmpz_init(e);
    fmpz_init(J);
    fmpz_init(z1);
    fmpz_init(alpha);

    u = mpoly_void_ring_elem_init(R);
    v = mpoly_void_ring_elem_init(R);
    w = mpoly_void_ring_elem_init(R);
    s = mpoly_void_ring_elem_init(R);

    i = FLINT_MAX(B->length, A->length);
    mpoly_univar_init2(C, i + 1, R);
    mpoly_univar_init2(D, i + 1, R);
    mpoly_univar_init2(H, i + 1, R);
    mpoly_univar_init2(T, i + 1, R);


    fmpz_sub(z1, B->exps + 0, A->exps + 0);
    POW_FMPZ(s, A->coeffs + 0, z1);

    mpoly_univar_prem(B, A, D, R);

looper:

    if (B->length < 1)
        goto done;

    mpoly_univar_set(seq[*len], B, R);
    (*len)++;

    fmpz_set(d, A->exps + 0);
    fmpz_set(e, B->exps + 0);


    fmpz_sub(z1, d, e);
    if (fmpz_is_one(z1))
    {
        if (fmpz_is_zero(e))
            goto done;

        /* D = (B[e]*A - A[e]*B)/A[d] */
        /*           i        j       */
        i = 1;
        j = 1;
        if (A->length > 1 && fmpz_equal(A->exps + 1, e))
            i++;
        else
            j = B->length;
        D->length = 0;
        while (i < A->length || j < B->length)
        {
            mpoly_univar_fit_length(D, D->length + 1, R);

            if (i < A->length && j < B->length &&
                                          fmpz_equal(A->exps + i, B->exps + j))
            {
                MUL(u, COEFF(A, i), COEFF(B, 0));
                MUL(v, COEFF(A, 1), COEFF(B, j));
                SUB(w, u, v);
                DIVEXACT(COEFF(D, D->length), w, COEFF(A, 0));
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
                j++;
            }
            else if (i < A->length && (j >= B->length ||
                                       fmpz_cmp(A->exps + i, B->exps + j) > 0))
            {
                MUL(u, COEFF(A, i), COEFF(B, 0));
                DIVEXACT(COEFF(D, D->length), u, COEFF(A, 0));
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
            }
            else
            {
                FLINT_ASSERT((j < B->length && (i >= A->length ||
                                     fmpz_cmp(B->exps + j, A->exps + i) > 0)));

                MUL(v, COEFF(A, 1), COEFF(B, j));
                DIVEXACT(COEFF(D, D->length), v, COEFF(A, 0));
                NEG(COEFF(D, D->length), COEFF(D, D->length));
                fmpz_set(D->exps + D->length, B->exps + j);
                j++;
            }

            D->length += !IS_ZERO(COEFF(D, D->length));
        }

        /* A = (B[e]*(D - B*x) + B[e-1]*B)/s */
        /*            i    j            k    */
        i = 0;
        fmpz_sub_ui(z1, e, 1);
        if (B->length > 1 && fmpz_equal(B->exps + 1, z1))
        {
            j = 2;
            k = 1;
        }
        else
        {
            j = 1;
            k = B->length;
        }

        A->length = 0;
        while (i < D->length || j < B->length || k < B->length)
        {
            fmpz * exp;

            mpoly_univar_fit_length(A, A->length + 1, R);

            exp = A->exps + A->length;

            fmpz_zero(exp);

            if (i < D->length)
                _fmpz_max(exp, exp, D->exps + i);

            if (j < B->length)
            {
                fmpz_add_ui(z1, B->exps + j, 1);
                _fmpz_max(exp, exp, z1);
            }

            if (k < B->length)
                _fmpz_max(exp, exp, B->exps + k);

            iexists = (i < D->length) && fmpz_equal(exp, D->exps + i);
            jexists = (j < B->length) && fmpz_equal(exp, z1);
            kexists = (k < B->length) && fmpz_equal(exp, B->exps + k);

            FLINT_ASSERT(iexists || jexists || kexists);

            if (iexists)
            {
                if (jexists)
                {
                    SUB(w, COEFF(D, i), COEFF(B, j));
                    MUL(u, COEFF(B, 0), w);
                }
                else
                {
                    MUL(u, COEFF(B, 0), COEFF(D, i));
                }

                if (kexists)
                {
                    MUL(v, COEFF(B, 1), COEFF(B, k));
                    ADD(w, u, v);
                    DIVEXACT(COEFF(A, A->length), w, s);
                }
                else
                {
                    DIVEXACT(COEFF(A, A->length), u, s);
                }
            }
            else
            {
                if (kexists)
                {
                    MUL(u, COEFF(B, 1), COEFF(B, k));
                    if (jexists)
                    {
                        MUL(v, COEFF(B, 0), COEFF(B, j));
                        SUB(w, u, v);
                        DIVEXACT(COEFF(A, A->length), w, s);
                    }
                    else
                    {
                        DIVEXACT(COEFF(A, A->length), u, s);
                    }
                }
                else
                {
                    MUL(u, COEFF(B, 0), COEFF(B, j));
                    DIVEXACT(COEFF(A, A->length), u, s);
                    NEG(COEFF(A, A->length), COEFF(A, A->length));
                }
            }

            A->length += !IS_ZERO(COEFF(A, A->length));

            i += iexists;
            j += jexists;
            k += kexists;
        }

        mpoly_univar_swap(A, B);
        SET(s, COEFF(A, 0));
    }
    else
    {
        fmpz_sub(n, d, e);
        fmpz_sub_ui(n, n, 1);
        fmpz_one(alpha);
        while (fmpz_add(z1, alpha, alpha), fmpz_cmp(z1, n) <= 0)
            fmpz_set(alpha, z1);

        SET(u, COEFF(B, 0));
        fmpz_sub(n, n, alpha);
        while (fmpz_cmp_ui(alpha, 1) > 0)
        {
            fmpz_tdiv_q_2exp(alpha, alpha, 1);
            MUL(v, u, u);
            DIVEXACT(u, v, s);
            if (fmpz_cmp(n, alpha) >= 0)
            {
                MUL(v, u, COEFF(B, 0));
                DIVEXACT(u, v, s);
                fmpz_sub(n, n, alpha);
            }
        }

        mpoly_univar_fit_length(C, B->length, R);
        for (i = 0; i < B->length; i++)
        {
            MUL(v, u, COEFF(B, i));
            DIVEXACT(COEFF(C, i), v, s);
            fmpz_set(C->exps + i, B->exps + i);
        }
        C->length = B->length;


        if (fmpz_is_zero(e))
            goto done;

        /* H = C - C[e]*x^e */
        mpoly_univar_fit_length(H, C->length, R);
        for (i = 1; i < C->length; i++)
        {
            SET(COEFF(H, i - 1), COEFF(C, i));
            fmpz_set(H->exps + i - 1, C->exps + i);
        }
        H->length = C->length - 1;

        /* D = C[e]*A - A[e]*H  (truncated to powers of x < e) */
        i = 0;
        j = H->length;
        ae = A->length;
        while (i < A->length && fmpz_cmp(A->exps + i, e) >= 0)
        {
            if (fmpz_equal(A->exps + i, e))
            {
                j = 0;
                ae = i;
            }
            i++;
        }
        D->length = 0;
        while (i < A->length || j < H->length)
        {
            mpoly_univar_fit_length(D, D->length + 1, R);

            if (i < A->length && j < H->length &&
                                          fmpz_equal(A->exps + i, H->exps + j))
            {
                MUL(u, COEFF(A, i), COEFF(C, 0));
                MUL(v, COEFF(A, ae), COEFF(H, j));
                SUB(COEFF(D, D->length), u, v);
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
                j++;
            }
            else if (i < A->length && (j >= H->length ||
                                       fmpz_cmp(A->exps + i, H->exps + j) > 0))
            {
                MUL(COEFF(D, D->length), COEFF(A, i), COEFF(C, 0));
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
            }
            else
            {
                FLINT_ASSERT(j < H->length && (i >= A->length ||
                                      fmpz_cmp(H->exps + j, A->exps + i) > 0));

                MUL(COEFF(D, D->length), COEFF(A, ae), COEFF(H, j));
                NEG(COEFF(D, D->length), COEFF(D, D->length));
                fmpz_set(D->exps + D->length, H->exps + j);
                j++;
            }

            D->length += !IS_ZERO(COEFF(D, D->length));
        }

        for (fmpz_add_ui(J, e, 1); fmpz_cmp(J, d) < 0; fmpz_add_ui(J, J, 1))
        {
            if (H->length < 1)
                break;

            /* H = H*x - H[e-1]*B/B[e] */
            fmpz_sub_ui(z1, e, 1);
            if (fmpz_equal(H->exps + 0, z1))
            {
                i = 1;
                j = 1;
                T->length = 0;
                while (i < H->length || j < B->length)
                {
                    mpoly_univar_fit_length(T, T->length + 1, R);

                    if (i < H->length)
                        fmpz_add_ui(z1, H->exps + i, 1);

                    if (i < H->length && j < B->length &&
                                                   fmpz_equal(z1, B->exps + j))
                    {
                        MUL(u, COEFF(H, 0), COEFF(B, j));
                        DIVEXACT(v, u, COEFF(B, 0));
                        SUB(COEFF(T, T->length), COEFF(H, i), v);
                        fmpz_set(T->exps + T->length, B->exps + j);
                        i++;
                        j++;
                    }
                    else if (i < H->length && (j >= B->length ||
                                                fmpz_cmp(z1, B->exps + j) > 0))
                    {
                        SET(COEFF(T, T->length), COEFF(H, i));
                        fmpz_set(T->exps + T->length, z1);
                        i++;
                    }
                    else
                    {
                        FLINT_ASSERT(j < B->length && (i >= H->length ||
                                               fmpz_cmp(z1, B->exps + j) < 0));

                        MUL(u, COEFF(H, 0), COEFF(B, j));
                        DIVEXACT(COEFF(T, T->length), u, COEFF(B, 0));
                        NEG(COEFF(T, T->length), COEFF(T, T->length));
                        fmpz_set(T->exps + T->length, B->exps + j);
                        j++;
                    }

                    T->length += !IS_ZERO(COEFF(T, T->length));
                }

                mpoly_univar_swap(H, T);
            }
            else
            {
                FLINT_ASSERT(fmpz_cmp(H->exps + 0, z1) < 0);
                for (i = 0; i < H->length; i++)
                    fmpz_add_ui(H->exps + i, H->exps + i, 1);
            }

            /* find coefficient of x^J in A */
            aJ = 0;
            while (aJ < A->length && !fmpz_equal(A->exps + aJ, J))
                aJ++;
            if (aJ >= A->length)
                continue;

            /* D = D - A[J]*H */
            i = 0;
            j = 0;
            T->length = 0;
            while (i < D->length || j < H->length)
            {
                mpoly_univar_fit_length(T, T->length + 1, R);

                if (i < D->length && j < H->length &&
                                          fmpz_equal(D->exps + i, H->exps + j))
                {
                    MUL(u, COEFF(H, j), COEFF(A, aJ));
                    SUB(COEFF(T, T->length), COEFF(D, i), u);
                    fmpz_set(T->exps + T->length, D->exps + i);
                    i++;
                    j++;
                }
                else if (i < D->length && (j >= H->length ||
                                       fmpz_cmp(D->exps + i, H->exps + j) > 0))
                {
                    SET(COEFF(T, T->length), COEFF(D, i));
                    fmpz_set(T->exps + T->length, D->exps + i);
                    i++;
                }
                else
                {
                    FLINT_ASSERT(j < H->length && (i >= D->length ||
                                      fmpz_cmp(D->exps + i, H->exps + j) < 0));

                    MUL(COEFF(T, T->length), COEFF(H, j), COEFF(A, aJ));
                    NEG(COEFF(T, T->length), COEFF(T, T->length));
                    fmpz_set(T->exps + T->length, H->exps + j);
                    j++;
                }

                T->length += !IS_ZERO(COEFF(T, T->length));
            }
            mpoly_univar_swap(D, T);
        }

        /* B = (-1)^(d-e+1) * (B[e]*(D/A[d] - H*x) +  H[e-1]*B)/s */
        i = 0;
        fmpz_sub_ui(z1, e, 1);
        if (H->length > 0 && fmpz_equal(H->exps + 0, z1))
        {
            j = 1;
            k = 1;
        }
        else
        {
            j = 0;
            k = B->length;
        }
        T->length = 0;
        while (i < D->length || j < H->length || k < B->length)
        {
            fmpz * exp;

            mpoly_univar_fit_length(T, T->length + 1, R);

            exp = T->exps + T->length;
            fmpz_zero(exp);

            if (i < D->length)
                _fmpz_max(exp, exp, D->exps + i);

            if (j < H->length)
            {
                fmpz_add_ui(z1, H->exps + j, 1);
                _fmpz_max(exp, exp, z1);
            }

            if (k < B->length)
                _fmpz_max(exp, exp, B->exps + k);

            iexists = (i < D->length && fmpz_equal(exp, D->exps + i));
            jexists = (j < H->length && fmpz_equal(exp, z1));
            kexists = (k < B->length && fmpz_equal(exp, B->exps + k));

            FLINT_ASSERT(iexists || jexists || kexists);

            if (iexists)
            {
                if (jexists)
                {
                    DIVEXACT(u, COEFF(D, i), COEFF(A, 0));
                    SUB(w, u, COEFF(H, j));
                    MUL(u, COEFF(B, 0), w);
                }
                else
                {
                    DIVEXACT(u, COEFF(D, i), COEFF(A, 0));
                    MUL(u, COEFF(B, 0), u);
                }
                if (kexists)
                {
                    MUL(v, COEFF(H, 0), COEFF(B, k));
                    ADD(w, u, v);
                    DIVEXACT(COEFF(T, T->length), w, s);
                }
                else
                {
                    DIVEXACT(COEFF(T, T->length), u, s);
                }
            }
            else
            {
                if (kexists)
                {
                    MUL(u, COEFF(H, 0), COEFF(B, k));
                    if (jexists)
                    {
                        MUL(v, COEFF(B, 0), COEFF(H, j));
                        SUB(w, u, v);
                        DIVEXACT(COEFF(T, T->length), w, s);
                    }
                    else
                    {
                        DIVEXACT(COEFF(T, T->length), u, s);
                    }
                }
                else
                {
                    MUL(u, COEFF(B, 0), COEFF(H, j));
                    DIVEXACT(COEFF(T, T->length), u, s);
                    NEG(COEFF(T, T->length), COEFF(T, T->length));
                }
            }

            if (((fmpz_get_ui(d) - fmpz_get_ui(e)) & 1) == 0)
                NEG(COEFF(T, T->length), COEFF(T, T->length));

            T->length += !IS_ZERO(COEFF(T, T->length));

            i += iexists;
            j += jexists;
            k += kexists;
        }

        mpoly_univar_swap(B, T);
        mpoly_univar_swap(A, C);
        SET(s, COEFF(A, 0));
    }

    goto looper;

done:


    fmpz_clear(n);
    fmpz_clear(d);
    fmpz_clear(e);
    fmpz_clear(J);
    fmpz_clear(z1);
    fmpz_clear(alpha);
    mpoly_void_ring_elem_clear(u, R);
    mpoly_void_ring_elem_clear(v, R);
    mpoly_void_ring_elem_clear(w, R);
    mpoly_void_ring_elem_clear(s, R);
    mpoly_univar_clear(C, R);
    mpoly_univar_clear(D, R);
    mpoly_univar_clear(H, R);
    mpoly_univar_clear(T, R);
    return 1;
}


int
fmpz_mpoly_subresultant_sequence(
    fmpz_mpoly_t * subres,
    slong * len,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    slong var,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    /*
        Conversion FLINT :
            fmpz_mpoly -> fmpz_mpoly_univar
    */

    fmpz_mpoly_univar_t U0, U1;

    fmpz_mpoly_univar_init(U0, ctx);
    fmpz_mpoly_univar_init(U1, ctx);

    fmpz_mpoly_to_univar(U0, A, var, ctx);
    fmpz_mpoly_to_univar(U1, B, var, ctx);

    /*
        deg(U0) >= deg(U1)
    */

    if (U0->length != 0 &&
        U1->length != 0 &&
        fmpz_cmp(U0->exps + 0, U1->exps + 0) < 0)
    {
        fmpz_mpoly_univar_swap(U0, U1, ctx);
    }

    /*
        Conversion :
            fmpz_mpoly_univar
                    ->
            mpoly_univar
    */

    mpoly_void_ring_t R;

    mpoly_void_ring_init_fmpz_mpoly_ctx(R, ctx);

    mpoly_univar_t Ax, Bx;

    mpoly_univar_init(Ax, R);
    mpoly_univar_init(Bx, R);

    mpoly_univar_set_fmpz_mpoly_univar(
        Ax,
        R,
        U0,
        ctx);

    mpoly_univar_set_fmpz_mpoly_univar(
        Bx,
        R,
        U1,
        ctx);

    /*
        Tableau des sous-résultants
        en représentation interne
    */

    mpoly_univar_t seq[20];

    for (i = 0; i < 20; i++)
    {
        mpoly_univar_init(seq[i], R);
    }

    /*
        Calcul Ducos
    */

    *len = 0;

    mpoly_univar_subresultant_ducos(
        seq,
        len,
        Ax,
        Bx,
        R);

    /*
        Conversion inverse :
            mpoly_univar
                    ->
            fmpz_mpoly
    */

    for (i = 0; i < *len; i++)
    {
        fmpz_mpoly_univar_t T;

        fmpz_mpoly_univar_init(T, ctx);

        mpoly_univar_swap_fmpz_mpoly_univar(
            seq[i],
            R,
            T,
            ctx);

        _fmpz_mpoly_from_univar(
            subres[i],
            FLINT_BITS,
            T,
            var,
            ctx);

        fmpz_mpoly_univar_clear(T, ctx);
    }

    /*
        nettoyage
    */

    for (i = 0; i < 20; i++)
    {
        mpoly_univar_clear(seq[i], R);
    }

    mpoly_univar_clear(Ax, R);
    mpoly_univar_clear(Bx, R);

    fmpz_mpoly_univar_clear(U0, ctx);
    fmpz_mpoly_univar_clear(U1, ctx);
    mpoly_void_ring_clear(R);

    return 1;
}

int main(void)
{
    fmpz_mpoly_ctx_t ctx;

    /*
        2 variables : x,y
    */

    fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);

    fmpz_mpoly_t A, B;

    fmpz_mpoly_init(A, ctx);
    fmpz_mpoly_init(B, ctx);

    /*
        A(x,y)
    */

    fmpz_mpoly_set_str_pretty(
        A,
        "x^3 + 2*x^2*y + x*y + y^2",
        (const char *[]){"x","y"},
        ctx);

    /*
        B(x,y)
    */

    fmpz_mpoly_set_str_pretty(
        B,
        "x^2 + x*y + 1",
        (const char *[]){"x","y"},
        ctx);

    printf("A = ");
    fmpz_mpoly_print_pretty(
        A,
        (const char *[]){"x","y"},
        ctx);
    printf("\n");

    printf("B = ");
    fmpz_mpoly_print_pretty(
        B,
        (const char *[]){"x","y"},
        ctx);
    printf("\n\n");

    /*
        tableau sous-resultants
    */

    fmpz_mpoly_t S[20];

    slong i;
    slong len;

    for (i = 0; i < 20; i++)
        fmpz_mpoly_init(S[i], ctx);

    /*
        variable principale :
            var = 0 <=> x
    */

    fmpz_mpoly_subresultant_sequence(
        S,
        &len,
        A,
        B,
        0,
        ctx);

    printf("Suite de sous-resultants :\n\n");

    for (i = 0; i < len; i++)
    {
        printf("S[%ld] = ", i);

        fmpz_mpoly_print_pretty(
            S[i],
            (const char *[]){"x","y"},
            ctx);

        printf("\n");
    }

    /*
        nettoyage
    */

    for (i = 0; i < 20; i++)
        fmpz_mpoly_clear(S[i], ctx);

    fmpz_mpoly_clear(A, ctx);
    fmpz_mpoly_clear(B, ctx);

    fmpz_mpoly_ctx_clear(ctx);

    flint_cleanup();

    return 0;
}