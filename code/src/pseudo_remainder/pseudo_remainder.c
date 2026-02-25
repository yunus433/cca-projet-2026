#include "pseudo_remainder.h"

static int poly_is_monomial(const fmpz_poly_t P)
{
    slong len = fmpz_poly_length(P);
    int nz = 0;

    for (slong i = 0; i < len; i++)
    {
        const fmpz *ci = fmpz_poly_get_coeff_ptr(P, i);
        if (ci && !fmpz_is_zero(ci))
            if (++nz > 1) return 0;
    }
    return nz == 1;
}

void fmpz_poly_pseudo_remainder(
    fmpz_poly_t R,
    ulong *d,
    const fmpz_poly_t A,
    const fmpz_poly_t B
) { //R_{k+1}=lcB*R_k - lcR*x^shift*B
    fmpz_poly_t Q, term, tmp, Bshift;
    fmpz_t lcB, lcR,lcBabs;

    fmpz_poly_init(Q);
    fmpz_poly_init(term);
    fmpz_poly_init(tmp);
    fmpz_poly_init(Bshift);
    fmpz_init(lcB);
    fmpz_init(lcR);
    fmpz_init(lcBabs);
    fmpz_poly_set(R, A);
    fmpz_poly_zero(Q);

    fmpz_set(lcB, fmpz_poly_lead(B));
    fmpz_abs(lcBabs,lcB);
    *d = 0;


    /* Cas trivial : si A = 0 => R = 0, d = 0 */
    if (fmpz_poly_is_zero(A))
    {
        fmpz_poly_zero(R);
        *d = 0;
        goto cleanup;
    }

    /* Cas trivial : si deg(A) < deg(B) => pas de pseudo-division */
    if (fmpz_poly_degree(A) < fmpz_poly_degree(B))
    {
        fmpz_poly_set(R, A);
        *d = 0;
        goto cleanup;
    }
        /* Cas monôme/constant : branche spéciale */
    if (poly_is_monomial(B))
    {
        slong k = fmpz_poly_degree(B);   // B = lcB * x^k

        /* si lcB = ±1 => d=0 (convention FLINT observée) */
        if (fmpz_is_one(lcB) || fmpz_cmp_si(lcB, -1) == 0)
        {
            *d = 0;
            /* R = A mod x^k (reste de la division par x^k) */
            fmpz_poly_set(R, A);
            fmpz_poly_truncate(R, k);
            goto cleanup;
        }

        /* sinon : d=1 et R = lcB*(A mod x^k) */
        *d = 1;
        fmpz_poly_set(R, A);
        fmpz_poly_truncate(R, k);
        fmpz_poly_scalar_mul_fmpz(R, R, lcB);
        goto cleanup;
    }
    const slong degB = fmpz_poly_degree(B);
    while (!fmpz_poly_is_zero(R) && fmpz_poly_degree(R) >= degB)
    {
        const slong degR = fmpz_poly_degree(R);
        const slong s = degR - degB;

        fmpz_set(lcR, fmpz_poly_lead(R));

        /* term = c*x^s */
        fmpz_poly_zero(term);
        fmpz_poly_set_coeff_fmpz(term, s, lcR);

        /* Q = lcB*Q + term */
        fmpz_poly_scalar_mul_fmpz(Q, Q, lcB);
        fmpz_poly_add(Q, Q, term);

        /* R = lcB*R - term*B */
        fmpz_poly_scalar_mul_fmpz(R, R, lcB);

        /* tmp = term * B = (c*x^s)*B */
        fmpz_poly_shift_left(Bshift, B, s);     
        fmpz_poly_scalar_mul_fmpz(tmp, Bshift, lcR); /* tmp = c*(B*x^s) */

        fmpz_poly_sub(R, R, tmp);

        (*d)++;
    }
    /* 1) Convention pour lcB = ±1 */
    if (fmpz_is_one(lcB) || fmpz_cmp_si(lcB, -1) == 0)
    {
        if (fmpz_cmp_si(lcB, -1) == 0 && ((*d) & 1UL))
            fmpz_poly_neg(R, R);
        *d = 0;
    }
    else if (!fmpz_poly_is_zero(R))
    {
        /* 2) Normalisation : si content(R) divisible par |lcB|, on divise par lcB et on décrémente d */
        fmpz_t cont;
        fmpz_init(cont);

        while (*d > 0)
        {
            fmpz_poly_content(cont, R);

            if (!fmpz_divisible(cont, lcBabs))
                break;

            fmpz_poly_scalar_divexact_fmpz(R, R, lcB);  /* division par lcB SIGNÉ */
            (*d)--;
        }

        fmpz_clear(cont);
    }
    cleanup:
    /* NB : si R==0, on ne touche pas à d ici */  
    fmpz_clear(lcB);
    fmpz_clear(lcR);

    fmpz_poly_clear(Q);
    fmpz_poly_clear(term);
    fmpz_poly_clear(tmp);
    fmpz_poly_clear(Bshift);
    fmpz_clear(lcBabs);
}

/*int main(void)
{
    fmpz_poly_t f, g, R, pseudorem;
    ulong d, a;

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_init(R);
    fmpz_poly_init(pseudorem);

    fmpz_poly_set_coeff_si(g, 0, -5);
    fmpz_poly_set_coeff_si(g, 1, 2);
    fmpz_poly_set_coeff_si(g, 2, 8);
    fmpz_poly_set_coeff_si(g, 3, -3);
    fmpz_poly_set_coeff_si(g, 4, -3);
    fmpz_poly_set_coeff_si(g, 6, 1);
    fmpz_poly_set_coeff_si(g, 8, 1);

    fmpz_poly_set_coeff_si(f, 0, 21);
    fmpz_poly_set_coeff_si(f, 1, -9);
    fmpz_poly_set_coeff_si(f, 2, -4);
    fmpz_poly_set_coeff_si(f, 4, 5);
    fmpz_poly_set_coeff_si(f, 6, 3);

    fmpz_poly_pseudo_rem_naive(R, &d, f, g);

    flint_printf("A(x) = "); fmpz_poly_print_pretty(f, "x"); flint_printf("\n");
    flint_printf("B(x) = "); fmpz_poly_print_pretty(g, "x"); flint_printf("\n");
    flint_printf("Pseudo-remainder naive R(x) = ");
    fmpz_poly_print_pretty(R, "x");
    flint_printf("\n");

    flint_printf("Vrai pseudo remainder (FLINT) :\n");
    fmpz_poly_pseudo_rem(pseudorem, &a, f, g);
    fmpz_poly_print_pretty(pseudorem, "x");
    flint_printf("\n");

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(R);
    fmpz_poly_clear(pseudorem);

    flint_cleanup();
    return 0;
}
*/