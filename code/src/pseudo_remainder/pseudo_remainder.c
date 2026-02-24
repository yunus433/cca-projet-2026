#include "pseudo_remainder.h"

void fmpz_poly_pseudo_remainder(
    fmpz_poly_t R,
    ulong *d,
    const fmpz_poly_t A,
    const fmpz_poly_t B
) { //R_{k+1}=lcB*R_k - lcR*x^shift*B
    fmpz_poly_t Q, term, tmp, Bshift;
    fmpz_t lcB, lcR;

    fmpz_poly_init(Q);
    fmpz_poly_init(term);
    fmpz_poly_init(tmp);
    fmpz_poly_init(Bshift);

    fmpz_init(lcB);
    fmpz_init(lcR);

    fmpz_poly_set(R, A);
    fmpz_poly_zero(Q);

    fmpz_set(lcB, fmpz_poly_lead(B));

    *d = 0;

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

    fmpz_clear(lcB);
    fmpz_clear(lcR);

    fmpz_poly_clear(Q);
    fmpz_poly_clear(term);
    fmpz_poly_clear(tmp);
    fmpz_poly_clear(Bshift);

    return 0;
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