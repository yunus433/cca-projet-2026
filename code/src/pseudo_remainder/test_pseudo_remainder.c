// gcc -O2 -Wall -Wextra -std=c11   pseudo_remainder.c test_pseudo_remainder.c -lflint -lgmp -lmpfr   -o test_pseudorem

#include <time.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

#include "pseudo_remainder.h"

static int test_fmpz_poly_pseudo_remainder(int degree, int count, int bits)
{
    flint_rand_t state;
    flint_randinit(state);
    flint_randseed(state, (ulong) time(NULL), (ulong) (time(NULL) ^ 0x9e3779b97f4a7c15ULL));

    fmpz_poly_t A, B;
    fmpz_poly_t R_my, R_ref;
    fmpz_poly_t Q_ref;                 // Q de la pseudo-division FLINT
    fmpz_poly_t lhs, rhs, tmp;         // pour l'identité (lc(B))^d * A = B*Q + R
    fmpz_t lcB, lcBpow;

    fmpz_poly_init(A);
    fmpz_poly_init(B);
    fmpz_poly_init(R_my);
    fmpz_poly_init(R_ref);
    fmpz_poly_init(Q_ref);
    fmpz_poly_init(lhs);
    fmpz_poly_init(rhs);
    fmpz_poly_init(tmp);

    fmpz_init(lcB);
    fmpz_init(lcBpow);

    for (int i = 0; i < count; i++)
    {
        flint_printf("Starting test %d...\n", i + 1);

        // A de degré <= degree (taille degree+1)
        fmpz_poly_randtest(A, state, degree + 1, bits);

        // B non nul
        do {
            fmpz_poly_randtest(B, state, degree + 1, bits);
        } while (fmpz_poly_is_zero(B));

        ulong d_my = 0;
        fmpz_poly_pseudo_remainder(R_my, &d_my, A, B);

        ulong d_ref = 0;
        fmpz_poly_pseudo_rem(R_ref, &d_ref, A, B);

        if (d_my != d_ref || !fmpz_poly_equal(R_my, R_ref))
        {
            flint_printf("Test %d FAILED (mismatch vs FLINT).\n", i + 1);
            flint_printf("A(x) = "); fmpz_poly_print_pretty(A, "x"); flint_printf("\n");
            flint_printf("B(x) = "); fmpz_poly_print_pretty(B, "x"); flint_printf("\n");
            flint_printf("d_my=%wu, d_ref=%wu\n", d_my, d_ref);
            flint_printf("R_my(x)  = "); fmpz_poly_print_pretty(R_my, "x"); flint_printf("\n");
            flint_printf("R_ref(x) = "); fmpz_poly_print_pretty(R_ref, "x"); flint_printf("\n");
            goto cleanup_fail;
        }

        // 2) Vérif identité de pseudo-division via FLINT (donne Q)
        //    fmpz_poly_pseudo_divrem(Q, R, &d, A, B)
        fmpz_poly_pseudo_divrem(Q_ref, R_ref, &d_ref, A, B);

        // lcBpow = (lc(B))^d_ref
        fmpz_set(lcB, fmpz_poly_lead(B));
        fmpz_pow_ui(lcBpow, lcB, d_ref);

        // lhs = lcBpow * A
        fmpz_poly_scalar_mul_fmpz(lhs, A, lcBpow);

        // rhs = B*Q_ref + R_ref
        fmpz_poly_mul(rhs, B, Q_ref);
        fmpz_poly_add(rhs, rhs, R_ref);

        if (!fmpz_poly_equal(lhs, rhs))
        {
            flint_printf("Test %d FAILED (pseudo-division identity).\n", i + 1);
            flint_printf("A(x) = "); fmpz_poly_print_pretty(A, "x"); flint_printf("\n");
            flint_printf("B(x) = "); fmpz_poly_print_pretty(B, "x"); flint_printf("\n");
            flint_printf("d=%wu\n", d_ref);
            flint_printf("Q_ref(x) = "); fmpz_poly_print_pretty(Q_ref, "x"); flint_printf("\n");
            flint_printf("R_ref(x) = "); fmpz_poly_print_pretty(R_ref, "x"); flint_printf("\n");
            flint_printf("(lcB)^d * A = "); fmpz_poly_print_pretty(lhs, "x"); flint_printf("\n");
            flint_printf("B*Q + R     = "); fmpz_poly_print_pretty(rhs, "x"); flint_printf("\n");
            goto cleanup_fail;
        }

        flint_printf("Test %d OK.\n", i + 1);
    }

    // OK
    fmpz_clear(lcB);
    fmpz_clear(lcBpow);
    fmpz_poly_clear(A);
    fmpz_poly_clear(B);
    fmpz_poly_clear(R_my);
    fmpz_poly_clear(R_ref);
    fmpz_poly_clear(Q_ref);
    fmpz_poly_clear(lhs);
    fmpz_poly_clear(rhs);
    fmpz_poly_clear(tmp);
    flint_randclear(state);
    return 0;

cleanup_fail:
    fmpz_clear(lcB);
    fmpz_clear(lcBpow);
    fmpz_poly_clear(A);
    fmpz_poly_clear(B);
    fmpz_poly_clear(R_my);
    fmpz_poly_clear(R_ref);
    fmpz_poly_clear(Q_ref);
    fmpz_poly_clear(lhs);
    fmpz_poly_clear(rhs);
    fmpz_poly_clear(tmp);
    flint_randclear(state);
    return 1;
}

int main(void)
{
    // ex: degree=12, 50 tests, coeff bits=10
    return test_fmpz_poly_pseudo_remainder(12, 500, 10);
}