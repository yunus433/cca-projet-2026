#include "subresultant_pseudo_remainder.h"

static void _clear_all(
    fmpz_poly_t rim1,
    fmpz_poly_t ri,
    fmpz_poly_t rip1,
    fmpz_poly_t prem,
    fmpz_t gammai,
    fmpz_t gammaim1,
    fmpz_t beta,
    fmpz_t psi,
    fmpz_t psim1,
    fmpz_t tmp1,
    fmpz_t tmp2,
    fmpz_t tmp3
) {
    fmpz_poly_clear(rim1);
    fmpz_poly_clear(ri);
    fmpz_poly_clear(rip1);
    fmpz_poly_clear(prem);

    fmpz_clear(gammai);
    fmpz_clear(gammaim1);
    fmpz_clear(beta);
    fmpz_clear(psi);
    fmpz_clear(psim1);

    fmpz_clear(tmp1);
    fmpz_clear(tmp2);
    fmpz_clear(tmp3);
}

int fmpz_poly_subresultant_pseudo_remainder(
    fmpz *subresultants,
    const fmpz_poly_t P,
    const fmpz_poly_t Q
) {
    fmpz_poly_t rim1, ri, rip1, prem;
    fmpz_t gammai, gammaim1, beta, psi, psim1;
    fmpz_t tmp1, tmp2, tmp3;
    slong max_deg, min_deg;
    slong deg_im1, deg_i, deg_ip1;
    ulong di = 0, dim1 = 0, d = 0;
    slong i = 1;

    if (fmpz_poly_is_zero(P) || fmpz_poly_is_zero(Q)) {
        return 1;
    }

    fmpz_poly_init(rim1);
    fmpz_poly_init(ri);
    fmpz_poly_init(rip1);
    fmpz_poly_init(prem);

    fmpz_init(gammai);
    fmpz_init(gammaim1);
    fmpz_init(beta);
    fmpz_init(psi);
    fmpz_init(psim1);

    fmpz_init(tmp1);
    fmpz_init(tmp2);
    fmpz_init(tmp3);

    fmpz_poly_set(rim1, P);
    fmpz_poly_set(ri, Q);

    if (fmpz_poly_degree(rim1) < fmpz_poly_degree(ri)) {
        fmpz_poly_swap(rim1, ri);
    }

    max_deg = fmpz_poly_degree(rim1);
    min_deg = fmpz_poly_degree(ri); // TODO: Update ça comme nécessaire!!!

    for (slong k = 0; k <= max_deg; k++) {
    fmpz_zero(subresultants + k);
    }

    {
        slong p = fmpz_poly_degree(rim1);
        slong q = fmpz_poly_degree(ri);

        if (q >= 0 && q <= min_deg) {
            ulong gap = (ulong)(p - q);
            fmpz_pow_ui(tmp1, fmpz_poly_lead(ri), gap);
            fmpz_set(subresultants + q, tmp1);
        }
    }
    

    while (!fmpz_poly_is_zero(ri)) {
        deg_im1 = fmpz_poly_degree(rim1);
        deg_i   = fmpz_poly_degree(ri);

        if (deg_im1 < deg_i) {
            _clear_all(rim1, ri, rip1, prem,
                       gammai, gammaim1, beta, psi, psim1,
                       tmp1, tmp2, tmp3);
            return 3;
        }

        di = (ulong)(deg_im1 - deg_i);
        fmpz_set(gammai, fmpz_poly_lead(ri));
        
        if (i == 1) {
            if (((di + 1UL) & 1UL) == 0UL) {
                fmpz_one(beta);
            } else {
                fmpz_set_si(beta, -1);
            }
            fmpz_set_si(psi, -1);
        } else {
            fmpz_neg(tmp1, gammaim1);
            fmpz_pow_ui(tmp2, tmp1, dim1);

            if (dim1 <= 1UL) {
                fmpz_one(tmp3);
            } else {
                fmpz_pow_ui(tmp3, psim1, dim1 - 1UL);
            }

            if (fmpz_is_zero(tmp3)) {
                _clear_all(rim1, ri, rip1, prem,
                           gammai, gammaim1, beta, psi, psim1,
                           tmp1, tmp2, tmp3);
                return 4;
            }

            fmpz_divexact(psi, tmp2, tmp3);

            fmpz_pow_ui(tmp1, psi, di);
            fmpz_mul(tmp2, gammaim1, tmp1);
            fmpz_neg(beta, tmp2);
        }

        d = 0;
        fmpz_poly_pseudo_rem(prem, &d, rim1, ri);

        {
            ulong expected = (ulong)(deg_im1 - deg_i + 1);

            if (d > expected) {
                _clear_all(rim1, ri, rip1, prem,
                           gammai, gammaim1, beta, psi, psim1,
                           tmp1, tmp2, tmp3);
                return 5;
            }

            {
                ulong missing = expected - d;
                if (missing > 0 && !fmpz_poly_is_zero(prem)) {
                    fmpz_pow_ui(tmp1, gammai, missing);
                    fmpz_poly_scalar_mul_fmpz(prem, prem, tmp1);
                }
            }
        }

        if (fmpz_poly_is_zero(prem)) {
            fmpz_poly_zero(rip1);
            deg_ip1 = -1;
        } else {
            if (fmpz_is_zero(beta)) {
                _clear_all(rim1, ri, rip1, prem,
                           gammai, gammaim1, beta, psi, psim1,
                           tmp1, tmp2, tmp3);
                return 6;
            }

            fmpz_poly_scalar_divexact_fmpz(rip1, prem, beta);
            deg_ip1 = fmpz_poly_degree(rip1);

            if (deg_ip1 >= deg_i) {
                _clear_all(rim1, ri, rip1, prem,
                           gammai, gammaim1, beta, psi, psim1,
                           tmp1, tmp2, tmp3);
                return 7;
            }
        }

        if (deg_ip1 >= 0 && deg_ip1 <= min_deg) {
            // if (deg_i == deg_im1 - 1)
            //     fmpz_set(subresultants + (deg_i - 1), fmpz_poly_lead(rip1));
            // else {
                fmpz_pow_ui(subresultants + (deg_ip1), fmpz_poly_lead(rip1), deg_i - deg_ip1);
            // }
        }

        // for (slong j = deg_ip1 + 1; j < deg_i - 1; j++) {
        //     if (j >= 0 && j <= min_deg) {
        //         fmpz_zero(subresultants + j);
        //     }
        // }
        fmpz_set(gammaim1, gammai);
        fmpz_set(psim1, psi);
        dim1 = di;

        fmpz_poly_set(rim1, ri);
        fmpz_poly_set(ri, rip1);

        i++;
    }

    _clear_all(rim1, ri, rip1, prem,
               gammai, gammaim1, beta, psi, psim1,
               tmp1, tmp2, tmp3);

    return 0;
}