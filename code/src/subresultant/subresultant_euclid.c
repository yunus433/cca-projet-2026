
#include "subresultant_euclid.h"

int fmpq_poly_subresultant_euclid(
    fmpq_t *subresultants,
    const fmpq_poly_t P,
    const fmpq_poly_t Q
) {
    if (fmpq_poly_is_zero(P) || fmpq_poly_is_zero(Q)){
        return -1;
    }
    fmpq_poly_t a, b;
    fmpq_poly_init(a);
    fmpq_poly_init(b);

    fmpq_t rho0;
    fmpq_init(rho0);
    slong dP=fmpq_poly_degree(P);
    fmpq_poly_get_coeff_fmpq(rho0,P,dP);
    fmpq_poly_make_monic(a, P);
    fmpq_t rho1;
    fmpq_init(rho1);
    slong d=fmpq_poly_degree(Q);
    fmpq_poly_get_coeff_fmpq(rho1,Q,d);
    fmpq_poly_make_monic(b, Q);

    if(fmpq_poly_degree(a)<fmpq_poly_degree(b)){
        fmpq_poly_swap(a,b); //a >= b
    }
    fmpq_t *rho;
    slong max_steps=d+1;

    rho=flint_malloc(sizeof(fmpq_t)*(max_steps+5));

    for(slong i=0;i<max_steps+5;i++){
        fmpq_init(rho[i]);
    }

    fmpq_set(rho[0],rho0); 
    fmpq_set(rho[1],rho1);
    slong rho_len=2;

    slong *n;
    n=flint_malloc(sizeof(slong)*(max_steps+5));
    n[0]=fmpq_poly_degree(a);
    n[1]=fmpq_poly_degree(b);
    slong n_len=2;

    fmpq_poly_t q,r1,r2; //r2 = monic version of r1
    fmpq_poly_init(q);
    fmpq_poly_init(r1);
    fmpq_poly_init(r2);


    while(!fmpq_poly_is_zero(b)){
        fmpq_poly_divrem(q,r1,a,b);
        slong deg_r1=fmpq_poly_degree(r1);
        fmpq_poly_get_coeff_fmpq(rho[rho_len],r1,deg_r1); // TODO degr1 peut etre -1 si r1 == 0
        rho_len++;
        fmpq_poly_make_monic(r2,r1);
        n[n_len]=fmpq_poly_degree(r2);
        n_len++;
        fmpq_poly_set(a,b);
        fmpq_poly_set(b,r2);
    }

    fmpq_t prod_rho, sigma_i, sigma_ip1, factor;
    fmpq_init(prod_rho);
    fmpq_init(sigma_i);
    fmpq_init(sigma_ip1);
    fmpq_init(factor);
    
    //initialisation subresultants à 0
    for (slong k = 0; k <= d; k++) fmpq_zero(subresultants[k]);

    fmpq_pow_si(subresultants[n[1]], rho[1], n[0] - n[1]);  // sigma_{n1}
    fmpq_set(sigma_i, subresultants[n[1]]);

    // prod_rho = rho0*rho1
    fmpq_mul(prod_rho, rho[0], rho[1]); //

    for (slong i = 1; i + 1 < n_len; i++) {
        slong n_im1 = n[i - 1];
        slong n_i   = n[i];
        slong n_ip1 = n[i + 1];

        // multiplier prod_rho par rho_{i+1}
        fmpq_mul(prod_rho, prod_rho, rho[i + 1]);

        slong exp_pow = n_i - n_ip1;  // >= 1 normalement

        // factor = prod_rho^(exp_pow)
        fmpq_pow_si(factor, prod_rho, exp_pow);

        // signe
        slong E = (n_im1 - n_i) * (n_i - n_ip1 + 1);
        if (E & 1) fmpq_neg(factor, factor);

        // sigma_{n_{i+1}} = factor * sigma_{n_i}
        fmpq_mul(sigma_ip1, factor, sigma_i);

        if (0 <= n_ip1 && n_ip1 <= d) {
            fmpq_set(subresultants[n_ip1], sigma_ip1);
        }

        fmpq_set(sigma_i, sigma_ip1);
    }

    fmpq_clear(prod_rho);
    fmpq_clear(sigma_i);
    fmpq_clear(sigma_ip1);
    fmpq_clear(factor);

    fmpq_poly_clear(q);
    fmpq_poly_clear(r1);
    fmpq_poly_clear(r2);

    for (slong i = 0; i < max_steps+5; i++) fmpq_clear(rho[i]);
    flint_free(rho);
    flint_free(n);

    fmpq_poly_clear(a);
    fmpq_poly_clear(b);
    fmpq_clear(rho0);
    fmpq_clear(rho1);
    return 0;
}
