#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_mat.h>
#include <flint/nmod_poly.h>

int fmpq_poly_subres_euclid(fmpq_t *subres,const fmpq_poly_t P, const fmpq_poly_t Q){
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

    for(slong i=0;i<max_steps;i++){
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

    fmpq_t lc;
    fmpq_init(lc);

    while(!fmpq_poly_is_zero(b)){
        fmpq_poly_divrem(q,r1,a,b);
        slong deg_r1=fmpq_poly_degree(r1);
        fmpq_poly_get_coeff_fmpq(rho[rho_len],r1,deg_r1);
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
    
    //initialisation subres à 0
    for (slong k = 0; k <= d; k++) fmpq_zero(subres[k]);

    fmpq_pow_si(subres[n[1]], rho[1], n[0] - n[1]);  // sigma_{n1}
    fmpq_set(sigma_i, subres[n[1]]);

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
            fmpq_set(subres[n_ip1], sigma_ip1);
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

    for (slong i = 0; i < max_steps; i++) fmpq_clear(rho[i]);
    flint_free(rho);
    flint_free(n);

    fmpq_poly_clear(a);
    fmpq_poly_clear(b);
    return 0;
}

int main(void)
{
    fmpq_poly_t P, Q;
    fmpq_poly_init(P);
    fmpq_poly_init(Q);

    // P = x^3 + 2x^2 + 3x + 4
    fmpq_poly_set_coeff_si(P, 3, -3);
    fmpq_poly_set_coeff_si(P, 2, 8);
    fmpq_poly_set_coeff_si(P, 1, 2);
    fmpq_poly_set_coeff_si(P, 0, -5);
    fmpq_poly_set_coeff_si(P, 8, 1);
    fmpq_poly_set_coeff_si(P, 6, 1);
    fmpq_poly_set_coeff_si(P, 4, -3);



    // Q = 3x^2 + 2x + 1
    fmpq_poly_set_coeff_si(Q, 2, -4);
    fmpq_poly_set_coeff_si(Q, 1, -9);
    fmpq_poly_set_coeff_si(Q, 0, 21);
    fmpq_poly_set_coeff_si(Q, 4, 5);
    fmpq_poly_set_coeff_si(Q, 6, 3);


    slong d = fmpq_poly_degree(Q);   // degré de g
    fmpq_t *subres = flint_malloc(sizeof(fmpq_t) * (d + 1));
    for (slong k = 0; k <= d; k++) {
        fmpq_init(subres[k]);
    }

    int ret = fmpq_poly_subres_euclid(subres, P, Q);
    if (ret != 0) {
        printf("Erreur return -1 \n");
        return 1;
    }

    for (slong k = d; k >= 0; k--) {
        printf("sigma_%ld = ", k);
        fmpq_print(subres[k]);
        printf("\n");
    }

    for (slong k = 0; k <= d; k++) {
        fmpq_clear(subres[k]);
    }
    flint_free(subres);

    fmpq_poly_clear(P);
    fmpq_poly_clear(Q);

    flint_cleanup();
    return 0;
}