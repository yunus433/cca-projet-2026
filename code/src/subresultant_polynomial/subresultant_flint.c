#include "subresultant_flint.h"

void fmpz_poly_pseudo_rem_cohen_2(
    fmpz_poly_t R,
    const fmpz_poly_t A,
    const fmpz_poly_t B)
{
    if (fmpz_poly_is_zero(B))
    {
        flint_throw(FLINT_DIVZERO, "Division by zero\n");
    }

    // si deg(A) < deg(B)
    if (fmpz_poly_degree(A) < fmpz_poly_degree(B))
    {
        fmpz_poly_set(R, A);
        return;
    }

    fmpz_poly_set(R, A);

    slong degA = fmpz_poly_degree(R);
    slong degB = fmpz_poly_degree(B);

    fmpz_t lcB, lcR;
    fmpz_init(lcB);
    fmpz_init(lcR);

    fmpz_poly_get_coeff_fmpz(lcB, B, degB);

    slong e = degA - degB + 1;

    fmpz_poly_t tmp;
    fmpz_poly_init(tmp);

    while (degA >= degB && !fmpz_poly_is_zero(R))
    {
        // lc(R)
        fmpz_poly_get_coeff_fmpz(lcR, R, degA);

        // R = lc(B) * R
        fmpz_poly_scalar_mul_fmpz(R, R, lcB);

        // tmp = lc(R) * B
        fmpz_poly_scalar_mul_fmpz(tmp, B, lcR);

        // décalage (x^(degA-degB))
        fmpz_poly_shift_left(tmp, tmp, degA - degB);

        // R = R - tmp
        fmpz_poly_sub(R, R, tmp);
        // après sub
        fmpz_poly_set_coeff_si(R, degA, 0);

        // normalisation
        fmpz_poly_normalize(R);
        degA = fmpz_poly_degree(R);

        e--;
    }

    // multiplication finale par lc(B)^e
    fmpz_t pow;
    fmpz_init(pow);
    fmpz_pow_ui(pow, lcB, e);
    fmpz_poly_scalar_mul_fmpz(R, R, pow);
    fmpz_clear(pow);
    

    fmpz_poly_clear(tmp);
    fmpz_clear(lcB);
    fmpz_clear(lcR);
}

int fmpz_subresultant_flint(fmpz_poly_t *res,const fmpz_poly_t P,const fmpz_poly_t Q){
    fmpz_t a,b,g,h;
    fmpz_poly_t A,B,R;

    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(g);
    fmpz_init(h);
    slong count = 0;

    fmpz_poly_init(A);
    fmpz_poly_init(B);
    fmpz_poly_init(R);

    fmpz_poly_content(a,P);
    fmpz_poly_content(b,Q);

    if (fmpz_is_zero(a)||fmpz_is_zero(b)) {
        fprintf(stderr,"ERREUR 1\n");
        fflush(stderr);
        goto clear;
    }
    fmpz_poly_scalar_divexact_fmpz(A,P,a);
    fmpz_poly_scalar_divexact_fmpz(B,Q,b);
    if (fmpz_poly_degree(A) < fmpz_poly_degree(B)) {
        fmpz_poly_swap(A, B);
    }   

    fmpz_one(g);
    fmpz_one(h);

    fmpz_poly_set(res[count++], A);
    fmpz_poly_set(res[count++], B);

    while(!fmpz_poly_is_zero(B)){
        const slong delta = fmpz_poly_degree(A) - fmpz_poly_degree(B);

        fmpz_poly_pseudo_rem_cohen_2(R,A,B);
        fmpz_poly_normalize(R);
        if(fmpz_poly_is_zero(R)){
            break;
        }
        fmpz_poly_swap(A,B);
        fmpz_poly_set(B,R);
        
        if(delta==1){
            if (fmpz_is_zero(g)||fmpz_is_zero(h)) {
                printf("ERREUR 4\n");
                break;
            }
            fmpz_mul(b,g,h);
            if (fmpz_is_zero(b)) {
                fprintf(stderr,"ERREUR 2\n");
                fflush(stderr);
                break;
            }
            fmpz_poly_scalar_divexact_fmpz(B,B,b);
            fmpz_poly_normalize(B);
            fmpz_poly_get_coeff_fmpz(g,A,fmpz_poly_degree(A));
            fmpz_set(h,g);
        }
        else{
            if(delta<0){
                printf("JESUIS0\n");
            }
            fmpz_pow_ui(a,h,delta);
            if (fmpz_is_zero(a)||fmpz_is_zero(g)) {
                printf("ERREUR 4\n");
                break;
            }
            fmpz_mul(b,g,a);
            if (fmpz_is_zero(a)||fmpz_is_zero(b)) {
                fprintf(stderr,"ERREUR 3\n");
                fflush(stderr);
                break;
            }
            fmpz_poly_scalar_divexact_fmpz(B,B,b);
            fmpz_poly_normalize(B);
            
            fmpz_t lcA;
            fmpz_init(lcA);
            fmpz_poly_get_coeff_fmpz(lcA,A,fmpz_poly_degree(A));
            fmpz_pow_ui(b,lcA,delta);
            if (fmpz_is_zero(h)||fmpz_is_zero(b)) {
                printf("ERREUR 4\n");
                break;
            }
            fmpz_mul(g,h,b);
            if (fmpz_is_zero(a)||fmpz_is_zero(b)) {
                fprintf(stderr,"ERREUR 4\n");
                fflush(stderr);
                break;
            }
            fmpz_divexact(h,g,a);
            fmpz_poly_get_coeff_fmpz(g,A,fmpz_poly_degree(A));
            fmpz_clear(lcA);
        }
        fmpz_poly_set(res[count++], B);
    }
    clear:
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(g);
    fmpz_clear(h);

    fmpz_poly_clear(A);
    fmpz_poly_clear(B);
    fmpz_poly_clear(R);

    return 0;
}