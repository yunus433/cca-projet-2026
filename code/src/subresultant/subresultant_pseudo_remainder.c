#include "subresultant_pseudo_remainder.h"

int fmpz_poly_subresultant_pseudo_remainder(
    fmpz_t *subresultants,
    const fmpz_poly_t P,
    const fmpz_poly_t Q
) {
    fmpz_poly_t rim1, ri,rip1,pseudorem; /*im1=i-1, ip1=i+1 i=i*/
    fmpz_t gammai,gammaim1,beta,psi,psim1;
    fmpz_t tmp1,tmp2,tmp3;
    ulong di=0,dim1=0,d=0;
    slong i=1;

    fmpz_poly_init(rim1);
    fmpz_poly_init(ri);
    fmpz_poly_init(rip1);
    fmpz_poly_init(pseudorem);
    
    fmpz_init(gammai);
    fmpz_init(gammaim1);
    fmpz_init(beta);
    fmpz_init(psi);
    fmpz_init(psim1);

    fmpz_init(tmp1);
    fmpz_init(tmp2);
    fmpz_init(tmp3);

    fmpz_poly_set(rim1,P);
    fmpz_poly_set(ri,Q);
    if(fmpz_poly_degree(rim1)<fmpz_poly_degree(ri)){ //deg(P)>=deg(Q)
        fmpz_poly_swap(rim1,ri);
    }

    while(!fmpz_poly_is_zero(ri)){
        slong deg_im1=fmpz_poly_degree(rim1);
        slong deg_i=fmpz_poly_degree(ri);
        di=deg_im1-deg_i;
        fmpz_set(gammai,fmpz_poly_lead(ri));
        if(i==1){
            if(((di+1)%2)==0){
                fmpz_one(beta);
            }
            else{
                fmpz_set_si(beta,-1);
            }
            fmpz_set_si(psi,-1);
        }
        else{
            fmpz_neg(tmp1,gammaim1); 
            fmpz_pow_ui(tmp2,tmp1,dim1);
            if(dim1<=1){
                fmpz_one(tmp3);
            }
            else{
                fmpz_pow_ui(tmp3,psim1,dim1-1);
            }
            fmpz_divexact(psi,tmp2,tmp3);
            fmpz_pow_ui(tmp1,psi,di);
            fmpz_mul(tmp2,gammaim1,tmp1);
            fmpz_neg(beta,tmp2);
        }
        d=0;
        fmpz_poly_pseudo_rem(pseudorem,&d,rim1,ri);
        
        fmpz_poly_scalar_divexact_fmpz(rip1, pseudorem, beta);
        fmpz_poly_primitive_part(rip1, rip1);

        if(deg_im1-1>=0){
            fmpz_set(subresultants[deg_im1-1],gammai);
        }
        if(deg_i<deg_im1-1){
            ulong e=(ulong)(deg_im1-deg_i-1);
            fmpz_pow_ui(tmp1,gammai,e);  
            fmpz_set(subresultants[deg_i],tmp1);
        }

        fmpz_set(gammaim1,gammai);
        fmpz_set(psim1,psi);
        dim1=di;
        fmpz_poly_set(rim1,ri);
        fmpz_poly_set(ri,rip1);
        i++;
    }
    fmpz_poly_clear(rim1);
    fmpz_poly_clear(ri);
    fmpz_poly_clear(rip1);
    fmpz_poly_clear(pseudorem);

    fmpz_clear(gammai);
    fmpz_clear(gammaim1);
    fmpz_clear(beta);
    fmpz_clear(psi);
    fmpz_clear(psim1);

    fmpz_clear(tmp1);
    fmpz_clear(tmp2);
    fmpz_clear(tmp3);

    return 0;
}