#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "euclidean-division/euclidean.h"
#include "gcd/euclid.h"
#include "pseudorem/pseudorem.h"
#include "resultant/subresultant2.h"
#include "resultant/subresultant.h"
#include "resultant/resultant.h"
#include "resultant/sylvester.h"

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_mat.h>
#include <flint/nmod_poly.h>

void generate_random_polys(fmpz_poly_t P,fmpz_poly_t Q,flint_rand_t state,slong degree,slong bits){
    fmpz_poly_randtest(P,state,degree+1,bits);
    do {
        fmpz_poly_randtest(Q,state,degree+1,bits);
    } while (fmpz_poly_is_zero(Q));
}

void generate_random_polysQ(fmpq_poly_t F,fmpq_poly_t G,flint_rand_t state,slong degree,slong bits){
    fmpq_poly_randtest(F,state,degree+1,bits);
    do {
        fmpq_poly_randtest(G,state,degree+1,bits);
    } while (fmpq_poly_is_zero(G));
}

int main(){
    int N=5000;
    fmpz_poly_t P,Q,q,r,R;
    fmpz_poly_init(P);
    fmpz_poly_init(Q);
    fmpz_poly_init(q);
    fmpz_poly_init(r);
    flint_rand_t state;
    flint_randinit(state);
    ulong s = (ulong) time(NULL);
    flint_randseed(state, s, s + 1UL);
    fmpz_poly_init(R);
    
    fmpq_poly_t F,G,q2,r2;
    fmpq_poly_init(F);
    fmpq_poly_init(G);
    fmpq_poly_init(q2);
    fmpq_poly_init(r2);

    slong degree=100;
    slong bits=50;
    clock_t start;
    clock_t end;
    clock_t add1=0;
    clock_t add2=0;
    fmpz_t resultant;
    fmpz_init(resultant);

    /*printf("TEST POUR LA FONCTION RESULTANTE\n");
    for(int i=0;i<N;i++){ //Environ 15 fois plus lent 
        generate_random_polys(P,Q,state,degree,bits);
        start=clock();
        fmpz_poly_resultant_naive(resultant,P,Q);
        end=clock();
        add1+=end-start;
        start=clock();
        fmpz_poly_resultant(resultant,P,Q);
        end=clock();
        add2+=end-start;
    }
    double time1=(double)add1/CLOCKS_PER_SEC;
    time1=time1/N;
    double time2=(double)add2/CLOCKS_PER_SEC;
    time2=time2/N;
    printf("Naive : %.6e sec\n", time1);
    printf("FLINT : %.6e sec\n", time2);
    printf("La fonction naive est %.2f fois plus lent que la fonction FLINT\n", time1 / time2);
    add1=0;
    add2=0;*/


    /*printf("TEST POUR LA FONCTION division euclidienne\n");
    for(int i=0;i<N;i++){ //ENVIRON 4,8 FOIS PLUS LENT
        generate_random_polysQ(F,G,state,degree,bits);
        start=clock();
        fmpq_poly_euclidean_division(q2,r2,F,G);
        end=clock();
        add1+=end-start;
        start=clock();
        fmpq_poly_divrem(q2,r2,F,G);
        end=clock();
        add2+=end-start;
    }
    double time1=(double)add1/CLOCKS_PER_SEC;
    time1=time1/N;
    double time2=(double)add2/CLOCKS_PER_SEC;
    time2=time2/N;
    printf("Naive : %.6e sec\n", time1);
    printf("FLINT : %.6e sec\n", time2);
    printf("La fonction naive est %.2f fois plus lent que la fonction FLINT\n", time1 / time2);
    add1=0;
    add2=0;*/
    
    printf("TEST POUR LA FONCTION PSEUDO REMAINDER\n"); //CHANGER ET FAIRE QUUNE ITERATION
    ulong d;
    for(int i=0;i<N;i++){ //sur les grandes instances on remarque que ca prend 3,75 fois plus de temps avec la fonction naïve
        generate_random_polys(P,Q,state,degree,bits);
        start=clock();
        fmpz_poly_pseudo_rem_naive(R,&d,P,Q);
        end=clock();
        add1+=end-start;
        start=clock();
        fmpz_poly_pseudo_rem(R,&d,P,Q);
        end=clock();
        add2+=end-start;
    }
    double time1=(double)add1/CLOCKS_PER_SEC;
    time1=time1/N;
    double time2=(double)add2/CLOCKS_PER_SEC;
    time2=time2/N;
    printf("Naive : %.6e sec\n", time1);
    printf("FLINT : %.6e sec\n", time2);
    printf("La fonction naive est %.2f fois plus lent que la fonction FLINT\n", time1 / time2);
    add1=0;
    add2=0;

    /*printf("TEST POUR LA FONCTION GCD\n"); //820 fois plus lent ???
    for(int i=0;i<N;i++){
        generate_random_polys(P,Q,state,degree,bits);
        start=clock();
        fmpz_poly_gcd_euclid(R,P,Q);
        end=clock();
        add1+=end-start;
        start=clock();
        fmpz_poly_gcd(R,P,Q);
        end=clock();
        add2+=end-start;
    }
    double time1=(double)add1/CLOCKS_PER_SEC;
    time1=time1/N;
    double time2=(double)add2/CLOCKS_PER_SEC;
    time2=time2/N;
    printf("Naive : %.6e sec\n", time1);
    printf("FLINT : %.6e sec\n", time2);
    printf("La fonction naive est %.2f fois plus lent que la fonction FLINT\n", time1 / time2);*/

    /*fmpz_poly_clear(P);
    fmpz_poly_clear(Q);
    fmpz_clear(resultant);
    flint_randclear(state);*/
    
    /*printf("TEST CORRECTNESS PSEUDO REM\n");

    generate_random_polys(P, Q, state, degree, bits);

    // éviter division par 0
    while (fmpz_poly_is_zero(Q))
        generate_random_polys(P, Q, state, degree, bits);

    fmpz_poly_t R1, R2;
    fmpz_poly_init(R1);
    fmpz_poly_init(R2);

    ulong d1, d2;

    fmpz_poly_pseudo_rem_naive(R1, &d1, P, Q);
    fmpz_poly_pseudo_rem(R2, &d2, P, Q);

    // afficher
    printf("Naive remainder:\n");
    fmpz_poly_print(R1);
    printf("\nd1 = %lu\n\n", d1);

    printf("FLINT remainder:\n");
    fmpz_poly_print(R2);
    printf("\nd2 = %lu\n\n", d2);

    // comparaison
    if (fmpz_poly_equal(R1, R2) && d1 == d2)
        printf("OK : meme resultat\n");
    else
        printf("DIFFERENT\n");

    fmpz_poly_clear(R1);
    fmpz_poly_clear(R2);*/

    return 0;
}