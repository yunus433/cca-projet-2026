#ifndef SUBRESULTANT_PSEUDO_REMAINDER
#define SUBRESULTANT_PSEUDO_REMAINDER

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_poly.h>
#include <flint/fmpz_poly.h>

int fmpz_subres_pseudo_rem(
    fmpz_t *subresultants,
    const fmpz_poly_t P,
    const fmpz_poly_t Q
);

#endif