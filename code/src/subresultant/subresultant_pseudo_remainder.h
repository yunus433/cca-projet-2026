#ifndef SUBRESULTANT_PSEUDO_REMAINDER_H
#define SUBRESULTANT_PSEUDO_REMAINDER_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_poly.h>
#include <flint/fmpz_poly.h>

int fmpz_poly_subresultant_pseudo_remainder(
    fmpz_t *subresultants,
    const fmpz_poly_t P,
    const fmpz_poly_t Q
);

#endif