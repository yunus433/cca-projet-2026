#ifndef PSEUDOREM_H
#define PSEUDOREM_H

#include <stdio.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

void fmpz_poly_pseudo_rem_naive(fmpz_poly_t R, ulong *d,const fmpz_poly_t A, const fmpz_poly_t B);
#endif