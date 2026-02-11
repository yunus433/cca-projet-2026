#ifndef SUBRES2_H
#define SUB_RES2_H

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_mat.h>
#include <flint/nmod_poly.h>

int fmpq_poly_subres_euclid(fmpq_t *subres,const fmpq_poly_t P, const fmpq_poly_t Q);
#endif