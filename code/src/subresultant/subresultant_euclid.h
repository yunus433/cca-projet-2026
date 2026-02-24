#ifndef SUBRESULTANT_EUCLID_H
#define SUBRESULTANT_EUCLID_H

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_mat.h>
#include <flint/nmod_poly.h>

int fmpq_poly_subresultant_euclid(
  fmpq_t *subresultants,
  const fmpq_poly_t P,
  const fmpq_poly_t Q
);

#endif