#ifndef EUCLIDEAN_H
#define EUCLIDEAN_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_mat.h>

int fmpq_poly_euclidean_division(
  fmpq_poly_t q,
  fmpq_poly_t r,
  const fmpq_poly_t a,
  const fmpq_poly_t b
);

#endif