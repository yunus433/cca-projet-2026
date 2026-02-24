#include "euclidean_division.h"

int fmpq_poly_euclidean_division(
  fmpq_poly_t q,
  fmpq_poly_t r,
  const fmpq_poly_t a,
  const fmpq_poly_t b
) {
  if (fmpq_poly_is_zero(b))
    return -1;

  fmpq_poly_zero(q);
  fmpq_poly_set(r, a);

  slong d, dr;
  d = fmpq_poly_degree(b);

  fmpq_t c, lc, coeff;
  fmpq_init(c);
  fmpq_init(lc);
  fmpq_init(coeff);

  fmpq_poly_t s;
  fmpq_poly_init(s);

  fmpq_poly_get_coeff_fmpq(c,b,d);

  while (fmpq_poly_degree(r) >= d && !fmpq_poly_is_zero(r)) {
    dr = fmpq_poly_degree(r);

    fmpq_poly_get_coeff_fmpq(lc,r,dr);
    fmpq_div(coeff,lc,c);
    
    fmpq_poly_set_coeff_fmpq(s, dr - d,  coeff);

    fmpq_poly_add(q, q, s);
    fmpq_poly_mul(s, s, b);
    fmpq_poly_sub(r, r, s);
  }

  fmpq_clear(c);
  fmpq_clear(lc);
  fmpq_clear(coeff);
  fmpq_poly_clear(s);

  return 0;
}
