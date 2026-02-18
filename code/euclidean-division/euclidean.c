
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
) {
  if (fmpq_poly_is_zero(b))
    return -1;

  fmpq_poly_zero(q);
  fmpq_poly_set(r, a);

  slong d;
  d = fmpq_poly_degree(b);

  fmpq_t c;
  fmpq_init(c);
  fmpq_poly_get_coeff_fmpq(c,b,d);

  while (fmpq_poly_degree(r) >= d && !fmpq_poly_is_zero(r)) {
    fmpq_poly_t s;
    fmpq_poly_init(s);
    slong dr = fmpq_poly_degree(r);
    fmpq_t lc,coeff;
    fmpq_init(lc);
    fmpq_init(coeff);
    fmpq_poly_get_coeff_fmpq(lc,r,dr);
    fmpq_div(coeff,lc,c);
    
    fmpq_poly_set_coeff_fmpq(s, dr - d,  coeff);

    fmpq_poly_add(q, q, s);
    fmpq_poly_mul(s, s, b);
    fmpq_poly_sub(r, r, s);

    fmpq_poly_clear(s);
    fmpq_clear(lc);
    fmpq_clear(coeff);
  }

  return 1;
}


/*int main() {
  fmpz_poly_t f, g, q, r;
  fmpz_poly_init(f);
  fmpz_poly_init(g);
  fmpz_poly_init(q);
  fmpz_poly_init(r);

  fmpz_poly_set_coeff_si(f, 0, 6);
  fmpz_poly_set_coeff_si(f, 1, 7);
  fmpz_poly_set_coeff_si(f, 2, 1);

  fmpz_poly_set_coeff_si(g, 0, -6);
  fmpz_poly_set_coeff_si(g, 1, -5);
  fmpz_poly_set_coeff_si(g, 2, 1);

  fmpz_poly_euclidean_division(q, r, f, g);

  fmpz_poly_print(q);
  printf("\n");
  fmpz_poly_print(r);
  printf("\n");
}*/