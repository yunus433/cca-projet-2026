
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

int fmpz_poly_euclidean_division(
  fmpz_poly_t q,
  fmpz_poly_t r,
  const fmpz_poly_t a,
  const fmpz_poly_t b
) {
  if (fmpz_poly_is_zero(b))
    return -1;

  fmpz_poly_zero(q);
  fmpz_poly_set(r, a);

  slong d, c;
  d = fmpz_poly_degree(b);
  c = * fmpz_poly_lead(b);

  while (fmpz_poly_degree(r) >= d) {
    fmpz_poly_t s;
    fmpz_poly_init(s);

    slong lc = * fmpz_poly_lead(r);
    
    fmpz_poly_set_coeff_si(s, fmpz_poly_degree(r) - d,  (lc / c));

    fmpz_poly_add(q, q, s);
    fmpz_poly_mul(s, s, b);
    fmpz_poly_sub(r, r, s);
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