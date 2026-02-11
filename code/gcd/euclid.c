#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

int fmpz_poly_gcd_euclid(
  fmpz_poly_t r,
  const fmpz_poly_t p,
  const fmpz_poly_t q
) {
  if (fmpz_poly_is_zero(p) && fmpz_poly_is_zero(q))
    return -1;

  fmpz_poly_t a, b, tmp;
  ulong d;
  fmpz_poly_init(a);
  fmpz_poly_init(b);
  fmpz_poly_init(tmp);

  fmpz_poly_set(a, p);
  fmpz_poly_set(b, q);
  int i=0;
  while (!fmpz_poly_is_zero(b)&&i<20) {
    fmpz_poly_pseudo_rem(tmp, &d, a, b);
    fmpz_poly_set(a, b);
    fmpz_poly_set(b, tmp);
    i++;
  }

  fmpz_poly_set(r, a);

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

  fmpz_poly_gcd_euclid(r, f, g);

  fmpz_poly_print(r);
  printf("\n");
}*/