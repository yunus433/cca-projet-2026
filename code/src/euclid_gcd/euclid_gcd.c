#include "euclid_gcd.h"

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

  while (!fmpz_poly_is_zero(b)) {
    fmpz_poly_pseudo_rem(tmp, &d, a, b);

    if (!fmpz_poly_is_zero(tmp)){
      fmpz_poly_primitive_part(tmp, tmp);
    }

    fmpz_poly_set(a, b);
    fmpz_poly_set(b, tmp);
  }

  fmpz_poly_set(r, a);

  fmpz_poly_clear(a);
  fmpz_poly_clear(b);
  fmpz_poly_clear(tmp);

  return 0;
}
