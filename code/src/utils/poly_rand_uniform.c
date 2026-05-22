// This file generates a random fmpz polynomial with uniform coefficients, meaning each coefficient has equal chances of getting 0.
// In normal fmpz_poly_randtest(), number of terms is also random, thus some polynomials have much smaller number of terms,
// even though the degree can be huge. This results in sparse polynomials and less efficient testing with subresultants.

// The function has the exact same signature as fmpz_poly_randtest().

#include "utils.h"

void fmpz_poly_randtest_uniform(
  fmpz_poly_struct *f,
  flint_rand_struct *state,
  slong len,
  flint_bitcnt_t bits
) {
  fmpz_t coeff;
  fmpz_init(coeff);

  fmpz_poly_zero(f); // In case it was preset

  for (int i = 0; i < len; i++) {
    fmpz_randbits(coeff, state, bits);
    fmpz_poly_set_coeff_fmpz(f, i, coeff);
  }

  fmpz_clear(coeff);
};