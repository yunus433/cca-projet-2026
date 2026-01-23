#ifndef RESULTANT_H
#define RESULTANT_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

#include "sylvester.h"

int fmpz_poly_resultant(
  fmpz_t resultant,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
) {
  fmpz_mat_t S;

  fmpz_poly_sylvester_matrix(S, P, Q);

  fmpz_mat_det(resultant, S);

  return 1;
}

#endif

// int main () {
//   fmpz_poly_t f, g;
//   fmpz_poly_init(f);
//   fmpz_poly_init(g);

//   /* Example: f = 1 + 2x + 3x^2 + 5x^3 + 7x^4, g = 5 + 7x + 9x^2 */
//   fmpz_poly_set_coeff_si(f, 0, 1);
//   fmpz_poly_set_coeff_si(f, 1, 2);
//   fmpz_poly_set_coeff_si(f, 2, 3);
//   fmpz_poly_set_coeff_si(f, 3, 5);
//   fmpz_poly_set_coeff_si(f, 4, 7);

//   fmpz_poly_set_coeff_si(g, 0, 5);
//   fmpz_poly_set_coeff_si(g, 1, 7);
//   fmpz_poly_set_coeff_si(g, 2, 9);

//   fmpz_t res;

//   fmpz_poly_resultant(res, f, g);

//   fmpz_print(res); // Prints 19359
//   printf("\n");

//   return 0;
// }