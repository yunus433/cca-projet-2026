#ifndef SYLV_H
#define SYLV_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

int fmpz_poly_sylvester_matrix(
  fmpz_mat_t S, // Not a pointer, already accessed by referenced in Sylvester types. No dynamic memory allocation either, so no need to free later on
  const fmpz_poly_t P,
  const fmpz_poly_t Q
) {
  fmpz_t tmp; // Temp values

  int n = fmpz_poly_degree(P);
  int m = fmpz_poly_degree(Q);

  int N = n + m; // The Sylvester matrix has size (n+m)x(n+m)

  if (m < 0 || n < 0) { // If a polynomial does not exist (degree -1 in Flint)
    fmpz_mat_init(S, 0, 0); // A matrix of 0 size?
    return -1;
  }

  fmpz_mat_init(S, N, N);
  fmpz_mat_zero(S);

  if (m <= n) {
    for (int i = n; i >= 0; i--) {
      fmpz_poly_get_coeff_fmpz(tmp, P, i); // The coefficient of poly

      for (int j = 0; j < m; j++) { // We shift m times
        fmpz_set(fmpz_mat_entry(S, j + (n - i), j), tmp);
      }
    }
    
    for (int i = m; i >= 0; i--) {
      fmpz_poly_get_coeff_fmpz(tmp, Q, i); // The coefficient of poly

      for (int j = 0; j < n; j++) { // We shift n times
        fmpz_set(fmpz_mat_entry(S, j + (m - i), j + m), tmp);
      }
    }
  } else {
    for (int i = m; i >= 0; i--) {
      fmpz_poly_get_coeff_fmpz(tmp, Q, i); // The coefficient of poly

      for (int j = 0; j < n; j++) { // We shift n times
        fmpz_set(fmpz_mat_entry(S, j + (m - i), j), tmp);
      }
    }

    for (int i = n; i >= 0; i--) {
      fmpz_poly_get_coeff_fmpz(tmp, P, i); // The coefficient of poly

      for (int j = 0; j < m; j++) { // We shift m times
        fmpz_set(fmpz_mat_entry(S, j + (n - i), j + n), tmp);
      }
    }
  }

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

//   fmpz_mat_t S;

//   fmpz_poly_sylvester_matrix(S, f, g);

//   fmpz_mat_print(S);

//   return 0;
// }