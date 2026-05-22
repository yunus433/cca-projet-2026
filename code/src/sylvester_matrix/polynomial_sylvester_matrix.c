#include "polynomial_sylvester_matrix.h"

int fmpz_poly_polynomial_sylvester_matrix(
  fmpz_poly_mat_t S,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
) {
  fmpz_t tmp; // Temp values
  fmpz_init(tmp);

  int n = fmpz_poly_degree(P);
  int m = fmpz_poly_degree(Q);

  if (m < 0 || n < 0) { // If a polynomial does not exist (degree -1 in Flint)
    fmpz_clear(tmp);
    return -1;
  }

  fmpz_poly_mat_zero(S);

  if (m <= n) {
    for (int i = n; i >= 0; i--) {
      fmpz_poly_get_coeff_fmpz(tmp, P, i); // The coefficient of poly

      for (int j = 0; j < m; j++) { // We shift m times
        fmpz_poly_set_coeff_fmpz(fmpz_poly_mat_entry(S, j + (n - i), j), 0, tmp);
      }
    }
    
    for (int i = m; i >= 0; i--) {
      fmpz_poly_get_coeff_fmpz(tmp, Q, i); // The coefficient of poly

      for (int j = 0; j < n; j++) { // We shift n times
        fmpz_poly_set_coeff_fmpz(fmpz_poly_mat_entry(S, j + (m - i), j + m), 0, tmp);
      }
    }
  } else {
    for (int i = m; i >= 0; i--) {
      fmpz_poly_get_coeff_fmpz(tmp, Q, i); // The coefficient of poly

      for (int j = 0; j < n; j++) { // We shift n times
        fmpz_poly_set_coeff_fmpz(fmpz_poly_mat_entry(S, j + (m - i), j), 0, tmp);
      }
    }

    for (int i = n; i >= 0; i--) {
      fmpz_poly_get_coeff_fmpz(tmp, P, i); // The coefficient of poly

      for (int j = 0; j < m; j++) { // We shift m times
        fmpz_poly_set_coeff_fmpz(fmpz_poly_mat_entry(S, j + (n - i), j + n), 0, tmp);
      }
    }
  }

  fmpz_clear(tmp);
  return 0;
}
