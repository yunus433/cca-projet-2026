#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>

#include "sylvester.h"

int fmpz_poly_subresultants_naive(
  fmpz_t * subresultants,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
) {
  fmpz_mat_t S, W;
  fmpz_poly_sylvester_matrix(S, P, Q);

  int n = fmpz_poly_degree(P);
  int m = fmpz_poly_degree(Q);
  int N = n + m, temp = n;

  if (m > n) { // m is always the smaller degree
    n = m;
    m = temp;
  }

  fmpz_mat_window_init(W, S, 0, 0, N, N); // The entire matrix

  int K = m;

  for (int k = 0; k <= K; k++) {

    fmpz_mat_det(subresultants[k], W); // The matrix W already represents the new submatrix

    for (int i = m - 1; i < m + n - 1; i++) {
      fmpz_mat_swap_cols(W, NULL, i, i + 1);
    }

    n -= 1;
    m -= 1;
    N -= 2; // The matrix is shrink by 2 columns / rows at each iteration

    // TODO regarder s'il y a window_clear
    fmpz_mat_window_init(W, W, 0, 0, N, N); // The new submatrix 
  }

  return 1;
}

/*int main () {
  fmpz_poly_t f, g;
  fmpz_poly_init(f);
  fmpz_poly_init(g);

  // Example: f = 824x^5 - 65x^4 -814x^3 - 741x^2 - 979x- 764, g = 216x^4 + 663x^3 + 880x^2 + 916x + 617
  fmpz_poly_set_coeff_si(g, 0, -5);
  fmpz_poly_set_coeff_si(g, 1, 2);
  fmpz_poly_set_coeff_si(g, 2, 8);
  fmpz_poly_set_coeff_si(g, 3, -3);
  fmpz_poly_set_coeff_si(g, 4, -3);
  fmpz_poly_set_coeff_si(g, 6, 1);
  fmpz_poly_set_coeff_si(g, 8, 1);

  fmpz_poly_set_coeff_si(f, 0, 21);
  fmpz_poly_set_coeff_si(f, 1, -9);
  fmpz_poly_set_coeff_si(f, 2, -4);
  fmpz_poly_set_coeff_si(f, 4, 5);
  fmpz_poly_set_coeff_si(f, 6, 3);

  fmpz_t subresultants[5];

  for (int i = 0; i < 5; i++)
    fmpz_init(subresultants[i]);

  fmpz_poly_subresultants_naive(subresultants, f, g);
  for (int i = 0; i < 5; i++) {
    fmpz_print(subresultants[i]);
    printf("\n");
  }

  return 0;
}*/