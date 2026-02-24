#include "subresultant_naive.h"

int fmpz_poly_subresultant_naive(
  fmpz_t * subresultants,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
) {
  int n = fmpz_poly_degree(P);
  int m = fmpz_poly_degree(Q);
  int N = n + m, tmp = n;

  if (m > n) { // m is always the smaller degree
    n = m;
    m = tmp;
  }

  fmpz_mat_t S, W;
  fmpz_mat_init(S, N, N);
  fmpz_poly_sylvester_matrix(S, P, Q);
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

  fmpz_mat_clear(S);
  fmpz_mat_window_clear(W);

  return 0;
}
