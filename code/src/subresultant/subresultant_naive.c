#include "subresultant_naive.h"

int fmpz_poly_subresultant_naive(
  fmpz *subresultants,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
) {
  int n = fmpz_poly_degree(P);
  int m = fmpz_poly_degree(Q);

  if (n < 0 || m < 0) {
    return 1;
  }

  if (m > n) {
    int tmp = n;
    n = m;
    m = tmp;
  }

  int N = n + m;
  int K = m;

  fmpz_mat_t S, W;
  fmpz_mat_init(S, N, N);
  fmpz_poly_sylvester_matrix(S, P, Q);
  fmpz_mat_window_init(W, S, 0, 0, N, N);

  for (int k = 0; k <= K; k++) {
    fmpz_mat_det(subresultants + k, W);

    if (m > 0) {
      for (int i = m - 1; i < m + n - 1; i++) {
        if (i + 1 >= N) return 2;
        fmpz_mat_swap_cols(W, NULL, i, i + 1);
      }
    }

    n--;
    m--;
    N -= 2;

    if (k < K && N > 0) {
      fmpz_mat_window_clear(W);
      fmpz_mat_window_init(W, S, 0, 0, N, N);
    }
  }

  fmpz_mat_window_clear(W);
  fmpz_mat_clear(S);
  return 0;
}