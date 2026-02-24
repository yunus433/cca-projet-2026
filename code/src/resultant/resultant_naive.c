#include "resultant_naive.h"

int fmpz_poly_resultant_naive(
  fmpz_t resultant,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
) {
  int n = fmpz_poly_degree(P);
  int m = fmpz_poly_degree(Q);
  int N = n + m;

  fmpz_mat_t S;
  fmpz_mat_init(S, N, N);

  fmpz_poly_sylvester_matrix(S, P, Q);

  fmpz_mat_det(resultant, S);

  fmpz_mat_clear(S);

  return 0;
}
