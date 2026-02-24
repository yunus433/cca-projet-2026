// gcc-15 -O2 subres-poly.c sylvester-poly.h $(pkg-config --cflags --libs flint)
// QUESTION: Comparison with wikipedia?

#include "subresultant_polynomial_naive.h"

void fmpz_poly_mat_swap_cols_inplace(fmpz_poly_mat_t A, slong c1, slong c2) {
    slong r, nrows = fmpz_poly_mat_nrows(A);

    if (c1 == c2) return;

    for (r = 0; r < nrows; r++)
    {
        fmpz_poly_swap(
            fmpz_poly_mat_entry(A, r, c1),
            fmpz_poly_mat_entry(A, r, c2)
        );
    }
}

int fmpz_subresultant_polynomials_naive(
  fmpz_poly_t * subresultant_polynomials,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
) {
  fmpz_poly_mat_t S, W;
  fmpz_poly_polynomial_sylvester_matrix(S, P, Q);

  int n = fmpz_poly_degree(P);
  int m = fmpz_poly_degree(Q);
  int N = n + m, temp = n;

  if (m > n) { // m is always the smaller degree
    n = m;
    m = temp;
  }

  fmpz_poly_mat_window_init(W, S, 0, 0, N, N); // The entire matrix

  // We build the matrix V, polynomial multiplication
  fmpz_poly_mat_t V;
  fmpz_poly_mat_init(V, N, N);

  fmpz_poly_mat_one(V); // The initial version is Id

  int K = m;
  int N_ = n + m;

  for (int k = 0; k <= K; k++) {
    slong nk = fmpz_poly_mat_nrows(V);
    slong mk = fmpz_poly_mat_ncols(W);

    fmpz_poly_mat_t T;
    fmpz_poly_mat_init(T, nk, mk);
    fmpz_poly_mat_mul(T, V, W);
    fmpz_poly_mat_print(V, "x");
    fmpz_poly_mat_print(W, "x");
    fmpz_poly_mat_print(T, "x");
    fmpz_poly_mat_det(subresultant_polynomials[k], T);

    for (int i = 0; i < n+m-1; i++) {
      fmpz_poly_set_coeff_si(fmpz_poly_mat_entry(V, i, i), 0, 0);
      fmpz_poly_set_coeff_si(fmpz_poly_mat_entry(V, i+1, i), 0, 1);
    }
    fmpz_poly_set_coeff_si(fmpz_poly_mat_entry(V, n+m-1, n+m-1), 0, 0);
    fmpz_poly_set_coeff_si(fmpz_poly_mat_entry(V, n+m-1, n+m-2), 0, 0);

    fmpz_poly_mat_window_init(V, V, 2, 1, n+m, N);

    for (int i = 0; i <= k+1; i++)
      fmpz_poly_set_coeff_si(fmpz_poly_mat_entry(V, n+m-3, N-2-i), i, 1);

    for (int i = m - 1; i < m + n - 1; i++) {
      fmpz_poly_mat_swap_cols_inplace(W, i, i + 1);
    }

    n -= 1;
    m -= 1;
    N -= 1; // The matrix is shrink by 1 row at each iteration

    fmpz_poly_mat_window_init(W, W, 0, 0, N, n+m); // The new submatrix 
  }

  return 1;
}

int main () {
  fmpz_poly_t f, g;
  fmpz_poly_init(f);
  fmpz_poly_init(g);

  // Example: f = 824x^5 - 65x^4 -814x^3 - 741x^2 - 979x- 764, g = 216x^4 + 663x^3 + 880x^2 + 916x + 617
  // fmpz_poly_set_coeff_si(g, 0, -764);
  // fmpz_poly_set_coeff_si(g, 1, -979);
  // fmpz_poly_set_coeff_si(g, 2, -741);
  // fmpz_poly_set_coeff_si(g, 3, -814);
  // fmpz_poly_set_coeff_si(g, 4, -65);
  // fmpz_poly_set_coeff_si(g, 5, 824);
  fmpz_poly_set_coeff_si(g, 0, -5);
  fmpz_poly_set_coeff_si(g, 1, 2);
  fmpz_poly_set_coeff_si(g, 2, 8);
  fmpz_poly_set_coeff_si(g, 3, -3);
  fmpz_poly_set_coeff_si(g, 4, -3);
  fmpz_poly_set_coeff_si(g, 6, 1);
  fmpz_poly_set_coeff_si(g, 8, 1);

  // fmpz_poly_set_coeff_si(f, 0, 617);
  // fmpz_poly_set_coeff_si(f, 1, 916);
  // fmpz_poly_set_coeff_si(f, 2, 880);
  // fmpz_poly_set_coeff_si(f, 3, 663);
  // fmpz_poly_set_coeff_si(f, 4, 216);
  fmpz_poly_set_coeff_si(f, 0, 21);
  fmpz_poly_set_coeff_si(f, 1, -9);
  fmpz_poly_set_coeff_si(f, 2, -4);
  fmpz_poly_set_coeff_si(f, 4, 5);
  fmpz_poly_set_coeff_si(f, 6, 3);
  fmpz_poly_set_coeff_si(f, 7, 1);
  
  fmpz_poly_t subresultant_polys[7];

  for (int i = 0; i < 7; i++)
    fmpz_poly_init(subresultant_polys[i]);

  fmpz_subresultant_polynomials_naive(subresultant_polys,  f, g);
  for (int i = 0; i < 7; i++) {
    fmpz_poly_print_pretty(subresultant_polys[i], "x");
    printf("\n");
  }

  return 0;
}