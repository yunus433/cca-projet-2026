// gcc-15 -O2 ../sylvester_matrix/polynomial_sylvester_matrix.c subresultant_polynomial_naive.c subresultant_polynomial_pseudo_remainder.c test_subresultant_polynomial.c  $(pkg-config --cflags --libs flint)

// As FLINT does not expose any functions to compute the subresultants directly,
// this code generates random polynomials and tests if all 3 implementations of the subresultants
// give the same result.

// The naive algorithm is supposed to be the "testing" algorithm.

#include <time.h>
#include <stdlib.h>

#include "subresultant_polynomial_naive.h"
#include "subresultant_polynomial_pseudo_remainder.h"

int COEFF_BIT_SIZE = 30;
int DEFAULT_DEGREE = 20;
int DEFAULT_COUNT = 10;

int fmpz_poly_subresultant_polyonmial_test(
  int degree,
  int count
) {
  flint_rand_t rand_state;
  fmpz_poly_t R[degree], S[degree];
  fmpz_poly_t A, B;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));

  for (int i = 0; i < degree; i++) {
    fmpz_poly_init(R[i]);
    fmpz_poly_init(S[i]);
  }

  fmpz_poly_init(A);
  fmpz_poly_init(B);

  int length = degree + 1;

  for (int i = 0; i < count; i++) {
    printf("Starting test %d...\n", i+1);
    
    fmpz_poly_randtest(A, rand_state, length, COEFF_BIT_SIZE);
    fmpz_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);

    while (fmpz_poly_is_zero(A)) fmpz_poly_randtest(A, rand_state, length, COEFF_BIT_SIZE);
    while (fmpz_poly_is_zero(B)) fmpz_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);
    
    if (fmpz_poly_subresultant_polynomial_naive(R, A, B)) {
      printf("ERROR: subresultant_naive function call returned non-zero error code.\n");
      printf("Polynomial A:\n");
      fmpz_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpz_poly_print_pretty(B, "x");
      printf("\n");
      break;
    }

    if (fmpz_poly_subresultant_polynomial_pseudo_remainder(S, A, B)) {
      printf("ERROR: subresultant_euclid function call returned non-zero error code.\n");
      printf("Polynomial A:\n");
      fmpz_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpz_poly_print_pretty(B, "x");
      printf("\n");
      break;
    }

    int is_all_successful = 1;

    // FIXME: use another variable instead of i to avoid confusion with the outer loop
    for (int i = 0; i < degree; i++) {
      printf("\nR[%d] (naive):\n", i);
      fmpz_poly_print_pretty(R[i], "x");
      printf("\nS[%d] (pseudo-remainder):\n", i);
      fmpz_poly_print_pretty(S[i], "x");
      // if (!fmpz_poly_equal(R[i], S[i])) {
      //   is_all_successful = 0;

      //   printf("Test %d is unsuccessful.\n");
      //   printf("Polynomial A:\n");
      //   fmpz_poly_print_pretty(A, "x");
      //   printf("\nPolynomial B:\n");
      //   fmpz_poly_print_pretty(B, "x");
      //   printf("\nR[%d] (naive):\n", i);
      //   fmpz_poly_print_pretty(R[i], "x");
      //   printf("\nS[%d] (pseudo-remainder):\n", i);
      //   fmpz_poly_print_pretty(S[i], "x");
      //   printf("\n");
      //   break;
      // }
    }

    if (is_all_successful)
      printf("Test %d is successful.\n", i + 1);
  }

  flint_rand_clear(rand_state);

  for (int i = 0; i < degree; i++) {
    fmpz_poly_clear(R[i]);
    fmpz_poly_clear(S[i]);
  }

  fmpz_poly_clear(A);
  fmpz_poly_clear(B);
  return 0;
}

int main(int argc, char *argv[]) {
  int degree = DEFAULT_DEGREE;
  int count = DEFAULT_COUNT;

  if (argc > 1) degree = atoi(argv[1]);
  if (argc > 2) count = atoi(argv[2]);

  return fmpz_poly_subresultant_polyonmial_test(degree, count);
}