// gcc-15 -O2 ../sylvester_matrix/sylvester_matrix.c subresultant_naive.c subresultant_pseudo_remainder.c test_subresultant.c  $(pkg-config --cflags --libs flint)

// As FLINT does not expose any functions to compute the subresultants directly,
// this code generates random polynomials and tests if both implementations of the subresultants
// give the same result.

// If you want to test subresultant codes against pre-generated test cases, use test_cases_subresultant.c.

// Note that since subresultant_euclid is on Q, it is not tested. The implementation is only a reference.

// The naive algorithm is supposed to be the "testing" algorithm.

// How it works?
// ./a.out {DEGREE} {COUNT}

// DEGREE := Degree of polynomials tested.
// COUNT := Number of tests to be performed.

// Note: The coefficient bit size is set to be COEFF_BIT_SIZE as a constant, you can change it in the code.

#include <time.h>
#include <stdlib.h>

#include "subresultant_naive.h"
#include "subresultant_pseudo_remainder.h"

int COEFF_BIT_SIZE = 30;
int DEFAULT_DEGREE = 20;
int DEFAULT_COUNT = 10;

int fmpz_poly_subresultant_test(
  int degree,
  int count
) {
  flint_rand_t rand_state;
  fmpz_t R[degree+1], S[degree+1];
  fmpz_poly_t A, B;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));

  for (int i = 0; i <= degree; i++) {
    fmpz_init(R[i]);
    fmpz_init(S[i]);
  }

  fmpz_poly_init(A);
  fmpz_poly_init(B);

  int length = degree + 1;

  for (int i = 0; i < count; i++) {
    printf("Starting test %d...\n", i+1);
    
    fmpz_poly_randtest(A, rand_state, length, COEFF_BIT_SIZE);
    fmpz_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);

    // TO DO: NOT USE THIS
    while (fmpz_poly_is_zero(A)) fmpz_poly_randtest(A, rand_state, length, COEFF_BIT_SIZE);
    while (fmpz_poly_is_zero(B)) fmpz_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);
    
    if (fmpz_poly_subresultant_naive(R, A, B)) {
      printf("ERROR: subresultant_naive function call returned non-zero error code.\n");
      printf("Polynomial A:\n");
      fmpz_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpz_poly_print_pretty(B, "x");
      printf("\n");
      break;
    }

    if (fmpz_poly_subresultant_pseudo_remainder(S, A, B)) {
      printf("ERROR: subresultant_naive function call returned non-zero error code.\n");
      printf("Polynomial A:\n");
      fmpz_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpz_poly_print_pretty(B, "x");
      printf("\n");
      break;
    }

    int is_all_successful = 1;

    // FIXME: use another variable instead of i here to avoid confusion with the outer loop
    for (int j = 0; j < degree; j++) {
      fmpz_abs(R[j], R[j]);
      fmpz_abs(S[j], S[j]);

      if (!fmpz_equal(R[j], S[j])) {
        is_all_successful = 0;

        break;
      }
    }

    if (!is_all_successful) {
      printf("Test %d is unsuccessful.\n", i+1);
      printf("Polynomial A:\n");
      fmpz_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpz_poly_print_pretty(B, "x");
      for (int j = 0; j < degree; j++) {
        printf("\nR[%d] (naive):\n", j);
        fmpz_print(R[j]);
         printf("\nS[%d] (pseudo remainder):\n", j);
        fmpz_print(S[j]);
        printf("\n");
      }
    }

    if (is_all_successful)
      printf("Test %d is successful.\n", i+1);
      for (int k = 0; k <= degree; k++) {
        fmpz_zero(R[k]);
        fmpz_zero(S[k]);
      }
  }

  flint_rand_clear(rand_state);

  for (int i = 0; i <= degree; i++) {
    fmpz_clear(R[i]);
    fmpz_clear(S[i]);
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

  return fmpz_poly_subresultant_test(degree, count);
}