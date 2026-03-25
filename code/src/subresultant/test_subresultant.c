// gcc-15 -O2 ../sylvester_matrix/sylvester_matrix.c subresultant_euclid.c subresultant_naive.c subresultant_pseudo_remainder.c test_subresultant.c  $(pkg-config --cflags --libs flint)

// As FLINT does not expose any functions to compute the subresultants directly,
// this code generates random polynomials and tests if all 3 implementations of the subresultants
// give the same result.

// The naive algorithm is supposed to be the "testing" algorithm.

#include <time.h>
#include <stdlib.h>

#include "subresultant_euclid.h"
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
  fmpz_t R[degree+1], S[degree+1], T[degree+1];
  fmpz_poly_t A, B;

  flint_randinit(rand_state);
  flint_randseed(rand_state, time(NULL), time(NULL));

  for (int i = 0; i <= degree; i++) {
    fmpz_init(R[i]);
    fmpz_init(S[i]);
    fmpz_init(T[i]);
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
    
    if (fmpz_poly_subresultant_naive(R, A, B)) {
      printf("ERROR: subresultant_naive function call returned non-zero error code.\n");
      printf("Polynomial A:\n");
      fmpz_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpz_poly_print_pretty(B, "x");
      printf("\n");
      break;
    }

    if (fmpz_poly_subresultant_naive(S, A, B)) {
      printf("ERROR: subresultant_euclid function call returned non-zero error code.\n");
      printf("Polynomial A:\n");
      fmpz_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpz_poly_print_pretty(B, "x");
      printf("\n");
      break;
    }

    if (fmpz_poly_subresultant_pseudo_remainder(T, A, B)) {
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
      fmpz_abs(T[j], T[j]);

      if (!fmpz_equal(R[j], S[j]) || !fmpz_equal(R[j], T[j])) {
        is_all_successful = 0;

        printf("Test %d is unsuccessful.\n",j+1);
        printf("Polynomial A:\n");
        fmpz_poly_print_pretty(A, "x");
        printf("\nPolynomial B:\n");
        fmpz_poly_print_pretty(B, "x");
        printf("\nR[%d] (naive):\n", j);
        fmpz_print(R[j]);
        printf("\nS[%d] (euclid):\n", j);
        fmpz_print(S[j]);
         printf("\nT[%d] (pseudo remainder):\n", j);
        fmpz_print(T[j]);
        printf("\n");
        break;
      }
    }

    if (is_all_successful)
      printf("Test %d is successful.\n",i+1);
      for (int k = 0; k <= degree; k++) {
        fmpz_zero(R[k]);
        fmpz_zero(S[k]);
        fmpz_zero(T[k]);
      }
  }

  flint_randclear(rand_state);

  for (int i = 0; i <= degree; i++) {
    fmpz_clear(R[i]);
    fmpz_clear(S[i]);
    fmpz_clear(T[i]);
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