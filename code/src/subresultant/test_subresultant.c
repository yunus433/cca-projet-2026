// gcc-15 -O2 subresultant_euclid.h subresultant_naive.h subresultant_pseudo_remainder.h test_subresultant.c  $(pkg-config --cflags --libs flint)

// As FLINT does not expose any functions to compute the subresultants directly,
// this code generates random polynomials and tests if all 3 implementations of the subresultants
// give the same result.

// The naive algorithm is supposed to be the "testing" algorithm.

#include <time.h>
#include <stdlib.h>

#include "subresultant_euclid.c"
#include "subresultant_naive.c"
#include "subresultant_pseudo_remainder.c"

int COEFF_BIT_SIZE = 30;
int DEFAULT_DEGREE = 20;
int DEFAULT_COUNT = 10;

int fmpz_poly_subresultant_test(
  int degree,
  int count
) {
  flint_rand_t rand_state;
  fmpz_t R[degree], S[degree], T[degree];
  fmpz_poly_t A, B;

  flint_randinit(rand_state);
  flint_randseed(rand_state, time(NULL), time(NULL));

  for (int i = 0; i < degree; i++) {
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

    for (int i = 0; i < degree; i++) {
      fmpz_abs(R[i], R[i]);
      fmpz_abs(S[i], S[i]);
      fmpz_abs(T[i], T[i]);

      if (!fmpz_equal(R[i], S[i]) || !fmpz_equal(R[i], T[i])) {
        is_all_successful = 0;

        printf("Test %d is unsuccessful.\n",i+1);
        printf("Polynomial A:\n");
        fmpz_poly_print_pretty(A, "x");
        printf("\nPolynomial B:\n");
        fmpz_poly_print_pretty(B, "x");
        printf("\nR[%d] (naive):\n", i);
        fmpz_print(R[i]);
        printf("\nS[%d] (euclid):\n", i);
        fmpz_print(S[i]);
         printf("\nT[%d] (pseudo remainder):\n", i);
        fmpz_print(T[i]);
        printf("\n");
        break;
      }
    }

    if (is_all_successful)
      printf("Test %d is successful.\n",i+1);
  }

  flint_randclear(rand_state);

  for (int i = 0; i < degree; i++) {
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