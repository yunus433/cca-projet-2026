// gcc-15 -O2 euclidean_division.c test_euclidean_division.c  $(pkg-config --cflags --libs flint)

#include <time.h>

#include "euclidean_division.h"

int COEFF_BIT_SIZE = 5;

int fmpq_poly_euclidean_division_test(
  int degree,
  int count
) {
  flint_rand_t rand_state;
  fmpq_poly_t Q, R, S, A, B;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));

  fmpq_poly_init(Q);
  fmpq_poly_init(R);
  fmpq_poly_init(S);
  fmpq_poly_init(A);
  fmpq_poly_init(B);

  int length = degree + 1;

  for (int i = 0; i < count; i++) {
    printf("Starting test %d...\n", i+1);
    
    fmpq_poly_randtest(A, rand_state, length + 1, COEFF_BIT_SIZE);
    fmpq_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);

    while (fmpq_poly_is_zero(B)) fmpq_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);
    
    if (fmpq_poly_euclidean_division(Q, R, A, B)) {
      printf("ERROR: Function call returned non-zero error code.\n");
      printf("Polynomial A:\n");
      fmpq_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpq_poly_print_pretty(B, "x");
      printf("\n");
      break;
    }

    fmpq_poly_mul(S, Q, B);
    fmpq_poly_add(S, S, R);
    
    if (!fmpq_poly_equal(S, A)) {
      printf("Test %d is unsuccessful.\n",i+1);
      printf("Polynomial A:\n");
      fmpq_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpq_poly_print_pretty(B, "x");
      printf("\nPolynomial Q:\n");
      fmpq_poly_print_pretty(Q, "x");
      printf("\nPolynomial R:\n");
      fmpq_poly_print_pretty(R, "x");
      printf("\nPolynomial S = Q * B + R:\n");
      fmpq_poly_print_pretty(S, "x");
      printf("\n");
      continue;
    }

    printf("Test %d is successful.\n",i+1);
  }

  flint_rand_clear(rand_state);
  fmpq_poly_clear(Q);
  fmpq_poly_clear(R);
  fmpq_poly_clear(S);
  fmpq_poly_clear(A);
  fmpq_poly_clear(B);
  return 0;
}

int main() {
  return fmpq_poly_euclidean_division_test(12, 10);
}