// gcc-15 -O2 resultant_naive.h test_resultant_naive.c  $(pkg-config --cflags --libs flint)

#include <time.h>
#include <stdlib.h>

#include "resultant_naive.c"

int COEFF_BIT_SIZE = 30;
int DEFAULT_DEGREE = 20;
int DEFAULT_COUNT = 10;

int fmpz_poly_resultant_naive_test(
  int degree,
  int count
) {
  flint_rand_t rand_state;
  fmpz_t R, S;
  fmpz_poly_t A, B;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));

  fmpz_init(R);
  fmpz_init(S);
  fmpz_poly_init(A);
  fmpz_poly_init(B);

  int length = degree + 1;

  for (int i = 0; i < count; i++) {
    printf("Starting test %d...\n", i+1);
    
    fmpz_poly_randtest(A, rand_state, length, COEFF_BIT_SIZE);
    fmpz_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);

    while (fmpz_poly_is_zero(A)) fmpz_poly_randtest(A, rand_state, length, COEFF_BIT_SIZE);
    while (fmpz_poly_is_zero(B)) fmpz_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);
    
    if (fmpz_poly_resultant_naive(R, A, B)) {
      printf("ERROR: Function call returned non-zero error code.\n");
      printf("Polynomial A:\n");
      fmpz_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpz_poly_print_pretty(B, "x");
      printf("\n");
      break;
    }

    fmpz_poly_resultant(S, A, B);
    
    fmpz_abs(R, R);
    fmpz_abs(S, S);

    if (!fmpz_equal(R, S)) {
      printf("Test %d is unsuccessful.\n");
      printf("Polynomial A:\n");
      fmpz_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpz_poly_print_pretty(B, "x");
      printf("\nR (implementation):\n");
      fmpz_print(R);
      printf("\nS (library):\n");
      fmpz_print(S);
      printf("\n");
      continue;
    }

    printf("Test %d is successful.\n");
  }

  flint_rand_clear(rand_state);
  fmpz_clear(R);
  fmpz_clear(S);
  fmpz_poly_clear(A);
  fmpz_poly_clear(B);
}

int main(int argc, char *argv[]) {
  int degree = DEFAULT_DEGREE;
  int count = DEFAULT_COUNT;

  if (argc > 1) degree = atoi(argv[1]);
  if (argc > 2) count = atoi(argv[2]);

  fmpz_poly_resultant_naive_test(degree, count);
}