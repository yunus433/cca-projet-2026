// gcc-15 -O2 euclid_gcd.h test_euclid_gcd.c  $(pkg-config --cflags --libs flint)

#include <time.h>
#include <stdlib.h>

#include "euclid_gcd.c"

int COEFF_BIT_SIZE = 30;
int DEFAULT_DEGREE = 20;
int DEFAULT_COUNT = 10;

int fmpz_poly_euclid_gcd_test(
  int degree,
  int count
) {
  flint_rand_t rand_state;
  fmpz_poly_t R, S, A, B;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));

  fmpz_poly_init(R);
  fmpz_poly_init(S);
  fmpz_poly_init(A);
  fmpz_poly_init(B);

  int length = degree + 1;

  for (int i = 0; i < count; i++) {
    printf("Starting test %d...\n", i+1);
    
    fmpz_poly_randtest(A, rand_state, length, COEFF_BIT_SIZE);
    fmpz_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);

    while (fmpz_poly_is_zero(A)) fmpz_poly_randtest(A, rand_state, length, COEFF_BIT_SIZE);
    while (fmpz_poly_is_zero(B)) fmpz_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);
    
    if (fmpz_poly_gcd_euclid(R, A, B)) {
      printf("ERROR: Function call returned non-zero error code.\n");
      printf("Polynomial A:\n");
      fmpz_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpz_poly_print_pretty(B, "x");
      printf("\n");
      break;
    }

    fmpz_poly_gcd(S, A, B);
    fmpz_poly_scalar_mul_fmpz(S, S, fmpz_poly_lead(R)); // Normalization of the monic polynomial
    
    if (!fmpz_poly_equal(R, S)) {
      printf("Test %d is unsuccessful.\n");
      printf("Polynomial A:\n");
      fmpz_poly_print_pretty(A, "x");
      printf("\nPolynomial B:\n");
      fmpz_poly_print_pretty(B, "x");
      printf("\nPolynomial R (implementation):\n");
      fmpz_poly_print_pretty(R, "x");
      printf("\nPolynomial S (library):\n");
      fmpz_poly_print_pretty(S, "x");
      printf("\n");
      continue;
    }

    printf("Test %d is successful.\n");
  }

  flint_rand_clear(rand_state);
  fmpz_poly_clear(R);
  fmpz_poly_clear(S);
  fmpz_poly_clear(A);
  fmpz_poly_clear(B);
  return 0;
}

int main(int argc, char *argv[]) {
  int degree = DEFAULT_DEGREE;
  int count = DEFAULT_COUNT;

  if (argc > 1) degree = atoi(argv[1]);
  if (argc > 2) count = atoi(argv[2]);

  return fmpz_poly_euclid_gcd_test(degree, count);
}