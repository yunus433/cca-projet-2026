// gcc-15 -O2 test.c $(pkg-config --cflags --libs flint)

// This file is to test polynomial interpolation in FLINT

#include <time.h>
#include <stdlib.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

int COEFF_BIT_SIZE = 30;
int DEFAULT_DEGREE = 20;
int DEFAULT_COUNT = 10;

int test_interpolation_flint(
  int degree,
  int count
) {
  int length = degree + 1;

  flint_rand_t rand_state;
  fmpz_poly_t P, Q;
  fmpz x[length], y[length];

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));

  fmpz_poly_init(P);
  fmpz_poly_init(Q);

  for (int i = 0; i < length; i++) {
    fmpz_init(x + i);
    fmpz_init(y + i);
  }

  for (int i = 0; i < count; i++) {
    printf("Starting test %d...\n", i+1);

    fmpz_poly_randtest(P, rand_state, length, COEFF_BIT_SIZE);
    while (fmpz_poly_is_zero(P)) fmpz_poly_randtest(P, rand_state, length, COEFF_BIT_SIZE);

    for (int j = 0; j < length; j++) {
      fmpz_randbits(x + j, rand_state, COEFF_BIT_SIZE);
      fmpz_poly_evaluate_fmpz(y + j, P, x + j);
    }

    if (!fmpz_poly_interpolate(Q, x, y, length)) {
      printf("ERROR: fmpz_poly_interpolate function call returned non-zero error code.\n");
      break;
    }

    if (!fmpz_poly_equal(P, Q)) {
      printf("Test %d is unsuccessful.\n", i+1);
      printf("Polynomial P:\n");
      fmpz_poly_print_pretty(P, "x");
      printf("\nPolynomial Q (evaluated):\n");
      fmpz_poly_print_pretty(Q, "x");
      break;
    }

    printf("Test %d is successful.\n", i+1);
  }

  flint_rand_clear(rand_state);
  fmpz_poly_clear(P);
  fmpz_poly_clear(Q);
  for (int i = 0; i < length; i++) {
    fmpz_clear(x + i);
    fmpz_clear(y + i);
  }

  return 0;
}

int main(int argc, char *argv[]) {
  int degree = DEFAULT_DEGREE;
  int count = DEFAULT_COUNT;

  if (argc > 1) degree = atoi(argv[1]);
  if (argc > 2) count = atoi(argv[2]);

  return test_interpolation_flint(degree, count);
}