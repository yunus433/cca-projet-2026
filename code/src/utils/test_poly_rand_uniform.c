// gcc-15 -O2 test_poly_rand_uniform.c utils.h $(pkg-config --cflags --libs flint)

#include <time.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

#include "poly_rand_uniform.c"

#define COEFF_BIT_SIZE 30
#define TEST_CASE_COUNT 10
#define LENGTH 10

int main() {
  flint_rand_t rand_state;
  fmpz_poly_t P;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));

  fmpz_poly_init(P);

  for (int i = 0; i < TEST_CASE_COUNT; i++) {
    fmpz_poly_randtest_uniform(P, rand_state, LENGTH, COEFF_BIT_SIZE);
    fmpz_poly_print_pretty(P, "x");
    printf("\n");
  }

  flint_rand_clear(rand_state);
  fmpz_poly_clear(P);
}