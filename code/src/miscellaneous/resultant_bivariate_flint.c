// gcc-15 -O2 resultant_bivariate_flint.c $(pkg-config --cflags --libs flint)

// This file is to test & understand how FLINT resultant works for bivariate polynomials.

#include <time.h>
#include <stdlib.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h> // Not poly, mpoly!

int COEFF_BIT_SIZE = 30;
int EXPO_BIT_SIZE = 5;
int DEFAULT_LENGTH = 20;
int DEFAULT_COUNT = 10;
int ELIMINATE_VAR = 1; // 0 for X

int test_resultant_bivariate_flint(
  int length,
  int count
) {
  const char *vars[] = {"x", "y"};

  flint_rand_t rand_state;
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_t P, Q, R;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpz_mpoly_init(P, ctx);
  fmpz_mpoly_init(Q, ctx);
  fmpz_mpoly_init(R, ctx);

  for (int i = 0; i < count; i++) {
    printf("Starting test %d...\n", i+1);

    fmpz_mpoly_randtest_bits(P, rand_state, length, COEFF_BIT_SIZE, EXPO_BIT_SIZE, ctx);
    fmpz_mpoly_randtest_bits(Q, rand_state, length, COEFF_BIT_SIZE, EXPO_BIT_SIZE, ctx);

    fmpz_mpoly_print_pretty(P, vars, ctx);
    printf("\n");
    fmpz_mpoly_print_pretty(Q, vars, ctx);
    printf("\n");

    if (!fmpz_mpoly_resultant(R, P, Q, ELIMINATE_VAR, ctx)) {
      printf("ERROR: fmpz_mpoly_resultant function call returned non-zero error code.\n");
      break;
    }

    fmpz_mpoly_print_pretty(R, vars, ctx);
    printf("\n");

    printf("Test %d is successful.\n", i+1);
  }

  flint_rand_clear(rand_state);
  fmpz_mpoly_clear(P, ctx);
  fmpz_mpoly_clear(Q, ctx);
  fmpz_mpoly_clear(R, ctx);
  fmpz_mpoly_ctx_clear(ctx);

  return 0;
}

int main(int argc, char ** argv) {
  int length = DEFAULT_LENGTH;
  int count = DEFAULT_COUNT;

  if (argc > 1) length = atoi(argv[1]);
  if (argc > 2) count = atoi(argv[2]);

  return test_resultant_bivariate_flint(length, count);
}