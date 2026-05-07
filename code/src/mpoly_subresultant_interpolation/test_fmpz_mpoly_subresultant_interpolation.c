// gcc-15 -O2 test_fmpz_mpoly_subresultant_interpolation.c fmpz_mpoly_subresultant_interpolation.c ../subresultant/subresultant_pseudo_remainder.c $(pkg-config --cflags --libs flint)
// This is no testing, just prints values. Run test_cases instead.

#include <time.h>

#include "fmpz_mpoly_subresultant_interpolation.h"

int COEFF_BIT_SIZE = 100;
int DEFAULT_EXPO_BIT_SIZE = 4;
int DEFAULT_LENGTH = 30;
int DEFAULT_COUNT = 1;

int fmpz_mpoly_subresultant_interpolation_test(
  int expo_bit_size,
  int length,
  int count
) {
  const char *vars[] = {"x", "y"};
  flint_rand_t rand_state;
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_t P, Q;

  clock_t start, end;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpz_mpoly_init(P, ctx);
  fmpz_mpoly_init(Q, ctx);

  int var_to_compute = 0, var_to_evaluate = 1;

  for (int test_case = 0; test_case < count; test_case++) {
    fmpz_mpoly_randtest_bits(P, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);
    fmpz_mpoly_randtest_bits(Q, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);

    slong degP_var = fmpz_mpoly_degree_si(P, var_to_compute, ctx);
    slong degQ_var = fmpz_mpoly_degree_si(Q, var_to_compute, ctx);

    slong K = degP_var > degQ_var ? degP_var : degQ_var;

    fmpz_mpoly_t *subresultants = FLINT_ARRAY_ALLOC(K + 1, fmpz_mpoly_t);
    if (subresultants == NULL) {
      printf("ERROR: subresultant memory could not be allocated.");
      goto cleanup;
    }

    for (slong i = 0; i < K + 1; i++) {
      fmpz_mpoly_init(subresultants[i], ctx);
    }

    start = clock();
    if (fmpz_mpoly_subresultant_interpolation(subresultants, P, Q, var_to_compute, ctx)) {
      printf("ERROR: fmpz_mpoly_resultant_interpolation function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Measure Time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    printf("Test case %d subresultants:\n", test_case+1);
    for (slong i = 0; i < K + 1; i++) {
      fmpz_mpoly_print_pretty(subresultants[i], vars, ctx);
      printf("\n");
      fmpz_mpoly_clear(subresultants[i], ctx);
    }

    flint_free(subresultants);
  }

cleanup:
  flint_rand_clear(rand_state);
  fmpz_mpoly_clear(P, ctx);
  fmpz_mpoly_clear(Q, ctx);
  fmpz_mpoly_ctx_clear(ctx);
  return 0;
}

int main(int argc, char *argv[]) {
  int expo_bit_size = DEFAULT_EXPO_BIT_SIZE;
  int length = DEFAULT_LENGTH;
  int count = DEFAULT_COUNT;

  if (argc > 1) expo_bit_size = atoi(argv[1]);
  if (argc > 2) length = atoi(argv[2]);
  if (argc > 3) count = atoi(argv[3]);
    
  return fmpz_mpoly_subresultant_interpolation_test(expo_bit_size, length, count);
}
