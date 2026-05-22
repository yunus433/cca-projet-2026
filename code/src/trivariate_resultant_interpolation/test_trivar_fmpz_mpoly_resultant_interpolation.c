// gcc-15 -O2 test_trivar_fmpz_mpoly_resultant_interpolation.c trivar_fmpz_mpoly_resultant_interpolation.c $(pkg-config --cflags --libs flint)

#include "trivar_fmpz_mpoly_resultant_interpolation.h"

int COEFF_BIT_SIZE = 8;
int DEFAULT_EXPO_BIT_SIZE = 4;
int DEFAULT_LENGTH = 10;
int DEFAULT_COUNT = 1;

int trivar_fmpz_mpoly_resultant_interpolation_test(
  int expo_bit_size,
  int length,
  int count
) {
  const char *vars[] = {"x", "y", "z"};
  flint_rand_t rand_state;
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_t P, Q, R1, R2;
  
  clock_t start, end;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpz_mpoly_ctx_init(ctx, 3, ORD_LEX); // 3 for trivariate
  fmpz_mpoly_init(P, ctx);
  fmpz_mpoly_init(Q, ctx);
  fmpz_mpoly_init(R1, ctx);
  fmpz_mpoly_init(R2, ctx);

  int var_to_compute = 0, var_to_evaluate = 1;

  for (int test_case = 0; test_case < count; test_case++) {
    fmpz_mpoly_randtest_bits(P, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);
    fmpz_mpoly_randtest_bits(Q, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);

    fmpz_mpoly_print_pretty(P, vars, ctx);
    printf("\n\n");
    fmpz_mpoly_print_pretty(Q, vars, ctx);
    printf("\n\n");

    start = clock();
    if (trivar_fmpz_mpoly_resultant_interpolation(R1, P, Q, var_to_compute, ctx)) {
      printf("ERROR: trivar_fmpz_mpoly_resultant_interpolation function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Interpolation: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    if (!fmpz_mpoly_resultant(R2, P, Q, var_to_compute, ctx)) {
      printf("ERROR: fmpz_mpoly_resultant function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Naive: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    printf("Real Degree: %ld\n", fmpz_mpoly_degree_si(R2, var_to_evaluate, ctx));

    if (fmpz_mpoly_equal(R1, R2, ctx))
      printf("Test %d is successful.\n", test_case+1);
    else
      printf("Test %d is failed.\n", test_case+1);

    // fmpz_mpoly_print_pretty(R1, vars, ctx);
    // printf("\n\n");
    // fmpz_mpoly_print_pretty(R2, vars, ctx);
    // printf("\n\n");
  }

cleanup:
  flint_rand_clear(rand_state);
  fmpz_mpoly_clear(P, ctx);
  fmpz_mpoly_clear(Q, ctx);
  fmpz_mpoly_clear(R1, ctx);
  fmpz_mpoly_clear(R2, ctx);
  fmpz_mpoly_ctx_clear(ctx);
  return 0;
}


int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("ERROR: Please specify the mode.\n");
    return -1;
  }

  char * mode = argv[1];

  int expo_bit_size = DEFAULT_EXPO_BIT_SIZE;
  int length = DEFAULT_LENGTH;
  int count = DEFAULT_COUNT;

  if (argc > 2) expo_bit_size = atoi(argv[2]);
  if (argc > 3) length = atoi(argv[3]);
  if (argc > 4) count = atoi(argv[4]);

  if (!strcmp("fmpz", mode)) {
    return trivar_fmpz_mpoly_resultant_interpolation_test(expo_bit_size, length, count);
  } else {
    printf("ERROR: Unknown mode specified.\n");
    return -1;
  }
}