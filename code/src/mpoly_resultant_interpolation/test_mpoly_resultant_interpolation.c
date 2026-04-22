// gcc-15 -O2 test_mpoly_resultant_interpolation.c fmpz_mpoly_resultant_interpolation.c fmpq_mpoly_resultant_interpolation.c nmod_mpoly_resultant_interpolation.c fmpz_mpoly_resultant_interpolation_multimod.c fmpq_mpoly_resultant_interpolation_multimod.c $(pkg-config --cflags --libs flint)

#include <time.h>

#include "mpoly_resultant_interpolation.h"

int COEFF_BIT_SIZE = 100;
int DEFAULT_EXPO_BIT_SIZE = 5;
int DEFAULT_LENGTH = 30;
int DEFAULT_COUNT = 10;
ulong DEFAULT_MODULUS = 65537;
flint_bitcnt_t DEFAULT_MULTIMOD_BITS = 0;
interpolation_mode_t DEFAULT_MULTIMOD_MODE = SMALL_ORDERED_POS;

int fmpz_mpoly_resultant_interpolation_test(
  int expo_bit_size,
  int length,
  int count
) {
  // const char *vars[] = {"x", "y"};
  flint_rand_t rand_state;
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_t P, Q, R1, R2;
  
  clock_t start, end;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpz_mpoly_init(P, ctx);
  fmpz_mpoly_init(Q, ctx);
  fmpz_mpoly_init(R1, ctx);
  fmpz_mpoly_init(R2, ctx);

  int var_to_compute = 0, var_to_evaluate = 1;

  for (int test_case = 0; test_case < count; test_case++) {
    fmpz_mpoly_randtest_bits(P, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);
    fmpz_mpoly_randtest_bits(Q, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);

    start = clock();
    if (fmpz_mpoly_resultant_interpolation_mode(R1, P, Q, var_to_compute, ctx, 0, SMALL_ORDERED)) {
      printf("ERROR: fmpz_mpoly_resultant_interpolation function call returned non-zero error code.\n");
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

    // fmpz_mpoly_print_pretty(P, vars, ctx);
    // printf("\n\n");
    // fmpz_mpoly_print_pretty(Q, vars, ctx);
    // printf("\n\n");
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

int fmpq_mpoly_resultant_interpolation_test(
  int expo_bit_size,
  int length,
  int count
) {
  // const char *vars[] = {"x", "y"};
  flint_rand_t rand_state;
  fmpq_mpoly_ctx_t ctx;
  fmpq_mpoly_t P, Q, R1, R2;
  
  clock_t start, end;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpq_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpq_mpoly_init(P, ctx);
  fmpq_mpoly_init(Q, ctx);
  fmpq_mpoly_init(R1, ctx);
  fmpq_mpoly_init(R2, ctx);

  int var_to_compute = 0, var_to_evaluate = 1;

  for (int test_case = 0; test_case < count; test_case++) {
    fmpq_mpoly_randtest_bits(P, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);
    fmpq_mpoly_randtest_bits(Q, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);

    start = clock();
    if (fmpq_mpoly_resultant_interpolation(R1, P, Q, var_to_compute, ctx)) {
      printf("ERROR: fmpq_mpoly_resultant_interpolation function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Interpolation: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    if (!fmpq_mpoly_resultant(R2, P, Q, var_to_compute, ctx)) {
      printf("ERROR: fmpq_mpoly_resultant function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Naive: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    printf("Real Degree: %ld\n", fmpq_mpoly_degree_si(R2, var_to_evaluate, ctx));

    if (fmpq_mpoly_equal(R1, R2, ctx))
      printf("Test %d is successful.\n", test_case+1);
    else
      printf("Test %d is failed.\n", test_case+1);

    // fmpq_mpoly_print_pretty(P, vars, ctx);
    // printf("\n\n");
    // fmpq_mpoly_print_pretty(Q, vars, ctx);
    // printf("\n\n");
    // fmpq_mpoly_print_pretty(R1, vars, ctx);
    // printf("\n\n");
    // fmpq_mpoly_print_pretty(R2, vars, ctx);
    // printf("\n\n");
  }

cleanup:
  flint_rand_clear(rand_state);
  fmpq_mpoly_clear(P, ctx);
  fmpq_mpoly_clear(Q, ctx);
  fmpq_mpoly_clear(R1, ctx);
  fmpq_mpoly_clear(R2, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  return 0;
}

int fmpz_mpoly_resultant_interpolation_multimod_test(
  int expo_bit_size,
  int length,
  int count,
  flint_bitcnt_t bits,
  interpolation_mode_t mode
) {
  flint_rand_t rand_state;
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_t P, Q, R1, R2, R3;

  clock_t start, end;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpz_mpoly_init(P, ctx);
  fmpz_mpoly_init(Q, ctx);
  fmpz_mpoly_init(R1, ctx);
  fmpz_mpoly_init(R2, ctx);
  fmpz_mpoly_init(R3, ctx);

  int var_to_compute = 0, var_to_evaluate = 1;

  for (int test_case = 0; test_case < count; test_case++) {
    fmpz_mpoly_randtest_bits(P, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);
    fmpz_mpoly_randtest_bits(Q, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);

    start = clock();
    if (fmpz_mpoly_resultant_interpolation_multimod_mode(R1, P, Q, var_to_compute, ctx, bits, mode)) {
      printf("ERROR: fmpz_mpoly_resultant_interpolation_multimod function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Interpolation Multimod: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    if (fmpz_mpoly_resultant_interpolation_mode(R3, P, Q, var_to_compute, ctx, bits, mode)) {
      printf("ERROR: fmpz_mpoly_resultant_interpolation function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Interpolation Pure: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

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
  }

cleanup:
  flint_rand_clear(rand_state);
  fmpz_mpoly_clear(P, ctx);
  fmpz_mpoly_clear(Q, ctx);
  fmpz_mpoly_clear(R1, ctx);
  fmpz_mpoly_clear(R2, ctx);
  fmpz_mpoly_clear(R3, ctx);
  fmpz_mpoly_ctx_clear(ctx);
  return 0;
}

int fmpq_mpoly_resultant_interpolation_multimod_test(
  int expo_bit_size,
  int length,
  int count,
  flint_bitcnt_t bits,
  interpolation_mode_t mode
) {
  flint_rand_t rand_state;
  fmpq_mpoly_ctx_t ctx;
  fmpq_mpoly_t P, Q, R1, R2, R3;

  clock_t start, end;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpq_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpq_mpoly_init(P, ctx);
  fmpq_mpoly_init(Q, ctx);
  fmpq_mpoly_init(R1, ctx);
  fmpq_mpoly_init(R2, ctx);
  fmpq_mpoly_init(R3, ctx);

  int var_to_compute = 0, var_to_evaluate = 1;

  for (int test_case = 0; test_case < count; test_case++) {
    fmpq_mpoly_randtest_bits(P, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);
    fmpq_mpoly_randtest_bits(Q, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);

    start = clock();
    if (fmpq_mpoly_resultant_interpolation_multimod_mode(R1, P, Q, var_to_compute, ctx, bits, mode)) {
      printf("ERROR: fmpq_mpoly_resultant_interpolation_multimod function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Interpolation Multimod: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    if (fmpq_mpoly_resultant_interpolation_mode(R3, P, Q, var_to_compute, ctx, bits, mode)) {
      printf("ERROR: fmpq_mpoly_resultant_interpolation function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Interpolation Pure: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    if (!fmpq_mpoly_resultant(R2, P, Q, var_to_compute, ctx)) {
      printf("ERROR: fmpq_mpoly_resultant function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Naive: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    printf("Real Degree: %ld\n", fmpq_mpoly_degree_si(R2, var_to_evaluate, ctx));

    if (fmpq_mpoly_equal(R1, R2, ctx))
      printf("Test %d is successful.\n", test_case+1);
    else
      printf("Test %d is failed.\n", test_case+1);
  }

cleanup:
  flint_rand_clear(rand_state);
  fmpq_mpoly_clear(P, ctx);
  fmpq_mpoly_clear(Q, ctx);
  fmpq_mpoly_clear(R1, ctx);
  fmpq_mpoly_clear(R2, ctx);
  fmpq_mpoly_clear(R3, ctx);
  fmpq_mpoly_ctx_clear(ctx);
  return 0;
}

int nmod_mpoly_resultant_interpolation_test(
  int expo_bit_size,
  int length,
  int count,
  ulong modulus
) {
  flint_rand_t rand_state;
  nmod_mpoly_ctx_t ctx;
  nmod_mpoly_t P, Q, R1, R2;

  clock_t start, end;

  if (modulus <= 1) {
    printf("ERROR: modulus must be > 1.\n");
    return -1;
  }

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, modulus);
  nmod_mpoly_init(P, ctx);
  nmod_mpoly_init(Q, ctx);
  nmod_mpoly_init(R1, ctx);
  nmod_mpoly_init(R2, ctx);

  int var_to_compute = 0, var_to_evaluate = 1;

  for (int test_case = 0; test_case < count; test_case++) {
    nmod_mpoly_randtest_bits(P, rand_state, length, expo_bit_size, ctx);
    nmod_mpoly_randtest_bits(Q, rand_state, length, expo_bit_size, ctx);

    start = clock();
    if (nmod_mpoly_resultant_interpolation_mode(R1, P, Q, var_to_compute, ctx, 0, SMALL_ORDERED_POS)) {
      printf("ERROR: nmod_mpoly_resultant_interpolation function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Interpolation: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    start = clock();
    if (!nmod_mpoly_resultant(R2, P, Q, var_to_compute, ctx)) {
      printf("ERROR: nmod_mpoly_resultant function call returned non-zero error code.\n");
      goto cleanup;
    }
    end = clock();

    printf("Naive: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    printf("Real Degree: %ld\n", nmod_mpoly_degree_si(R2, var_to_evaluate, ctx));

    if (nmod_mpoly_equal(R1, R2, ctx))
      printf("Test %d is successful.\n", test_case+1);
    else
      printf("Test %d is failed.\n", test_case+1);
  }

cleanup:
  flint_rand_clear(rand_state);
  nmod_mpoly_clear(P, ctx);
  nmod_mpoly_clear(Q, ctx);
  nmod_mpoly_clear(R1, ctx);
  nmod_mpoly_clear(R2, ctx);
  nmod_mpoly_ctx_clear(ctx);
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
  ulong modulus = DEFAULT_MODULUS;
  flint_bitcnt_t multimod_bits = DEFAULT_MULTIMOD_BITS;
  interpolation_mode_t multimod_mode = DEFAULT_MULTIMOD_MODE;

  if (argc > 2) expo_bit_size = atoi(argv[2]);
  if (argc > 3) length = atoi(argv[3]);
  if (argc > 4) count = atoi(argv[4]);
  if (!strcmp("nmod", mode)) {
    if (argc > 5) modulus = strtoul(argv[5], NULL, 10);
  }
  if (!strcmp("fmpz_multimod", mode) || !strcmp("fmpq_multimod", mode)) {
    if (argc > 5) multimod_bits = (flint_bitcnt_t) strtoul(argv[5], NULL, 10);
    if (argc > 6) multimod_mode = (interpolation_mode_t) atoi(argv[6]);
  }

  if (!strcmp("fmpz", mode)) {
    return fmpz_mpoly_resultant_interpolation_test(expo_bit_size, length, count);
  } else if (!strcmp("fmpz_multimod", mode)) {
    return fmpz_mpoly_resultant_interpolation_multimod_test(
      expo_bit_size, length, count, multimod_bits, multimod_mode);
  } else if (!strcmp("fmpq", mode)) {
    return fmpq_mpoly_resultant_interpolation_test(expo_bit_size, length, count);
  } else if (!strcmp("fmpq_multimod", mode)) {
    return fmpq_mpoly_resultant_interpolation_multimod_test(
      expo_bit_size, length, count, multimod_bits, multimod_mode);
  } else if (!strcmp("nmod", mode)) {
    return nmod_mpoly_resultant_interpolation_test(expo_bit_size, length, count, modulus);
  } else {
    printf("ERROR: Unknown mode specified.\n");
    return -1;
  }
}
