// gcc-15 -O2 fmpz_benchmark.c ../fmpz_mpoly_resultant_interpolation.c $(pkg-config --cflags --libs flint)

#include "../mpoly_resultant_interpolation.h"

int COEFF_BIT_SIZE = 100;
int DEFAULT_EXPO_BIT_SIZE_START = 2;
int DEFAULT_EXPO_BIT_SIZE_END = 5;
int DEFAULT_LENGTH_START = 4;
int DEFAULT_LENGTH_END = 25;
int DEFAULT_LENGTH_INC = 1;
int DEFAULT_COUNT = 20;
int OUTPUT_NAME_LENGTH = 100;

void benchmark(
  int expo_bit_size_start,
  int expo_bit_size_end,
  int length_start,
  int length_end,
  int length_increment,
  int count
) {
  int var_to_compute = 0, var_to_evoluate = 1;

  clock_t start, end;
  flint_rand_t rand_state;
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_t P, Q, R_lib, R_inter;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpz_mpoly_init(P, ctx);
  fmpz_mpoly_init(Q, ctx);
  fmpz_mpoly_init(R_lib, ctx);
  fmpz_mpoly_init(R_inter, ctx);

  for (int expo_bit_size = expo_bit_size_start; expo_bit_size <= expo_bit_size_end; expo_bit_size++) {
    printf("Starting expo bit size %d...\n", expo_bit_size);

    char output_name[OUTPUT_NAME_LENGTH];
    sprintf(output_name, "fmpz_expo_bit_size_%d", expo_bit_size);
    FILE *output = fopen(output_name, "w");

    fprintf(output, "length,lib,rand,rand_pos,rand_small,rand_small_pos,order,order_pos\n");

    for (int length = length_start; length <= length_end; length += length_increment) {
      double t_lib, t_rand, t_rand_pos, t_rand_small, t_rand_small_pos, t_order, t_order_pos;
      t_lib = t_rand = t_rand_pos = t_rand_small = t_rand_small_pos = t_order = t_order_pos = 0.0;

      for (int i = 0; i < count; i++) {
        fmpz_mpoly_randtest_bits(P, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);
        fmpz_mpoly_randtest_bits(Q, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);

        start = clock();
        if (!fmpz_mpoly_resultant(R_lib, P, Q, var_to_compute, ctx)) {
          printf("ERROR: fmpz_mpoly_resultant function call returned non-zero error code.\n");
          goto cleanup;
        }
        end = clock();
        printf("Library call completed.\n");

        // const char *vars[] = {"x", "y"};
        // fmpz_mpoly_print_pretty(R_lib, vars, ctx);
        // printf("\n");
        // fmpz_mpoly_print_pretty(R_inter, vars, ctx);
        // printf("\n");

        slong degree = fmpz_mpoly_degree_si(R_lib, var_to_evoluate, ctx);

        if (degree <= 0) {
          i--;
          continue;
        }

        t_lib += (double)(end - start) / CLOCKS_PER_SEC;

        int bitsize = (degree > 0 ? log2(degree) : 0) + 4;
        printf("Degree: %ld, Bitsize: %d\n", degree, bitsize);

        start = clock();
        if (fmpz_mpoly_resultant_interpolation_mode(R_inter, P, Q, var_to_compute, ctx, COEFF_BIT_SIZE, BITSIZE_RANDOM)) {
          printf("ERROR: fmpz_mpoly_resultant_interpolation function call returned non-zero error code.\n");
          goto cleanup;
        }
        end = clock();
        printf("Bitsize random call completed.\n");

        t_rand += (double)(end - start) / CLOCKS_PER_SEC;

        if (!fmpz_mpoly_equal(R_lib, R_inter, ctx)) {
          printf("ERROR: Two implementations are not equal, mode = BITSIZE_RANDOM\n");
          goto cleanup;
        }

        start = clock();
        if (fmpz_mpoly_resultant_interpolation_mode(R_inter, P, Q, var_to_compute, ctx, COEFF_BIT_SIZE, BITSIZE_RANDOM_POS)) {
          printf("ERROR: fmpz_mpoly_resultant_interpolation function call returned non-zero error code.\n");
          goto cleanup;
        }
        end = clock();
        printf("Bitsize random positive small call completed.\n");

        t_rand_pos += (double)(end - start) / CLOCKS_PER_SEC;

        if (!fmpz_mpoly_equal(R_lib, R_inter, ctx)) {
          printf("ERROR: Two implementations are not equal, mode = BITSIZE_RANDOM_POS\n");
          goto cleanup;
        }

        start = clock();
        if (fmpz_mpoly_resultant_interpolation_mode(R_inter, P, Q, var_to_compute, ctx, bitsize, BITSIZE_RANDOM)) {
          printf("ERROR: fmpz_mpoly_resultant_interpolation function call returned non-zero error code.\n");
          goto cleanup;
        }
        end = clock();
        printf("Bitsize random small call completed.\n");

        t_rand_small += (double)(end - start) / CLOCKS_PER_SEC;

        if (!fmpz_mpoly_equal(R_lib, R_inter, ctx)) {
          printf("ERROR: Two implementations are not equal, mode = BITSIZE_RANDOM_SMALL\n");
          goto cleanup;
        }

        start = clock();
        if (fmpz_mpoly_resultant_interpolation_mode(R_inter, P, Q, var_to_compute, ctx, bitsize, BITSIZE_RANDOM_POS)) {
          printf("ERROR: fmpz_mpoly_resultant_interpolation function call returned non-zero error code.\n");
          goto cleanup;
        }
        end = clock();
        printf("Bitsize random small positive call completed.\n");

        t_rand_small_pos += (double)(end - start) / CLOCKS_PER_SEC;

        if (!fmpz_mpoly_equal(R_lib, R_inter, ctx)) {
          printf("ERROR: Two implementations are not equal, mode = BITSIZE_RANDOM_SMALL_POS\n");
          goto cleanup;
        }

        start = clock();
        if (fmpz_mpoly_resultant_interpolation_mode(R_inter, P, Q, var_to_compute, ctx, 0, SMALL_ORDERED)) {
          printf("ERROR: fmpz_mpoly_resultant_interpolation function call returned non-zero error code.\n");
          goto cleanup;
        }
        end = clock();
        printf("Small ordered call completed.\n");

        t_order += (double)(end - start) / CLOCKS_PER_SEC;

        if (!fmpz_mpoly_equal(R_lib, R_inter, ctx)) {
          printf("ERROR: Two implementations are not equal, mode = SMALL_ORDERED\n");
          goto cleanup;
        }

        start = clock();
        if (fmpz_mpoly_resultant_interpolation_mode(R_inter, P, Q, var_to_compute, ctx, 0, SMALL_ORDERED_POS)) {
          printf("ERROR: fmpz_mpoly_resultant_interpolation function call returned non-zero error code.\n");
          goto cleanup;
        }
        end = clock();
        printf("Small ordered positive call completed.\n");

        t_order_pos += (double)(end - start) / CLOCKS_PER_SEC;

        if (!fmpz_mpoly_equal(R_lib, R_inter, ctx)) {
          printf("ERROR: Two implementations are not equal, mode = SMALL_ORDERED_POS\n");
          goto cleanup;
        }
      }

      t_lib = t_lib / (double)count;
      t_rand = t_rand / (double)count;
      t_rand_pos = t_rand_pos / (double)count;
      t_rand_small = t_rand_small / (double)count;
      t_rand_small_pos = t_rand_small_pos / (double)count;
      t_order = t_order / (double)count;
      t_order_pos = t_order_pos / (double)count;

      fprintf(output, "%d,%.7lf,%.7lf,%.7lf,%.7lf,%.7lf,%.7lf,%.7lf\n", length, t_lib, t_rand, t_rand_pos, t_rand_small, t_rand_small_pos, t_order, t_order_pos);
    }

    fclose(output);
  }

cleanup:
  flint_rand_clear(rand_state);
  fmpz_mpoly_clear(P, ctx);
  fmpz_mpoly_clear(Q, ctx);
  fmpz_mpoly_clear(R_lib, ctx);
  fmpz_mpoly_clear(R_inter, ctx);
  fmpz_mpoly_ctx_clear(ctx);
}

int main(int argc, char *argv[]) {
  int expo_bit_size_start = DEFAULT_EXPO_BIT_SIZE_START;
  int expo_bit_size_end = DEFAULT_EXPO_BIT_SIZE_END;
  int length_start = DEFAULT_LENGTH_START;
  int length_end = DEFAULT_LENGTH_END;
  int length_inc = DEFAULT_LENGTH_INC;
  int count = DEFAULT_COUNT;

  if (argc > 1) expo_bit_size_start = atoi(argv[1]);
  if (argc > 2) expo_bit_size_end = atoi(argv[2]);
  if (argc > 3) length_start = atoi(argv[3]);
  if (argc > 4) length_end = atoi(argv[4]);
  if (argc > 5) length_inc = atoi(argv[5]);
  if (argc > 6) count = atoi(argv[6]);

  benchmark(expo_bit_size_start, expo_bit_size_end, length_start, length_end, length_inc, count);
}