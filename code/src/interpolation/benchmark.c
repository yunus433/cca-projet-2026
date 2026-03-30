// gcc-15 -O2 benchmark.c ../utils/utils.h $(pkg-config --cflags --libs flint)

// This file is to test polynomial interpolation in FLINT

// How it works?
// ./a.out {MODE} {OPTIONS}

// All OPTIONS are optional, only MODE is required

// MODE 1: ./a.out DEPENDENCE {DEGREE_START} {DEGREE_END} {COUNT}
// Checks any dependence exist in the format T = a_d * N^k_d for all degrees between [DEGREE_START, DEGREE_END].
// For each (degree, number of points) pair, COUNT polynomials are generated and the mean time is recorded.
// Results are written top data-dependence.csv, can be run with plot.py DEPENDENCE.

// DEGREE_START := Starting degree of polynomials benchmarked.
// DEGREE_END := Ending degree of polynomials benchmarked.
// COUNT := Number of cases to generate.

// Mode 2: ./a.out POW2 {DEGREE_START} {DEGREE_END} {COUNT}
// Checks any dependence exist between the minimum number of points (degree) and the closest pow
// for all degrees between [DEGREE_START, DEGREE_END].
// For each degree, COUNT polynomials are generated for each number of points and the mean time is recorded.
// The ratio of time for power 2 over normal time is recorded for each degree.
// Results are written top data-pow-2.csv, can be run with plot.py POW2.

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

#include "../utils/poly_rand_uniform.c"

int COEFF_BIT_SIZE = 30;
int DEFAULT_DEGREE_START = 20;
int DEFAULT_DEGREE_END = 40;
int DEFAULT_COUNT = 30;

int benchmark_interpolation_flint_dependance(
  int degree_start,
  int degree_end,
  int count
) {
  FILE * data = fopen("./data-dependence.csv", "w");
  clock_t start, end;

  flint_rand_t rand_state;
  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));

  fprintf(data, "degree,n,time\n");

  for (int degree = degree_start; degree < degree_end; degree++) {
    printf("Degree: %d\n", degree);

    fmpz_poly_t P, Q;
    fmpz_poly_init(P);
    fmpz_poly_init(Q);

    int length = degree + 1;

    double totals[degree]; // Degree equals number of points to try, since we try from degree to 2*degree points
    memset(totals, 0, sizeof(totals));

    for (int i = 0; i < count; i++) {
      printf("Count: %d\n", i);

      fmpz_poly_randtest_uniform(P, rand_state, length, COEFF_BIT_SIZE);
      while (fmpz_poly_is_zero(P)) fmpz_poly_randtest_uniform(P, rand_state, length, COEFF_BIT_SIZE);

      for (int n = degree + 1; n < 2 * degree + 1; n++) {
        fmpz x[n], y[n];
        for (int j = 0; j < n; j++) {
          fmpz_init(x + j);
          fmpz_init(y + j);
          fmpz_randbits(x + j, rand_state, COEFF_BIT_SIZE);
          fmpz_poly_evaluate_fmpz(y + j, P, x + j);
        }

        start = clock();
        int return_code = fmpz_poly_interpolate(Q, x, y, n);
        end = clock();

        if (!return_code) {
          printf("ERROR: fmpz_poly_interpolate function call returned non-zero error code.\n");
          break;
        }

        if (!fmpz_poly_equal(P, Q)) {
          printf("ERROR: Interpolation failed.\n", i+1);
          break;
        }

        totals[n-degree-1] += (double)(end - start) / CLOCKS_PER_SEC;
      
        for (int j = 0; j < n; j++) {
          fmpz_clear(x + j);
          fmpz_clear(y + j);
        }
      }
    }

    for (int i = 0; i < degree; i++)
      fprintf(data, "%d,%d,%.18lf\n", degree, degree+i+1, totals[i] / ((double)count));
  
    fmpz_poly_clear(P);
    fmpz_poly_clear(Q);

    printf("\n");
  }

  flint_rand_clear(rand_state);
  fclose(data);

  return 0;
}

// For internal usage, returns the smallest 2^k such that 2^k >= n
int _closest_pow_2(int n) {
  int pow2 = log2(n);
  int check = pow(2, pow2);

  if (check == n)
    return check;
  else
    return 2 * check;
}

int benchmark_interpolation_flint_pow_2(
  int degree_start,
  int degree_end,
  int count
) {
  FILE * data = fopen("./data-pow-2.csv", "w");
  clock_t start, end;

  flint_rand_t rand_state;
  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));

  fprintf(data, "degree,ratio\n");

  for (int degree = degree_start; degree < degree_end; degree++) {
    printf("Degree: %d\n", degree);

    fmpz_poly_t P, Q;
    fmpz_poly_init(P);
    fmpz_poly_init(Q);

    int length = degree + 1;

    int n[2] = {degree + 1, _closest_pow_2(degree + 1)};
    double totals[2] = {0.0, 0.0};

    for (int i = 0; i < count; i++) {
      printf("Count: %d\n", i);

      fmpz_poly_randtest_uniform(P, rand_state, length, COEFF_BIT_SIZE);
      while (fmpz_poly_is_zero(P)) fmpz_poly_randtest_uniform(P, rand_state, length, COEFF_BIT_SIZE);

      for (int j = 0; j < 2; j++) {
        fmpz x[n[j]], y[n[j]];
        for (int k = 0; k < n[j]; k++) {
          fmpz_init(x + k);
          fmpz_init(y + k);
          fmpz_randbits(x + k, rand_state, COEFF_BIT_SIZE);
          fmpz_poly_evaluate_fmpz(y + k, P, x + k);
        }

        start = clock();
        int return_code = fmpz_poly_interpolate(Q, x, y, n[j]);
        end = clock();

        if (!return_code) {
          printf("ERROR: fmpz_poly_interpolate function call returned non-zero error code.\n");
          break;
        }

        if (!fmpz_poly_equal(P, Q)) {
          printf("ERROR: Interpolation failed.\n", i+1);
          break;
        }

        totals[j] += (double)(end - start) / CLOCKS_PER_SEC;
      
        for (int k = 0; k < n[j]; k++) {
          fmpz_clear(x + k);
          fmpz_clear(y + k);
        }
      }
    }

    for (int i = 0; i < 1; i++)
      fprintf(data, "%d,%.18lf\n", degree, totals[1] / totals[0]);
  
    fmpz_poly_clear(P);
    fmpz_poly_clear(Q);

    printf("\n");
  }

  flint_rand_clear(rand_state);
  fclose(data);

  return 0;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("Invalid mode, please see the comment on top of the file to learn how to specify options.\n");
    return -1;
  }

  char * mode = argv[1];

  if (!strcmp("DEPENDENCE", mode)) {
    int degree_start = DEFAULT_DEGREE_START;
    int degree_end = DEFAULT_DEGREE_END;
    int count = DEFAULT_COUNT;

    if (argc > 2) degree_start = atoi(argv[2]);
    if (argc > 3) degree_end = atoi(argv[3]);
    if (argc > 4) count = atoi(argv[4]);

    return benchmark_interpolation_flint_dependance(degree_start, degree_end, count);
  } else if (!strcmp("POW2", mode)) {
    int degree_start = DEFAULT_DEGREE_START;
    int degree_end = DEFAULT_DEGREE_END;
    int count = DEFAULT_COUNT;

    if (argc > 2) degree_start = atoi(argv[2]);
    if (argc > 3) degree_end = atoi(argv[3]);
    if (argc > 4) count = atoi(argv[4]);

    return benchmark_interpolation_flint_pow_2(degree_start, degree_end, count);
  } else {
    printf("Unknown mode specified, please see the comments.");
    return -1;
  }
}