// gcc-15 -O2 test_cases_fmpz_mpoly_subresultant_interpolation.c fmpz_mpoly_subresultant_interpolation.c ../subresultant/subresultant_pseudo_remainder.c $(pkg-config --cflags --libs flint)

// As FLINT does not expose any functionality to generate subresultants of bivariate polynomials,
// random tests can only be performed only if the results are pre-known.
// Thus, this code runs all test cases in "test_cases.txt" and checks if the interpolation method gets results correct. 

// Format:
// For each test case:
// First line: The 1st polynomial, use the same format as FLINT if you set by hand (otherwise use GENERATE mode, see below).
// Second line: The 2nd polynomial, use the same format as FLINT if you set by hand (otherwise use GENERATE mode, see below).
// Third line: Subresultants of the above two polynomials, seperate with a comma and space. Note that subresultans will be univariate polynomials.
// The file ends with a new line.

// ATTENTION
// While copying from Mathematica, coefficients do not have a multiplication sign (*) before the variable,
// which is not handled by FLINT. Please add multiplication signs wherever necessary.

// How it works?
// ./a.out {MODE} {OPTIONS}

// All OPTIONS are optional, only MODE is required

// MODE 1: ./a.out GENERATE {DEGREE} {COUNT}
// Generates COUNT new polynomials with the specified DEGREE and writes them at the end fof test_cases.txt
// After running this mode, 2nd line for all new test cases are TO_BE_FILLED by default.
// You should use Wolfram or a similar tool to generate the subresultant of these polynomials and add in the test_cases.txt.

// DEGREE := Degree of polynomials tested.
// COUNT := Number of tests to be performed.

// Mode 2: ./a.out RUN
// Runs all test cases inside test_cases.txt.
// No options provided.

#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "fmpz_mpoly_subresultant_interpolation.h"

#define COEFF_BIT_SIZE 30
#define DEFAULT_COUNT 5
#define DEFAULT_EXPO_BIT_SIZE 2
#define DEFAULT_LENGTH 5
#define MAX_LINE_LENGTH 100000 // This is the maximum char count per line
#define MAX_SUBRESULTANT_COUNT 512 // This is the max nb of subresultants assumed for any polynomial.

int generate_random_polynomials(int length, int expo_bit_size, int count) {
  const char *vars[] = {"x", "y"};

  FILE *test_cases = fopen("./test_cases.txt", "w");

  flint_rand_t rand_state;
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_t A, B;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);

  fmpz_mpoly_init(A, ctx);
  fmpz_mpoly_init(B, ctx);

  printf("Starting generation...\n");

  for (int i = 0; i < count; i++) {
    fmpz_mpoly_randtest_bits(A, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);
    while (fmpz_mpoly_is_zero(A, ctx))
      fmpz_mpoly_randtest_bits(A, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);

    fmpz_mpoly_randtest_bits(B, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);
    while (fmpz_mpoly_is_zero(B, ctx))
      fmpz_mpoly_randtest_bits(B, rand_state, length, COEFF_BIT_SIZE, expo_bit_size, ctx);

    fmpz_mpoly_fprint_pretty(test_cases, A, vars, ctx);
    fprintf(test_cases, "\n");
    fmpz_mpoly_fprint_pretty(test_cases, B, vars, ctx);
    fprintf(test_cases, "\n");
    fprintf(test_cases, "TO BE FILLED: subresultant1,subresultant2,subresultant3,...\n");
  }

  printf("Generation completed\n");

  flint_rand_clear(rand_state);
  fmpz_mpoly_clear(A, ctx);
  fmpz_mpoly_clear(B, ctx);
  fmpz_mpoly_ctx_clear(ctx);

  fclose(test_cases);

  return 0;
}

void remove_endline(char * line) {
  for (int i = 0; line[i] != '\0'; i++)
    if (line[i] == '\n')
      line[i] = ' ';
}

int run_test_cases() {
  const char *vars[] = {"x", "y"};
  slong var_to_compute = 0, var_to_evaluate = 1;

  FILE *test_cases = fopen("./test_cases.txt", "r");

  char line[MAX_LINE_LENGTH];
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_t P, Q;
  fmpz_mpoly_t subresultants[MAX_SUBRESULTANT_COUNT];
  fmpz_mpoly_t checks[MAX_SUBRESULTANT_COUNT];

  fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
  fmpz_mpoly_init(P, ctx);
  fmpz_mpoly_init(Q, ctx);

  for (int i = 0; i < MAX_SUBRESULTANT_COUNT; i++) {
    fmpz_mpoly_init(subresultants[i], ctx);
    fmpz_mpoly_init(checks[i], ctx);
  }

  int test_count = 0;
  int successful_test_count = 0;
  int successful_test_count_pseudo_remainder = 0;

  while (fgets(line, sizeof(line), test_cases)) {
    remove_endline(line);
    if (fmpz_mpoly_set_str_pretty(P, line, vars, ctx)) {
      printf("ERROR: Invalid input file format.\n");
      goto safe_exit;
    }

    if (!fgets(line, sizeof(line), test_cases)) {
      printf("ERROR: Invalid input file format.\n");
      goto safe_exit;
    }
    
    remove_endline(line);
    if (fmpz_mpoly_set_str_pretty(Q, line, vars, ctx)) {
      printf("ERROR: Invalid input file format.\n");
      goto safe_exit;
    }

    if (!fgets(line, sizeof(line), test_cases)) {
      printf("ERROR: Invalid input file format.\n");
      goto safe_exit;
    }

    char *subresultant = strtok(line, ",\n");
    int subresultant_count = 0;
    while (subresultant != NULL) {
      if (fmpz_mpoly_set_str_pretty(checks[subresultant_count], subresultant, vars, ctx)) {
        printf("ERROR: Invalid input file format.\n");
        goto safe_exit;
      }

      subresultant = strtok(NULL, ",\n");
      subresultant_count++;
    }

    test_count++;

    for (int i = 0; i < MAX_SUBRESULTANT_COUNT; i++)
      fmpz_mpoly_zero(subresultants[i], ctx);

    if (fmpz_mpoly_subresultant_interpolation(subresultants, P, Q, var_to_compute, ctx)) {
      printf("ERROR: fmpz_mpoly_subresultant_interpolation function call returned non-zero error code on test case %d.\n", test_count);
    } else {
      int is_all_successful = 1;

      for (int i = 0; i < subresultant_count && is_all_successful; i++)
        if (!fmpz_mpoly_equal(subresultants[i], checks[i], ctx)) {
          fmpz_mpoly_print_pretty(subresultants[i], vars, ctx);
          printf("\n");
          fmpz_mpoly_print_pretty(checks[i], vars, ctx);
          printf("\n");
          is_all_successful = 0;
        }
          

      if (!is_all_successful) {
        printf("FAIL: Test case %d\n", test_count);
      } else {
        printf("SUCCESS: Test case %d\n", test_count);
        successful_test_count++;
      }
    }
  }

  printf("\n------------\n");
  printf("Results: %d / %d tests passed.\n", successful_test_count, test_count);

safe_exit:
  fmpz_mpoly_clear(P, ctx);
  fmpz_mpoly_clear(Q, ctx);
  
  for (int i = 0; i < MAX_SUBRESULTANT_COUNT; i++) {
    fmpz_mpoly_clear(subresultants[i], ctx);
    fmpz_mpoly_clear(checks[i], ctx);
  }

  fmpz_mpoly_ctx_clear(ctx);
  fclose(test_cases);

  return 0;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    printf("Invalid mode, please see the comment on top of the file to learn how to specify options.\n");
    return -1;
  }

  char * mode = argv[1];

  if (!strcmp("GENERATE", mode)) {
    int length = DEFAULT_LENGTH;
    int expo_bit_size = DEFAULT_EXPO_BIT_SIZE;
    int count = DEFAULT_COUNT;

    if (argc > 2) length = atoi(argv[2]);
    if (argc > 3) expo_bit_size = atoi(argv[3]);
    if (argc > 4) count = atoi(argv[4]);

    return generate_random_polynomials(length, expo_bit_size, count);
  } else if  (!strcmp("RUN", mode)) {
    return run_test_cases();
  } else {
    printf("Unknown mode specified, please see the comments.");
    return -1;
  }
}