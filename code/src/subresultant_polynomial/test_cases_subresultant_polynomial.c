// gcc-15 -O2 ../sylvester_matrix/polynomial_sylvester_matrix.c subresultant_polynomial_naive.c subresultant_polynomial_pseudo_remainder.c test_cases_subresultant_polynomial.c  $(pkg-config --cflags --libs flint)

// As FLINT does not expose any functionality to generate subresultant polynomials, random tests can only be
// performed only if the results are pre-known.
// Thus, this code runs all test cases in "test_cases.txt" and checks if all subresultant polynomial algorithms
// get results correct. 

// Format:
// For each test case:
// First line: The 1st polynomial, use the same format as FLINT if you set by hand (otherwise use GENERATE mode, see below).
// Second line: The 2nd polynomial, use the same format as FLINT if you set by hand (otherwise use GENERATE mode, see below).
// Third line: Subresultant polynomials of the above two polynomials, seperate with a comma and space.
// The file ends with a new line.

// How it works?
// ./a.out {MODE} {OPTIONS}

// All OPTIONS are optional, only MODE is required

// MODE 1: ./a.out GENERATE {DEGREE} {COUNT}
// Generates COUNT new polynomials with the specified DEGREE and writes them at the end of test_cases.txt
// After running this mode, 2nd line for all new test cases are TO_BE_FILLED by default.
// You should use Wolfram or a similar tool to generate the subresultant polynomials of these polynomials and add in the test_cases.txt.

// DEGREE := Degree of polynomials tested.
// COUNT := Number of tests to be performed.

// Mode 2: ./a.out RUN
// Runs all test cases inside test_cases.txt.
// No options provided.

#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "../utils/poly_set_str_pretty.c"

#include "subresultant_polynomial_naive.h"
#include "subresultant_polynomial_pseudo_remainder.h"

#define COEFF_BIT_SIZE 30
#define DEFAULT_COUNT 10
#define DEFAULT_DEGREE 20
#define MAX_LINE_LENGTH 200000 // This is the maximum char count per line
#define MAX_SUBRESULTANT_COUNT 512 // This is the max nb of subresultant polynomials assumed for any polynomial.
#define NUMBER_BASE 10 // All numbers are in base 10

void fmpz_poly_normalize(fmpz_poly_t poly) {
    if (fmpz_poly_is_zero(poly))
        return;

    fmpz_t content;
    fmpz_init(content);

    // Compute GCD of coefficients
    fmpz_poly_content(content, poly);

    // Divide by content if not 1
    if (!fmpz_is_one(content))
        fmpz_poly_scalar_divexact_fmpz(poly, poly, content);

    // Make leading coefficient positive
    if (fmpz_sgn(fmpz_poly_lead(poly)) < 0)
        fmpz_poly_neg(poly, poly);

    fmpz_clear(content);
}

int generate_random_polynomials(int degree, int count) {
  FILE *test_cases = fopen("./test_cases.txt", "w");

  flint_rand_t rand_state;
  fmpz_poly_t A, B;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));

  fmpz_poly_init(A);
  fmpz_poly_init(B);

  int length = degree + 1;

  printf("Starting generation...\n");

  for (int i = 0; i < count; i++) {
    fmpz_poly_randtest(A, rand_state, length, COEFF_BIT_SIZE);
    while (fmpz_poly_is_zero(A))
      fmpz_poly_randtest(A, rand_state, length, COEFF_BIT_SIZE);

    fmpz_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);
    while (fmpz_poly_is_zero(B))
      fmpz_poly_randtest(B, rand_state, length, COEFF_BIT_SIZE);

    fmpz_poly_fprint_pretty(test_cases, A, "x");
    fprintf(test_cases, "\n");
    fmpz_poly_fprint_pretty(test_cases, B, "x");
    fprintf(test_cases, "\n");
    fprintf(test_cases, "TO BE FILLED: subresultant_polynomial_1,subresultant_polynomial_2,subresultant_polynomial_3,...\n");
  }

  printf("Generation completed\n");

  flint_rand_clear(rand_state);
  fmpz_poly_clear(A);
  fmpz_poly_clear(B);

  fclose(test_cases);

  return 0;
}

int run_test_cases() {
  FILE *test_cases = fopen("./test_cases.txt", "r");

  char line[MAX_LINE_LENGTH];
  fmpz_poly_t P, Q, TEMP;
  fmpz_poly_t subresultant_polynomials[MAX_SUBRESULTANT_COUNT];
  fmpz_poly_t checks[MAX_SUBRESULTANT_COUNT];

  fmpz_poly_init(P);
  fmpz_poly_init(Q);

    fmpz_poly_init(TEMP);

  for (int i = 0; i < MAX_SUBRESULTANT_COUNT; i++) {
    fmpz_poly_init(subresultant_polynomials[i]);
    fmpz_poly_init(checks[i]);
  }

  int test_count = 0;
  int successful_test_count_naive = 0;
  int successful_test_count_pseudo_remainder = 0;

  while (fgets(line, sizeof(line), test_cases)) {
    if (fmpz_poly_set_str_pretty(P, line, "x")) {
      printf("ERROR: Invalid input file format.\n");
      goto safe_exit;
    }

    if (!fgets(line, sizeof(line), test_cases)) {
      printf("ERROR: Invalid input file format.\n");
      goto safe_exit;
    }
    
    if (fmpz_poly_set_str_pretty(Q, line, "x")) {
      printf("ERROR: Invalid input file format.\n");
      goto safe_exit;
    }

    if (!fgets(line, sizeof(line), test_cases)) {
      printf("ERROR: Invalid input file format.\n");
      goto safe_exit;
    }


    char *subresultant_polynomial = strtok(line, ",\n");
    int subresultant_polynomial_count = 0;
    while (subresultant_polynomial != NULL) {
      if (fmpz_poly_set_str_pretty(checks[subresultant_polynomial_count], subresultant_polynomial, "x")) {
        printf("ERROR: Invalid input file format.\n");
        goto safe_exit;
      }
      subresultant_polynomial = strtok(NULL, ",\n");
      subresultant_polynomial_count++;
    }

    test_count++;

    for (int i = 0; i < MAX_SUBRESULTANT_COUNT; i++)
      fmpz_poly_zero(subresultant_polynomials[i]);

    if (fmpz_poly_subresultant_polynomial_naive(subresultant_polynomials, P, Q)) {
      printf("ERROR: subresultant_naive function call returned non-zero error code on test case %d.\n", test_count);
    } else {
      int is_all_successful = 1;

      for (int i = 0; i < subresultant_polynomial_count; i++) {        
        fmpz_poly_normalize(subresultant_polynomials[i]);

        if (!fmpz_poly_equal(subresultant_polynomials[i], checks[i])) {
          printf("HERE %d\n", i);
          fmpz_poly_print_pretty(subresultant_polynomials[i], "x");
          printf("\n\n");
          fmpz_poly_print_pretty(checks[i], "x");
          fmpz_poly_div(TEMP, checks[i], subresultant_polynomials[i]);
          printf("\n\n");
          fmpz_poly_print_pretty(TEMP, "x");
          fmpz_poly_rem(TEMP, checks[i], subresultant_polynomials[i]);
          printf("\n\n");
          fmpz_poly_print_pretty(TEMP, "x");
          printf("\n\n\n");
          is_all_successful = 0;
        }
      }
          

      if (!is_all_successful) {
        printf("FAIL: Naive - Test case %d\n", test_count);
      } else {
        printf("SUCCESS: Naive - Test case %d\n", test_count);
        successful_test_count_naive++;
      }
    }

    for (int i = 0; i < MAX_SUBRESULTANT_COUNT; i++)
      fmpz_poly_zero(subresultant_polynomials[i]);

    // if (fmpz_poly_subresultant_polynomial_pseudo_remainder(subresultant_polynomials, P, Q)) {
    //   printf("ERROR: subresultant_pseudo_remainder function call returned non-zero error code on test case %d.\n", test_count);
    // } else {
    //   int is_all_successful = 1;

    //   for (int i = 0; i < subresultant_polynomial_count && is_all_successful; i++)
    //     if (!fmpz_poly_equal(subresultant_polynomials[i], checks[i])) {
    //       fmpz_poly_print_pretty(subresultant_polynomials[i], "x");
    //       printf("\n");
    //       fmpz_poly_print_pretty(checks[i], "x");
    //       printf("\n");
    //       is_all_successful = 0;
    //     }
          

    //   if (!is_all_successful) {
    //     printf("FAIL: Pseudo Remainder - Test case %d\n", test_count);
    //   } else {
    //     printf("SUCCESS: Pseudo Remainder - Test case %d\n", test_count);
    //     successful_test_count_pseudo_remainder++;
    //   }
    // }
  }

  printf("\n------------\n");
  printf("Naive: %d / %d tests passed.\n", successful_test_count_naive, test_count);
  printf("Pseudo Remainder: %d / %d tests passed.\n", successful_test_count_pseudo_remainder, test_count);

safe_exit:

  fmpz_poly_clear(P);
  fmpz_poly_clear(Q);
    fmpz_poly_clear(TEMP);
  
  for (int i = 0; i < MAX_SUBRESULTANT_COUNT; i++) {
    fmpz_poly_clear(subresultant_polynomials[i]);
    fmpz_poly_clear(checks[i]);
  }

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
    int degree = DEFAULT_DEGREE;
    int count = DEFAULT_COUNT;

    if (argc > 2) degree = atoi(argv[1]);
    if (argc > 3) count = atoi(argv[2]);

    return generate_random_polynomials(degree, count);
  } else if  (!strcmp("RUN", mode)) {
    return run_test_cases();
  } else {
    printf("Unknown mode specified, please see the comments.");
    return -1;
  }
}