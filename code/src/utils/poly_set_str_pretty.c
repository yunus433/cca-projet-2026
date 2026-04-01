// FLINT does not have fmpz_poly_set_str_pretty() function, which is definitely a bug since fmpz_mpoly_set_str_pretty exist.
// Thus, we make a work around that will read in the pretty format from a temporary file.
// This function has the "supposed" signature for fmpz_poly_set_str_pretty().
// It is mainly used in test_cases since we need to be able to read pretty polynomials from files, as they are outputed by Wolfram.

#include "utils.h"

char* DEFAULT_VARIABLE = "x";

int fmpz_poly_set_str_pretty(
  fmpz_poly_struct *poly,
  const char *str,
  const char *x
) {
  FILE *fp;
  const char *expected_var = x == NULL ? DEFAULT_VARIABLE : x;
  char *returned_var = NULL;
  int ret;

  fp = fmemopen((void *) str, strlen(str), "r");
  if (fp == NULL) return -1;

  ret = fmpz_poly_fread_pretty(fp, poly, &returned_var);
  fclose(fp);

  if (ret <= 0) {
    fmpz_poly_zero(poly);
    return -1;
  }

  if (*returned_var && *returned_var != *expected_var) {
    fmpz_poly_zero(poly);
    return -1;
  }

  return 0;
}