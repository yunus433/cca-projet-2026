// gcc-15 -O2 test_poly_set_Str_pretty.c utils.h $(pkg-config --cflags --libs flint)

#include <time.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

#include "poly_set_str_pretty.c"

#define TEST_CASE_COUNT 4

char * TEST_CASES[TEST_CASE_COUNT] = {
  "18460856*x^6-195399*x^5+364791*x^4+312704*x^3",
  "-157055*x^20-327541*x^16+402*x^15-1407*x^13+4079*x^9-20*x^8+65536*x^5-4607*x^4+x^2",
  "21*x^20-3043650*x^19+x^18-524289*x^17-182683*x^16+84737*x^15-483024315*x^14+250574*x^13-37*x^12+39036798*x^9+8*x^8+232190*x^7-2178*x^6-17*x^5-2097151*x^4-1",
  "256*x^15"
};

int main() {
  fmpz_poly_t P;

  fmpz_poly_init(P);

  for (int i = 0; i < TEST_CASE_COUNT; i++) {
    if (fmpz_poly_set_str_pretty(P, TEST_CASES[i], "x")) {
      printf("ERROR: Test case %d returned non-zero error code.\n", i+1);
      break;
    }

    printf("Test Case %d:\n%s\nPolynomial:\n", i+1, TEST_CASES[i]);
    fmpz_poly_print_pretty(P, "x");
    printf("\n");
  }

  fmpz_poly_clear(P);
}