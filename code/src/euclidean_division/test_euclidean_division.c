#include "euclidean_division.h"

int _fmpq_poly_euclidean_division_test(
  int degree,
  int count
) {
  flint_rand_t state;
  fmpq_poly_t A, B;

  flint_rand_init(state);

  for (int i = 0; i < count; i++) {
    printf("Starting test %d...\n", i+1);
    
    fmpq_poly_init(B);
    fmpq_poly_randtest(B, state, /*len*/ 10, /*bits*/ 30);   // Q[x]
    fmpq_poly_print(B); flint_printf("\n");
  }

  flint_rand_clear(state);
}

int main() {

}