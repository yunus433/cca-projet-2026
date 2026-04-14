#include "mpoly_resultant_interpolation.h"

int fmpq_mpoly_resultant_interpolation(
  fmpq_mpoly_t R,
  const fmpq_mpoly_t P,
  const fmpq_mpoly_t Q,
  slong var,
  const fmpq_mpoly_ctx_struct *ctx,
  const flint_bitcnt_t bits
) {
  int code = 0;
  slong var_to_compute = var;
  slong var_to_evoluate = var_to_compute ? 0 : 1;

  flint_rand_t rand_state;
  fmpq_mpoly_t TempM;
  fmpq_poly_t Temp, Temp2, Resultant;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpq_mpoly_init(TempM, ctx);
  fmpq_poly_init(Temp);
  fmpq_poly_init(Temp2);
  fmpq_poly_init(Resultant);

  slong number_of_points = fmpq_mpoly_degree_si(P, var_to_evoluate, ctx) * fmpq_mpoly_degree_si(Q, var_to_compute, ctx) + fmpq_mpoly_degree_si(Q, var_to_evoluate, ctx) * fmpq_mpoly_degree_si(P, var_to_compute, ctx);
  printf("Number of Points: %ld\n", number_of_points);
  fmpq x[number_of_points], y[number_of_points];

  for (slong i = 0; i < number_of_points; i++) {
    fmpq_init(x + i);
    fmpq_init(y + i);

    // fmpq_randbits(x + i, rand_state, bits);
    fmpq_set_si(x + i, i + 1, 1); // 1 2 1/2 3 1/3 3/2 2/3 ...
  }

  for (slong i = 0; i < number_of_points; i++) {
    if (!fmpq_mpoly_evaluate_one_fmpq(TempM, P, var_to_evoluate, x + i, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_evaluate_one_fmpq returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpq_mpoly_get_fmpq_poly(Temp, TempM, var_to_compute, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_get_fmpq_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    if (!fmpq_mpoly_evaluate_one_fmpq(TempM, Q, var_to_evoluate, x + i, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_evaluate_one_fmpq returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpq_mpoly_get_fmpq_poly(Temp2, TempM, var_to_compute, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_get_fmpq_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    fmpq_poly_resultant(y + i, Temp, Temp2);

    printf("%d\n", i);
  }

  fmpq_poly_interpolate_fast(Resultant, x, y, number_of_points);
  fmpq_mpoly_set_fmpq_poly(R, Resultant, var_to_evoluate, ctx);

cleanup:
  flint_rand_clear(rand_state);
  fmpq_mpoly_clear(TempM, ctx);
  fmpq_poly_clear(Temp);
  fmpq_poly_clear(Temp2);
  fmpq_poly_clear(Resultant);

  return code;
}
