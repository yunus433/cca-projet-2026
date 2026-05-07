#include "fmpz_mpoly_subresultant_interpolation.h"

int fmpz_mpoly_subresultant_interpolation(
  fmpz_mpoly_t * subresultants,
  const fmpz_mpoly_t P,
  const fmpz_mpoly_t Q,
  slong var,
  const fmpz_mpoly_ctx_struct *ctx
) {
  if (fmpz_mpoly_ctx_nvars(ctx) != 2 || var < 0 || var >= 2)
    return -1;

  int code = 0;
  slong var_to_compute = var;
  slong var_to_evaluate = var_to_compute ? 0 : 1;
  slong vars[1];
  ulong exps[1];

  fmpz_mpoly_t temp_mpoly;
  fmpz_poly_t lead_p, lead_q, temp_poly, temp_poly_2, resultant;

  slong degP_var = fmpz_mpoly_degree_si(P, var_to_compute, ctx);
  slong degQ_var = fmpz_mpoly_degree_si(Q, var_to_compute, ctx);
  slong degP_eval = fmpz_mpoly_degree_si(P, var_to_evaluate, ctx);
  slong degQ_eval = fmpz_mpoly_degree_si(Q, var_to_evaluate, ctx);
  if (degP_var < 0 || degQ_var < 0 || degP_eval < 0 || degQ_eval < 0)
    return -1;

  int K = (degP_var > degQ_var ? degP_var : degQ_var) + 1; // number of subresultants

  if (fmpz_mpoly_is_zero(P, ctx) || fmpz_mpoly_is_zero(Q, ctx)) {
    for (int k = 0; k < K; k++)
      fmpz_mpoly_zero(subresultants[k], ctx);
    return 0;
  }

  slong number_of_points = degP_eval * degQ_var + degQ_eval * degP_var + 1;

  fmpz_mpoly_init(temp_mpoly, ctx);
  fmpz_poly_init(lead_p);
  fmpz_poly_init(lead_q);
  fmpz_poly_init(temp_poly);
  fmpz_poly_init(temp_poly_2);
  fmpz_poly_init(resultant);

  // printf("Number of Points: %ld\n", number_of_points);
  fmpz_t temp;
  fmpz_init(temp);

  fmpz *x = FLINT_ARRAY_ALLOC(number_of_points, fmpz);
  if (x == NULL) {
    printf("ERROR in fmpz_mpoly_subresultant_interpolation: FLINT_ARRAY_ALLOC returned NULL pointer, memory allocation error.\n");
    code = -1;
    goto cleanup;
  }
  fmpz *y = FLINT_ARRAY_ALLOC(K * number_of_points, fmpz);
  if (y == NULL) {
    printf("ERROR in fmpz_mpoly_subresultant_interpolation: FLINT_ARRAY_ALLOC returned NULL pointer, memory allocation error.\n");
    code = -1;
    goto cleanup;
  }
  fmpz *_y = FLINT_ARRAY_ALLOC(number_of_points, fmpz); // _y is the flattened cyclic version of y, used for interpolation variable assignment.
  if (_y == NULL) {
    printf("ERROR in fmpz_mpoly_subresultant_interpolation: FLINT_ARRAY_ALLOC returned NULL pointer, memory allocation error.\n");
    code = -1;
    goto cleanup;
  }

  for (slong i = 0; i < number_of_points; i++) {
    fmpz_init(x + i);
    fmpz_init(_y + i);
  }
  for (slong i = 0; i < (K + 1) * number_of_points; i++) {
    fmpz_init(y + i);
  }

  vars[0] = var_to_compute;
  exps[0] = degP_var;
  fmpz_mpoly_get_coeff_vars_ui(temp_mpoly, P, vars, exps, 1, ctx);
  if (!fmpz_mpoly_get_fmpz_poly(lead_p, temp_mpoly, var_to_evaluate, ctx)) {
    printf("ERROR in fmpz_mpoly_subresultant_interpolation: fmpz_mpoly_get_fmpz_poly returned non-zero error code.\n");
    code = -1;
    goto cleanup;
  }

  vars[0] = var_to_compute;
  exps[0] = degQ_var;
  fmpz_mpoly_get_coeff_vars_ui(temp_mpoly, Q, vars, exps, 1, ctx);
  if (!fmpz_mpoly_get_fmpz_poly(lead_q, temp_mpoly, var_to_evaluate, ctx)) {
    printf("ERROR in fmpz_mpoly_subresultant_interpolation: fmpz_mpoly_get_fmpz_poly returned non-zero error code.\n");
    code = -1;
    goto cleanup;
  }

  int value = 1;
  for (slong i = 0; i < number_of_points; i++, value++) {
    fmpz_set_si(x + i, value);

    fmpz_poly_evaluate_fmpz(temp, lead_p, x + i);
    if (fmpz_is_zero(temp)) {
      i--;
      continue;
    }
    fmpz_poly_evaluate_fmpz(temp, lead_q, x + i);
    if (fmpz_is_zero(temp)) {
      i--;
      continue;
    }
  }

  for (slong i = 0; i < number_of_points; i++) {
    if (!fmpz_mpoly_evaluate_one_fmpz(temp_mpoly, P, var_to_evaluate, x + i, ctx)) {
      printf("ERROR in fmpz_mpoly_subresultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.\n");
      code = -1;
      goto cleanup;
    }
    if (!fmpz_mpoly_get_fmpz_poly(temp_poly, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in fmpz_mpoly_subresultant_interpolation: fmpz_mpoly_get_fmpz_poly returned non-zero error code.\n");
      code = -1;
      goto cleanup;
    }

    if (!fmpz_mpoly_evaluate_one_fmpz(temp_mpoly, Q, var_to_evaluate, x + i, ctx)) {
      printf("ERROR in fmpz_mpoly_subresultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.\n");
      code = -1;
      goto cleanup;
    }
    if (!fmpz_mpoly_get_fmpz_poly(temp_poly_2, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in fmpz_mpoly_subresultant_interpolation: fmpz_mpoly_get_fmpz_poly returned non-zero error code.\n");
      code = -1;
      goto cleanup;
    }

    fmpz_poly_subresultant_pseudo_remainder(y + (K * i), temp_poly, temp_poly_2);
  }

  for (int k = 0; k < K; k++) {
    for (int i = k, j = 0; i < K * number_of_points; i += K, j++) {
      fmpz_set(_y + j, y + i);
    }

    if (!fmpz_poly_interpolate(resultant, x, _y, number_of_points)) {
      printf("ERROR in fmpz_mpoly_subresultant_interpolation: fmpz_poly_interpolate returned non-zero error code.\n");
      code = -1;
      goto cleanup;
    }

    fmpz_mpoly_set_fmpz_poly(subresultants[k], resultant, var_to_evaluate, ctx);
  }

cleanup:
  fmpz_mpoly_clear(temp_mpoly, ctx);
  fmpz_poly_clear(lead_p);
  fmpz_poly_clear(lead_q);
  fmpz_poly_clear(temp_poly);
  fmpz_poly_clear(temp_poly_2);
  fmpz_poly_clear(resultant);
  fmpz_clear(temp);

  if (x != NULL) {
    for (slong i = 0; i < number_of_points; i++) {
      fmpz_clear(x + i);
    }
    flint_free(x);
  }
  if (y != NULL) {
    for (slong i = 0; i < K * number_of_points; i++) {
      fmpz_clear(y + i);
    }
    flint_free(y);
  }
  if (_y != NULL) {
    for (slong i = 0; i < number_of_points; i++) {
      fmpz_clear(_y + i);
    }
    flint_free(_y);
  }

  return code;
}
