#include "mpoly_resultant_interpolation.h"

int fmpq_mpoly_resultant_interpolation(
  fmpq_mpoly_t R,
  const fmpq_mpoly_t P,
  const fmpq_mpoly_t Q,
  slong var,
  const fmpq_mpoly_ctx_struct *ctx
) {
  int code = 0;
  slong var_to_compute = var;
  slong var_to_evaluate = var_to_compute ? 0 : 1;
  slong vars[1], exps[1];

  fmpq_mpoly_t temp_mpoly;
  fmpq_poly_t lead_p, lead_q, temp_poly, temp_poly_2, resultant;

  fmpq_mpoly_init(temp_mpoly, ctx);
  fmpq_poly_init(lead_p);
  fmpq_poly_init(lead_q);
  fmpq_poly_init(temp_poly);
  fmpq_poly_init(temp_poly_2);
  fmpq_poly_init(resultant);

  slong number_of_points = (fmpq_mpoly_degree_si(P, var_to_evaluate, ctx) * fmpq_mpoly_degree_si(Q, var_to_compute, ctx) + fmpq_mpoly_degree_si(Q, var_to_evaluate, ctx) * fmpq_mpoly_degree_si(P, var_to_compute, ctx)) + 1;
  // printf("Number of Points: %ld\n", number_of_points);
  fmpq_t value, temp;
  fmpq x[number_of_points], y[number_of_points];

  fmpq_init(value);
  fmpq_init(temp);

  vars[0] = var_to_compute;
  exps[0] = fmpq_mpoly_degree_si(P, var_to_compute, ctx);
  fmpq_mpoly_get_coeff_vars_ui(temp_mpoly, P, vars, exps, 1, ctx);
  fmpq_mpoly_get_fmpq_poly(lead_p, temp_mpoly, var_to_evaluate, ctx);

  vars[0] = var_to_compute;
  exps[0] = fmpq_mpoly_degree_si(Q, var_to_compute, ctx);
  fmpq_mpoly_get_coeff_vars_ui(temp_mpoly, Q, vars, exps, 1, ctx);
  fmpq_mpoly_get_fmpq_poly(lead_q, temp_mpoly, var_to_evaluate, ctx);

  fmpq_set_si(value, 1, 1);
  for (slong i = 0; i < number_of_points; i++) {
    fmpq_init(x + i);
    fmpq_init(y + i);

    fmpq_set(x + i, value);
    fmpq_next_signed_calkin_wilf(value, value);

    fmpq_poly_evaluate_fmpq(temp, lead_p, x + i);
    if (fmpq_is_zero(temp)) {
      i--;
      continue;
    }
    fmpq_poly_evaluate_fmpq(temp, lead_q, x + i);
    if (fmpq_is_zero(temp)) {
      i--;
      continue;
    }
  }

  for (slong i = 0; i < number_of_points; i++) {
    if (!fmpq_mpoly_evaluate_one_fmpq(temp_mpoly, P, var_to_evaluate, x + i, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_evaluate_one_fmpq returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpq_mpoly_get_fmpq_poly(temp_poly, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_get_fmpq_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    if (!fmpq_mpoly_evaluate_one_fmpq(temp_mpoly, Q, var_to_evaluate, x + i, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_evaluate_one_fmpq returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpq_mpoly_get_fmpq_poly(temp_poly_2, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_get_fmpq_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    fmpq_poly_resultant(y + i, temp_poly, temp_poly_2);
  }

  fmpq_poly_interpolate_fast(resultant, x, y, number_of_points);
  fmpq_mpoly_set_fmpq_poly(R, resultant, var_to_evaluate, ctx);

cleanup:
  fmpq_mpoly_clear(temp_mpoly, ctx);
  fmpq_poly_clear(lead_p);
  fmpq_poly_clear(lead_q);
  fmpq_poly_clear(temp_poly);
  fmpq_poly_clear(temp_poly_2);
  fmpq_poly_clear(resultant);
  fmpq_clear(value);
  fmpq_clear(temp);

  for (slong i = 0; i < number_of_points; i++) {
    fmpq_clear(x + i);
    fmpq_clear(y + i);
  }

  return code;
}

int fmpq_mpoly_resultant_interpolation_mode(
  fmpq_mpoly_t R,
  const fmpq_mpoly_t P,
  const fmpq_mpoly_t Q,
  slong var,
  const fmpq_mpoly_ctx_struct *ctx,
  const flint_bitcnt_t bits, // This is only used for mode BITSIZE_RANDOM
  interpolation_mode_t mode
) {
  srand(time(NULL));

  int code = 0;
  slong var_to_compute = var;
  slong var_to_evaluate = var_to_compute ? 0 : 1;
  slong vars[1], exps[1];

  flint_rand_t rand_state;
  fmpq_mpoly_t temp_mpoly;
  fmpq_poly_t lead_p, lead_q, temp_poly, temp_poly_2, resultant;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpq_mpoly_init(temp_mpoly, ctx);
  fmpq_poly_init(lead_p);
  fmpq_poly_init(lead_q);
  fmpq_poly_init(temp_poly);
  fmpq_poly_init(temp_poly_2);
  fmpq_poly_init(resultant);

  slong number_of_points = (fmpq_mpoly_degree_si(P, var_to_evaluate, ctx) * fmpq_mpoly_degree_si(Q, var_to_compute, ctx) + fmpq_mpoly_degree_si(Q, var_to_evaluate, ctx) * fmpq_mpoly_degree_si(P, var_to_compute, ctx)) + 1;
  // printf("Number of Points: %ld\n", number_of_points);
  fmpq_t value, temp;
  fmpq x[number_of_points], y[number_of_points];

  fmpq_init(value);
  fmpq_init(temp);

  vars[0] = var_to_compute;
  exps[0] = fmpq_mpoly_degree_si(P, var_to_compute, ctx);
  fmpq_mpoly_get_coeff_vars_ui(temp_mpoly, P, vars, exps, 1, ctx);
  fmpq_mpoly_get_fmpq_poly(lead_p, temp_mpoly, var_to_evaluate, ctx);

  vars[0] = var_to_compute;
  exps[0] = fmpq_mpoly_degree_si(Q, var_to_compute, ctx);
  fmpq_mpoly_get_coeff_vars_ui(temp_mpoly, Q, vars, exps, 1, ctx);
  fmpq_mpoly_get_fmpq_poly(lead_q, temp_mpoly, var_to_evaluate, ctx);

  fmpq_set_si(value, 1, 1);
  for (slong i = 0; i < number_of_points; i++) {
    fmpq_init(x + i);
    fmpq_init(y + i);

    if (mode == BITSIZE_RANDOM) {
      fmpq_randbits(x + i, rand_state, bits);

      while (1) {
        int is_all_unique = 1;

        for (int j = 0; j < i; j++)
          if (fmpq_equal(x + j, x + i))
            is_all_unique = 0;

        if (is_all_unique) break;
        fmpq_randbits(x + i, rand_state, bits);
      }
    } else if (mode == BITSIZE_RANDOM_POS) { // Not defined properly
      fmpq_randbits(x + i, rand_state, bits);

      while (1) {
        int is_all_unique = 1;

        for (int j = 0; j < i; j++)
          if (fmpq_equal(x + j, x + i))
            is_all_unique = 0;

        if (is_all_unique) break;
        fmpq_randbits(x + i, rand_state, bits);
      }
    } else if (mode == SMALL_ORDERED) {
      fmpq_set(x + i, value);
      fmpq_next_signed_calkin_wilf(value, value);
    } else if (mode == SMALL_ORDERED_POS) {
      fmpq_set(x + i, value);
      fmpq_next_calkin_wilf(value, value);
    } else {
      code = -1;
      goto cleanup;
    }

    fmpq_poly_evaluate_fmpq(temp, lead_p, x + i);
    if (fmpq_is_zero(temp)) {
      i--;
      continue;
    }
    fmpq_poly_evaluate_fmpq(temp, lead_q, x + i);
    if (fmpq_is_zero(temp)) {
      i--;
      continue;
    }
  }

  for (slong i = 0; i < number_of_points; i++) {
    if (!fmpq_mpoly_evaluate_one_fmpq(temp_mpoly, P, var_to_evaluate, x + i, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_evaluate_one_fmpq returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpq_mpoly_get_fmpq_poly(temp_poly, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_get_fmpq_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    if (!fmpq_mpoly_evaluate_one_fmpq(temp_mpoly, Q, var_to_evaluate, x + i, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_evaluate_one_fmpq returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpq_mpoly_get_fmpq_poly(temp_poly_2, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in fmpq_mpoly_resultant_interpolation: fmpq_mpoly_get_fmpq_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    fmpq_poly_resultant(y + i, temp_poly, temp_poly_2);
  }

  fmpq_poly_interpolate_fast(resultant, x, y, number_of_points);
  fmpq_mpoly_set_fmpq_poly(R, resultant, var_to_evaluate, ctx);

cleanup:
  flint_rand_clear(rand_state);
  fmpq_mpoly_clear(temp_mpoly, ctx);
  fmpq_poly_clear(lead_p);
  fmpq_poly_clear(lead_q);
  fmpq_poly_clear(temp_poly);
  fmpq_poly_clear(temp_poly_2);
  fmpq_poly_clear(resultant);
  fmpq_clear(value);
  fmpq_clear(temp);

  for (slong i = 0; i < number_of_points; i++) {
    fmpq_clear(x + i);
    fmpq_clear(y + i);
  }

  return code;
}