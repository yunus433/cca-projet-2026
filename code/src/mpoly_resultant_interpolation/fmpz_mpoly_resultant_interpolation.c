#include "mpoly_resultant_interpolation.h"

int fmpz_mpoly_resultant_interpolation(
  fmpz_mpoly_t R,
  const fmpz_mpoly_t P,
  const fmpz_mpoly_t Q,
  slong var,
  const fmpz_mpoly_ctx_struct *ctx
) {
  int code = 0, degP, degQ;
  slong var_to_compute = var;
  slong var_to_evoluate = var_to_compute ? 0 : 1;
  slong vars[1], exps[1];

  fmpz_mpoly_t temp_mpoly;
  fmpz_poly_t lead_p, lead_q, temp_poly, temp_poly_2, resultant;

  fmpz_mpoly_init(temp_mpoly, ctx);
  fmpz_poly_init(lead_p);
  fmpz_poly_init(lead_q);
  fmpz_poly_init(temp_poly);
  fmpz_poly_init(temp_poly_2);
  fmpz_poly_init(resultant);

  slong number_of_points = (fmpz_mpoly_degree_si(P, var_to_evoluate, ctx) * fmpz_mpoly_degree_si(Q, var_to_compute, ctx) + fmpz_mpoly_degree_si(Q, var_to_evoluate, ctx) * fmpz_mpoly_degree_si(P, var_to_compute, ctx)) + 1;
  // printf("Number of Points: %ld\n", number_of_points);
  fmpz_t temp;
  fmpz x[number_of_points], y[number_of_points];

  fmpz_init(temp);

  vars[0] = var_to_compute;
  exps[0] = fmpz_mpoly_degree_si(P, var_to_compute, ctx);
  fmpz_mpoly_get_coeff_vars_ui(temp_mpoly, P, vars, exps, 1, ctx);
  fmpz_mpoly_get_fmpz_poly(lead_p, temp_mpoly, var_to_evoluate, ctx);

  vars[0] = var_to_compute;
  exps[0] = fmpz_mpoly_degree_si(Q, var_to_compute, ctx);
  fmpz_mpoly_get_coeff_vars_ui(temp_mpoly, Q, vars, exps, 1, ctx);
  fmpz_mpoly_get_fmpz_poly(lead_q, temp_mpoly, var_to_evoluate, ctx);

  int value = 1;
  for (slong i = 0; i < number_of_points; i++, value++) {
    fmpz_init(x + i);
    fmpz_init(y + i);

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
    if (!fmpz_mpoly_evaluate_one_fmpz(temp_mpoly, P, var_to_evoluate, x + i, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpz_mpoly_get_fmpz_poly(temp_poly, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_get_fmpz_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    if (!fmpz_mpoly_evaluate_one_fmpz(temp_mpoly, Q, var_to_evoluate, x + i, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpz_mpoly_get_fmpz_poly(temp_poly_2, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_get_fmpz_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    fmpz_poly_resultant(y + i, temp_poly, temp_poly_2);
  }

  if (!fmpz_poly_interpolate(resultant, x, y, number_of_points)) {
    printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_poly_interpolate returned non-zero error code.");
    code = -1;
    goto cleanup;
  }

  fmpz_mpoly_set_fmpz_poly(R, resultant, var_to_evoluate, ctx);

cleanup:
  fmpz_mpoly_clear(temp_mpoly, ctx);
  fmpz_poly_clear(lead_p);
  fmpz_poly_clear(lead_q);
  fmpz_poly_clear(temp_poly);
  fmpz_poly_clear(temp_poly_2);
  fmpz_poly_clear(resultant);
  fmpz_clear(temp);

  for (slong i = 0; i < number_of_points; i++) {
    fmpz_clear(x + i);
    fmpz_clear(y + i);
  }

  return code;
}

int fmpz_mpoly_resultant_interpolation_mode(
  fmpz_mpoly_t R,
  const fmpz_mpoly_t P,
  const fmpz_mpoly_t Q,
  slong var,
  const fmpz_mpoly_ctx_struct *ctx,
  const flint_bitcnt_t bits, // This is only used for mode BITSIZE_RANDOM
  interpolation_mode_t mode
) {
  srand(time(NULL));

  int code = 0, degP, degQ;
  slong var_to_compute = var;
  slong var_to_evoluate = var_to_compute ? 0 : 1;
  slong vars[1], exps[1];

  flint_rand_t rand_state;
  fmpz_mpoly_t temp_mpoly;
  fmpz_poly_t lead_p, lead_q, temp_poly, temp_poly_2, resultant;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpz_mpoly_init(temp_mpoly, ctx);
  fmpz_poly_init(lead_p);
  fmpz_poly_init(lead_q);
  fmpz_poly_init(temp_poly);
  fmpz_poly_init(temp_poly_2);
  fmpz_poly_init(resultant);

  slong number_of_points = (fmpz_mpoly_degree_si(P, var_to_evoluate, ctx) * fmpz_mpoly_degree_si(Q, var_to_compute, ctx) + fmpz_mpoly_degree_si(Q, var_to_evoluate, ctx) * fmpz_mpoly_degree_si(P, var_to_compute, ctx)) + 1;
  // printf("Number of Points: %ld\n", number_of_points);
  fmpz_t temp;
  fmpz x[number_of_points], y[number_of_points];

  fmpz_init(temp);

  vars[0] = var_to_compute;
  exps[0] = fmpz_mpoly_degree_si(P, var_to_compute, ctx);
  fmpz_mpoly_get_coeff_vars_ui(temp_mpoly, P, vars, exps, 1, ctx);
  fmpz_mpoly_get_fmpz_poly(lead_p, temp_mpoly, var_to_evoluate, ctx);

  vars[0] = var_to_compute;
  exps[0] = fmpz_mpoly_degree_si(Q, var_to_compute, ctx);
  fmpz_mpoly_get_coeff_vars_ui(temp_mpoly, Q, vars, exps, 1, ctx);
  fmpz_mpoly_get_fmpz_poly(lead_q, temp_mpoly, var_to_evoluate, ctx);

  slong value = 1;
  for (slong i = 0; i < number_of_points; i++) {
    fmpz_init(x + i);
    fmpz_init(y + i);

    if (mode == BITSIZE_RANDOM) {
      fmpz_randbits(x + i, rand_state, bits);

      while (1) {
        int is_all_unique = 1;

        for (int j = 0; j < i; j++)
          if (fmpz_equal(x + j, x + i))
            is_all_unique = 0;

        if (is_all_unique) break;
        fmpz_randbits(x + i, rand_state, bits);
      }
    } else if (mode == BITSIZE_RANDOM_POS) {
      fmpz_randbits_unsigned(x + i, rand_state, bits);

      while (1) {
        int is_all_unique = 1;

        for (int j = 0; j < i; j++)
          if (fmpz_equal(x + j, x + i))
            is_all_unique = 0;

        if (is_all_unique) break;
        fmpz_randbits(x + i, rand_state, bits);
      }
    } else if (mode == SMALL_ORDERED) {
      fmpz_set_si(x + i, value);
      if (value > 0) value *= -1;
      else value = (value * -1) + 1;
    } else if (mode == SMALL_ORDERED_POS) {
      fmpz_set_si(x + i, value);
      value++;
    } else {
      code = -1;
      goto cleanup;
    }

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
    if (!fmpz_mpoly_evaluate_one_fmpz(temp_mpoly, P, var_to_evoluate, x + i, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpz_mpoly_get_fmpz_poly(temp_poly, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_get_fmpz_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    if (!fmpz_mpoly_evaluate_one_fmpz(temp_mpoly, Q, var_to_evoluate, x + i, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpz_mpoly_get_fmpz_poly(temp_poly_2, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_get_fmpz_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    fmpz_poly_resultant(y + i, temp_poly, temp_poly_2);
  }

  if (!fmpz_poly_interpolate(resultant, x, y, number_of_points)) {
    printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_poly_interpolate returned non-zero error code.");
    code = -1;
    goto cleanup;
  }

  fmpz_mpoly_set_fmpz_poly(R, resultant, var_to_evoluate, ctx);

cleanup:
  flint_rand_clear(rand_state);
  fmpz_mpoly_clear(temp_mpoly, ctx);
  fmpz_poly_clear(lead_p);
  fmpz_poly_clear(lead_q);
  fmpz_poly_clear(temp_poly);
  fmpz_poly_clear(temp_poly_2);
  fmpz_poly_clear(resultant);
  fmpz_clear(temp);

  for (slong i = 0; i < number_of_points; i++) {
    fmpz_clear(x + i);
    fmpz_clear(y + i);
  }

  return code;
}