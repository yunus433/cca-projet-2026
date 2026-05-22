#include "mpoly_resultant_interpolation.h"

int nmod_mpoly_resultant_interpolation(
  nmod_mpoly_t R,
  const nmod_mpoly_t P,
  const nmod_mpoly_t Q,
  slong var,
  const nmod_mpoly_ctx_struct *ctx
) {
  if (nmod_mpoly_ctx_nvars(ctx) != 2 || var < 0 || var >= 2)
    return -1;

  int code = 0;
  slong var_to_compute = var;
  slong var_to_evaluate = var_to_compute ? 0 : 1;
  slong vars[1];
  ulong exps[1];

  nmod_mpoly_t temp_mpoly;
  nmod_poly_t lead_p, lead_q, temp_poly, temp_poly_2, resultant;

  if (nmod_mpoly_is_zero(P, ctx) || nmod_mpoly_is_zero(Q, ctx)) {
    nmod_mpoly_zero(R, ctx);
    return 0;
  }

  slong degP_var = nmod_mpoly_degree_si(P, var_to_compute, ctx);
  slong degQ_var = nmod_mpoly_degree_si(Q, var_to_compute, ctx);
  slong degP_eval = nmod_mpoly_degree_si(P, var_to_evaluate, ctx);
  slong degQ_eval = nmod_mpoly_degree_si(Q, var_to_evaluate, ctx);
  if (degP_var < 0 || degQ_var < 0 || degP_eval < 0 || degQ_eval < 0)
    return -1;

  slong number_of_points = degP_eval * degQ_var + degQ_eval * degP_var + 1;

  ulong modulus = nmod_mpoly_ctx_modulus(ctx);
  if (number_of_points > (slong) modulus)
    return -1;

  nmod_mpoly_init(temp_mpoly, ctx);
  nmod_poly_init(lead_p, modulus);
  nmod_poly_init(lead_q, modulus);
  nmod_poly_init(temp_poly, modulus);
  nmod_poly_init(temp_poly_2, modulus);
  nmod_poly_init(resultant, modulus);

  // printf("Number of Points: %ld\n", number_of_points);
  ulong temp = 0;

  ulong *x = NULL;
  ulong *y = NULL;

  x = FLINT_ARRAY_ALLOC(number_of_points, ulong);
  if (x == NULL) {
    code = -1;
    goto cleanup;
  }
  y = FLINT_ARRAY_ALLOC(number_of_points, ulong);
  if (y == NULL) {
    flint_free(x);
    x = NULL;
    code = -1;
    goto cleanup;
  }

  vars[0] = var_to_compute;
  exps[0] = degP_var;
  nmod_mpoly_get_coeff_vars_ui(temp_mpoly, P, vars, exps, 1, ctx);
  if (!nmod_mpoly_get_nmod_poly(lead_p, temp_mpoly, var_to_evaluate, ctx)) {
    printf("ERROR in nmod_mpoly_resultant_interpolation: nmod_mpoly_get_nmod_poly returned non-zero error code.");
    code = -1;
    goto cleanup;
  }

  vars[0] = var_to_compute;
  exps[0] = degQ_var;
  nmod_mpoly_get_coeff_vars_ui(temp_mpoly, Q, vars, exps, 1, ctx);
  if (!nmod_mpoly_get_nmod_poly(lead_q, temp_mpoly, var_to_evaluate, ctx)) {
    printf("ERROR in nmod_mpoly_resultant_interpolation: nmod_mpoly_get_nmod_poly returned non-zero error code.");
    code = -1;
    goto cleanup;
  }

  int value = 1;
  for (slong i = 0; i < number_of_points; i++, value++) {
    x[i] = ((ulong) value) % modulus;

    temp = nmod_poly_evaluate_nmod(lead_p, x[i]);
    if (temp == 0) {
      i--;
      continue;
    }
    temp = nmod_poly_evaluate_nmod(lead_q, x[i]);
    if (temp == 0) {
      i--;
      continue;
    }
  }

  for (slong i = 0; i < number_of_points; i++) {
    nmod_mpoly_evaluate_one_ui(temp_mpoly, P, var_to_evaluate, x[i], ctx);
    if (!nmod_mpoly_get_nmod_poly(temp_poly, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in nmod_mpoly_resultant_interpolation: nmod_mpoly_get_nmod_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    nmod_mpoly_evaluate_one_ui(temp_mpoly, Q, var_to_evaluate, x[i], ctx);
    if (!nmod_mpoly_get_nmod_poly(temp_poly_2, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in nmod_mpoly_resultant_interpolation: nmod_mpoly_get_nmod_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    y[i] = nmod_poly_resultant(temp_poly, temp_poly_2);
  }

  nmod_poly_interpolate_nmod_vec(resultant, x, y, number_of_points);
  nmod_mpoly_set_nmod_poly(R, resultant, var_to_evaluate, ctx);

cleanup:
  nmod_mpoly_clear(temp_mpoly, ctx);
  nmod_poly_clear(lead_p);
  nmod_poly_clear(lead_q);
  nmod_poly_clear(temp_poly);
  nmod_poly_clear(temp_poly_2);
  nmod_poly_clear(resultant);

  if (x != NULL) {
    flint_free(x);
  }
  if (y != NULL) {
    flint_free(y);
  }

  return code;
}

int nmod_mpoly_resultant_interpolation_mode(
  nmod_mpoly_t R,
  const nmod_mpoly_t P,
  const nmod_mpoly_t Q,
  slong var,
  const nmod_mpoly_ctx_struct *ctx,
  const flint_bitcnt_t bits, // This is only used for mode BITSIZE_RANDOM
  interpolation_mode_t mode
) {
  if (nmod_mpoly_ctx_nvars(ctx) != 2 || var < 0 || var >= 2)
    return -1;

  int code = 0;
  slong var_to_compute = var;
  slong var_to_evaluate = var_to_compute ? 0 : 1;
  slong vars[1];
  ulong exps[1];

  flint_rand_t rand_state;
  nmod_mpoly_t temp_mpoly;
  nmod_poly_t lead_p, lead_q, temp_poly, temp_poly_2, resultant;

  if (nmod_mpoly_is_zero(P, ctx) || nmod_mpoly_is_zero(Q, ctx)) {
    nmod_mpoly_zero(R, ctx);
    return 0;
  }

  slong degP_var = nmod_mpoly_degree_si(P, var_to_compute, ctx);
  slong degQ_var = nmod_mpoly_degree_si(Q, var_to_compute, ctx);
  slong degP_eval = nmod_mpoly_degree_si(P, var_to_evaluate, ctx);
  slong degQ_eval = nmod_mpoly_degree_si(Q, var_to_evaluate, ctx);
  if (degP_var < 0 || degQ_var < 0 || degP_eval < 0 || degQ_eval < 0)
    return -1;

  slong number_of_points = degP_eval * degQ_var + degQ_eval * degP_var + 1;

  ulong modulus = nmod_mpoly_ctx_modulus(ctx);
  if (number_of_points > (slong) modulus)
    return -1;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  nmod_mpoly_init(temp_mpoly, ctx);
  nmod_poly_init(lead_p, modulus);
  nmod_poly_init(lead_q, modulus);
  nmod_poly_init(temp_poly, modulus);
  nmod_poly_init(temp_poly_2, modulus);
  nmod_poly_init(resultant, modulus);

  // printf("Number of Points: %ld\n", number_of_points);
  ulong temp = 0;

  ulong *x = NULL;
  ulong *y = NULL;

  x = FLINT_ARRAY_ALLOC(number_of_points, ulong);
  if (x == NULL) {
    code = -1;
    goto cleanup;
  }
  y = FLINT_ARRAY_ALLOC(number_of_points, ulong);
  if (y == NULL) {
    flint_free(x);
    x = NULL;
    code = -1;
    goto cleanup;
  }

  vars[0] = var_to_compute;
  exps[0] = degP_var;
  nmod_mpoly_get_coeff_vars_ui(temp_mpoly, P, vars, exps, 1, ctx);
  if (!nmod_mpoly_get_nmod_poly(lead_p, temp_mpoly, var_to_evaluate, ctx)) {
    printf("ERROR in nmod_mpoly_resultant_interpolation: nmod_mpoly_get_nmod_poly returned non-zero error code.");
    code = -1;
    goto cleanup;
  }

  vars[0] = var_to_compute;
  exps[0] = degQ_var;
  nmod_mpoly_get_coeff_vars_ui(temp_mpoly, Q, vars, exps, 1, ctx);
  if (!nmod_mpoly_get_nmod_poly(lead_q, temp_mpoly, var_to_evaluate, ctx)) {
    printf("ERROR in nmod_mpoly_resultant_interpolation: nmod_mpoly_get_nmod_poly returned non-zero error code.");
    code = -1;
    goto cleanup;
  }

  slong value = 1;
  for (slong i = 0; i < number_of_points; i++) {
    if (mode == BITSIZE_RANDOM) {
      if (bits > 0 && bits < FLINT_BITS)
        x[i] = n_randlimb(rand_state) & ((UWORD(1) << bits) - 1);
      else
        x[i] = n_randlimb(rand_state);
      x[i] %= modulus;

      while (1) {
        int is_all_unique = 1;

        for (int j = 0; j < i; j++)
          if (x[j] == x[i])
            is_all_unique = 0;

        if (is_all_unique) break;

        if (bits > 0 && bits < FLINT_BITS)
          x[i] = n_randlimb(rand_state) & ((UWORD(1) << bits) - 1);
        else
          x[i] = n_randlimb(rand_state);
        x[i] %= modulus;
      }
    } else if (mode == BITSIZE_RANDOM_POS) {
      if (modulus > 1)
        x[i] = 1 + (n_randlimb(rand_state) % (modulus - 1));
      else
        x[i] = 0;

      while (1) {
        int is_all_unique = 1;

        for (int j = 0; j < i; j++)
          if (x[j] == x[i])
            is_all_unique = 0;

        if (is_all_unique) break;

        if (modulus > 1)
          x[i] = 1 + (n_randlimb(rand_state) % (modulus - 1));
        else
          x[i] = 0;
      }
    } else if (mode == SMALL_ORDERED) {
      if (value >= 0) {
        x[i] = ((ulong) value) % modulus;
      } else {
        ulong abs_value = ((ulong) (-value)) % modulus;
        x[i] = (abs_value == 0) ? 0 : modulus - abs_value;
      }

      if (value > 0) value *= -1;
      else value = (value * -1) + 1;
    } else if (mode == SMALL_ORDERED_POS) {
      x[i] = ((ulong) value) % modulus;
      value++;
    } else {
      code = -1;
      goto cleanup;
    }

    temp = nmod_poly_evaluate_nmod(lead_p, x[i]);
    if (temp == 0) {
      i--;
      continue;
    }
    temp = nmod_poly_evaluate_nmod(lead_q, x[i]);
    if (temp == 0) {
      i--;
      continue;
    }
  }

  for (slong i = 0; i < number_of_points; i++) {
    nmod_mpoly_evaluate_one_ui(temp_mpoly, P, var_to_evaluate, x[i], ctx);
    if (!nmod_mpoly_get_nmod_poly(temp_poly, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in nmod_mpoly_resultant_interpolation: nmod_mpoly_get_nmod_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    nmod_mpoly_evaluate_one_ui(temp_mpoly, Q, var_to_evaluate, x[i], ctx);
    if (!nmod_mpoly_get_nmod_poly(temp_poly_2, temp_mpoly, var_to_compute, ctx)) {
      printf("ERROR in nmod_mpoly_resultant_interpolation: nmod_mpoly_get_nmod_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    y[i] = nmod_poly_resultant(temp_poly, temp_poly_2);
  }

  nmod_poly_interpolate_nmod_vec(resultant, x, y, number_of_points);
  nmod_mpoly_set_nmod_poly(R, resultant, var_to_evaluate, ctx);

cleanup:
  flint_rand_clear(rand_state);
  nmod_mpoly_clear(temp_mpoly, ctx);
  nmod_poly_clear(lead_p);
  nmod_poly_clear(lead_q);
  nmod_poly_clear(temp_poly);
  nmod_poly_clear(temp_poly_2);
  nmod_poly_clear(resultant);

  if (x != NULL) {
    flint_free(x);
  }
  if (y != NULL) {
    flint_free(y);
  }

  return code;
}
