// DOES NOT WORK - not complete

#include "trivar_fmpz_mpoly_resultant_interpolation.h"

/*
    xs length = dx + 1
    ys length = dy + 1

    vals[b][a] = P(xs[a], ys[b])
    so rows are fixed y, varying x.
*/
int interpolate_bivariate_fmpz(
    fmpz_mpoly_t P,
    const fmpz *xs,
    const fmpz *ys,
    fmpz **vals,
    slong dx,
    slong dy,
    const fmpz_mpoly_ctx_t ctx
) {
    slong i, j;
    int ok = 1;

    fmpz_poly_t *row_poly;
    fmpz_poly_t col_poly;
    fmpz_t coeff;
    ulong exp[2];

    row_poly = flint_malloc((dy + 1) * sizeof(fmpz_poly_t));

    for (j = 0; j <= dy; j++)
    {
        fmpz_poly_init(row_poly[j]);

        /* Interpolate P(x, ys[j]) in x */
        if (!fmpz_poly_interpolate(row_poly[j], xs, vals[j], dx + 1))
        {
            ok = 0;
            goto cleanup_rows;
        }
    }

    fmpz_poly_init(col_poly);
    fmpz_init(coeff);

    fmpz_mpoly_zero(P, ctx);

    for (i = 0; i <= dx; i++)
    {
        fmpz *col_vals = _fmpz_vec_init(dy + 1);

        for (j = 0; j <= dy; j++)
        {
            /* coefficient of x^i in row_poly[j] */
            fmpz_poly_get_coeff_fmpz(col_vals + j, row_poly[j], i);
        }

        /* Interpolate coefficient of x^i as a polynomial in y */
        if (!fmpz_poly_interpolate(col_poly, ys, col_vals, dy + 1))
        {
            _fmpz_vec_clear(col_vals, dy + 1);
            ok = 0;
            goto cleanup;
        }

        for (j = 0; j <= dy; j++)
        {
            fmpz_poly_get_coeff_fmpz(coeff, col_poly, j);

            if (!fmpz_is_zero(coeff))
            {
                exp[0] = i;   /* x exponent */
                exp[1] = j;   /* y exponent */

                fmpz_mpoly_set_coeff_fmpz_ui(P, coeff, exp, ctx);
            }
        }

        _fmpz_vec_clear(col_vals, dy + 1);
    }

cleanup:
    fmpz_clear(coeff);
    fmpz_poly_clear(col_poly);

cleanup_rows:
    for (j = 0; j <= dy; j++)
        fmpz_poly_clear(row_poly[j]);

    flint_free(row_poly);

    return ok;
}

int trivar_fmpz_mpoly_resultant_interpolation(
  fmpz_mpoly_t R,
  const fmpz_mpoly_t P,
  const fmpz_mpoly_t Q,
  slong var,
  const fmpz_mpoly_ctx_struct *ctx
) {
  if (fmpz_mpoly_ctx_nvars(ctx) != 3 || var < 0 || var >= 3)
    return -1;

  if (fmpz_mpoly_is_zero(P, ctx) || fmpz_mpoly_is_zero(Q, ctx))
  {
    fmpz_mpoly_zero(R, ctx);
    return 0;
  }

  int code = 0;
  slong var_to_compute = var;
  slong var_to_evaluate[2] = {0, 0};
  if (var_to_compute == 0)
  {
    var_to_evaluate[0] = 1;
    var_to_evaluate[1] = 2;
  }
  else if (var_to_compute == 1)
  {
    var_to_evaluate[1] = 2;
  }
  else if (var_to_compute == 2) {
    var_to_evaluate[1] = 1;
  }
  slong vars[2];
  ulong exps[2];

  fmpz_mpoly_t lead_p, lead_q, temp_0, temp_1, temp_2;

  slong degP_var = fmpz_mpoly_degree_si(P, var_to_compute, ctx);
  slong degQ_var = fmpz_mpoly_degree_si(Q, var_to_compute, ctx);
  slong degP_eval_0 = fmpz_mpoly_degree_si(P, var_to_evaluate[0], ctx);
  slong degQ_eval_0 = fmpz_mpoly_degree_si(Q, var_to_evaluate[0], ctx);
  slong degP_eval_1 = fmpz_mpoly_degree_si(P, var_to_evaluate[1], ctx);
  slong degQ_eval_1 = fmpz_mpoly_degree_si(Q, var_to_evaluate[1], ctx);
  if (degP_var < 0 || degQ_var < 0 || degP_eval_0 < 0 || degQ_eval_0 < 0 || degP_eval_1 < 0 || degQ_eval_1 < 0)
    return -1;

  slong number_of_points_0 = degP_eval_0 * degQ_var + degQ_eval_0 * degP_var + 1;
  slong number_of_points_1 = degP_eval_1 * degQ_var + degQ_eval_1 * degP_var + 1;
  slong number_of_points = number_of_points_0 > number_of_points_1 ? number_of_points_0 : number_of_points_1;
  number_of_points = number_of_points * number_of_points; // We take 2nd power since this is a bivariate interpolation

  fmpz_mpoly_init(lead_p, ctx);
  fmpz_mpoly_init(lead_q, ctx);
  fmpz_mpoly_init(resultant, ctx);
  fmpz_mpoly_init(temp_0, ctx);
  fmpz_mpoly_init(temp_1, ctx);
  fmpz_mpoly_init(temp_2, ctx);

  fmpz_t temp;
  fmpz_init(temp);

  fmpz *x = FLINT_ARRAY_ALLOC(number_of_points, fmpz);
  if (x == NULL) {
    code = -1;
    goto cleanup;
  }
  fmpz *y = FLINT_ARRAY_ALLOC(number_of_points, fmpz);
  if (y == NULL) {
    flint_free(x);
    x = NULL;
    code = -1;
    goto cleanup;
  }
  fmpz_poly_t *z = FLINT_ARRAY_ALLOC(number_of_points, fmpz_poly_t);
  if (z == NULL) {
    flint_free(x);
    x = NULL;
    flint_free(y);
    y = NULL;
    code = -1;
    goto cleanup;
  }
  for (slong i = 0; i < number_of_points; i++) {
    fmpz_init(x + i);
    fmpz_init(y + i);
    fmpz_poly_init(z[i]);
  }

  vars[0] = var_to_compute;
  exps[0] = degP_var;
  fmpz_mpoly_get_coeff_vars_ui(lead_p, P, vars, exps, 1, ctx);
  exps[0] = degQ_var;
  fmpz_mpoly_get_coeff_vars_ui(lead_q, Q, vars, exps, 1, ctx);

  // const char *vars_[] = {"x", "y", "z"};
  // fmpz_mpoly_print_pretty(lead_p, vars_, ctx);
  // printf("\n");
  // fmpz_mpoly_print_pretty(lead_q, vars_, ctx);
  // printf("\n");

  slong value_bound = (slong) n_sqrt((ulong) number_of_points) + 1;
  slong value_0 = 1, value_1 = 1;
  for (slong i = 0; i < number_of_points; i++) {
    fmpz_set_si(x + i, value_0);
    fmpz_set_si(y + i, value_1);
    fmpz *vals[2] = {x + i, y + i};

    if (value_1 < value_bound)
      value_1++;
    else {
      value_0++;
      value_1 = 1;
    }

    fmpz_mpoly_evaluate_all_fmpz(temp, lead_p, vals, ctx);
    if (fmpz_is_zero(temp)) {
      i--;
      continue;
    }
    fmpz_mpoly_evaluate_all_fmpz(temp, lead_q, vals, ctx);
    if (fmpz_is_zero(temp)) {
      i--;
      continue;
    }
  }

  for (slong i = 0; i < number_of_points; i++) {
    if (!fmpz_mpoly_evaluate_one_fmpz(temp_0, P, var_to_evaluate[0], x + i, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpz_mpoly_evaluate_one_fmpz(temp_0, temp_0, var_to_evaluate[1], y + i, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    if (!fmpz_mpoly_evaluate_one_fmpz(temp_1, Q, var_to_evaluate[0], x + i, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpz_mpoly_evaluate_one_fmpz(temp_1, temp_1, var_to_evaluate[1], y + i, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    fmpz_mpoly_resultant(z[i], temp_0, temp_1, var_to_compute, ctx);
  }

  if (!interpolate_bivariate_fmpz(R, x, y, number_of_points)) {
    printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_poly_interpolate returned non-zero error code.");
    code = -1;
    goto cleanup;
  }

cleanup:
  fmpz_mpoly_clear(lead_p, ctx);
  fmpz_mpoly_clear(lead_q, ctx);
  fmpz_mpoly_clear(temp_0, ctx);
  fmpz_mpoly_clear(temp_1, ctx);
  fmpz_mpoly_clear(temp_2, ctx);
  fmpz_clear(temp);

  if (x != NULL) {
    for (slong i = 0; i < number_of_points; i++) {
      fmpz_clear(x + i);
    }
    flint_free(x);
  }
  if (y != NULL) {
    for (slong i = 0; i < number_of_points; i++) {
      fmpz_clear(y + i);
    }
    flint_free(y);
  }
  if (z != NULL) {
    for (slong i = 0; i < number_of_points; i++) {
      fmpz_poly_clear(z[i]);
    }
    flint_free(z);
  }

  return code;
}