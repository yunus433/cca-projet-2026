#include "mpoly_resultant_interpolation.h"

int fmpz_mpoly_resultant_interpolation(
  fmpz_mpoly_t R,
  const fmpz_mpoly_t P,
  const fmpz_mpoly_t Q,
  slong var,
  const fmpz_mpoly_ctx_struct *ctx,
  const flint_bitcnt_t bits
) {
  int code = 0;
  slong var_to_compute = var;
  slong var_to_evoluate = var_to_compute ? 0 : 1;

  flint_rand_t rand_state;
  fmpz_mpoly_t TempM;
  fmpz_poly_t Temp, Temp2, Resultant;

  flint_rand_init(rand_state);
  flint_rand_set_seed(rand_state, time(NULL), time(NULL));
  fmpz_mpoly_init(TempM, ctx);
  fmpz_poly_init(Temp);
  fmpz_poly_init(Temp2);
  fmpz_poly_init(Resultant);

  slong number_of_points = fmpz_mpoly_degree_si(P, var_to_evoluate, ctx) * (fmpz_mpoly_length(Q, ctx) + 1) + fmpz_mpoly_degree_si(Q, var_to_evoluate, ctx) * (fmpz_mpoly_length(P, ctx) + 1);
  printf("Number of Points: %ld\n", number_of_points);
  fmpz x[number_of_points], y[number_of_points];

  for (slong i = 0; i < number_of_points; i++) {
    fmpz_init(x + i);
    fmpz_init(y + i);

    fmpz_randbits(x + i, rand_state, bits);
    fmpz_randbits(y + i, rand_state, bits);
  }

  for (slong i = 0; i < number_of_points; i++) {
    if (!fmpz_mpoly_evaluate_one_fmpz(TempM, P, var_to_evoluate, x + i, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpz_mpoly_get_fmpz_poly(Temp, TempM, var_to_compute, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_get_fmpz_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    if (!fmpz_mpoly_evaluate_one_fmpz(TempM, Q, var_to_evoluate, x + i, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
      code = -1;
      goto cleanup;
    }
    if (!fmpz_mpoly_get_fmpz_poly(Temp2, TempM, var_to_compute, ctx)) {
      printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_mpoly_get_fmpz_poly returned non-zero error code.");
      code = -1;
      goto cleanup;
    }

    fmpz_poly_resultant(y + i, Temp, Temp2);
  }

  if (!fmpz_poly_interpolate(Resultant, x, y, number_of_points)) {
    printf("ERROR in fmpz_mpoly_resultant_interpolation: fmpz_poly_interpolate returned non-zero error code.");
    code = -1;
    goto cleanup;
  }

  fmpz_mpoly_set_fmpz_poly(R, Resultant, var_to_evoluate, ctx);

cleanup:
  flint_rand_clear(rand_state);
  fmpz_mpoly_clear(TempM, ctx);
  fmpz_poly_clear(Temp);
  fmpz_poly_clear(Temp2);
  fmpz_poly_clear(Resultant);

  return code;
}

// void resultant_from_interpolation() {
//   const char *vars[] = {"x", "y"};
//   int var_to_evoluate = 0, var_to_compute = 1;

//   flint_rand_t rand_state;
//   fmpz_mpoly_ctx_t ctx;
//   fmpz_mpoly_t P, Q, R, TempM;
//   fmpz_poly_t Temp, Temp2, Resultant, Resultant2;

//   flint_rand_init(rand_state);
//   flint_rand_set_seed(rand_state, time(NULL), time(NULL));
//   fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
//   fmpz_mpoly_init(P, ctx);
//   fmpz_mpoly_init(Q, ctx);
//   fmpz_mpoly_init(R, ctx);
//   fmpz_mpoly_init(TempM, ctx);
//   fmpz_poly_init(Temp);
//   fmpz_poly_init(Temp2);
//   fmpz_poly_init(Resultant);
//   fmpz_poly_init(Resultant2);

//   fmpz_mpoly_randtest_bits(P, rand_state, DEFAULT_LENGTH, COEFF_BIT_SIZE, EXPO_BIT_SIZE, ctx);
//   fmpz_mpoly_randtest_bits(Q, rand_state, DEFAULT_LENGTH, COEFF_BIT_SIZE, EXPO_BIT_SIZE, ctx);

//   slong number_of_points = fmpz_mpoly_degree_si(P, var_to_compute, ctx) * fmpz_mpoly_length(Q, ctx) + fmpz_mpoly_degree_si(Q, var_to_compute, ctx) * fmpz_mpoly_length(P, ctx);
//   printf("Number of Points: %ld\n", number_of_points);
//   fmpz x[number_of_points], y[number_of_points];

//   for (slong i = 0; i < number_of_points; i++) {
//     fmpz_init(x + i);
//     fmpz_init(y + i);

//     fmpz_randbits(x + i, rand_state, COEFF_BIT_SIZE);
//     fmpz_randbits(y + i, rand_state, COEFF_BIT_SIZE);
//   }

//   for (slong i = 0; i < number_of_points; i++) {
//     if (!fmpz_mpoly_evaluate_one_fmpz(TempM, P, var_to_evoluate, x + i, ctx)) {
//       printf("ERROR: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
//       goto cleanup;
//     }
//     if (!fmpz_mpoly_get_fmpz_poly(Temp, TempM, var_to_compute, ctx)) {
//       printf("ERROR: fmpz_mpoly_get_fmpz_poly returned non-zero error code.");
//       goto cleanup;
//     }

//     if (!fmpz_mpoly_evaluate_one_fmpz(TempM, Q, var_to_evoluate, x + i, ctx)) {
//       printf("ERROR: fmpz_mpoly_evaluate_one_fmpz returned non-zero error code.");
//       goto cleanup;
//     }
//     if (!fmpz_mpoly_get_fmpz_poly(Temp2, TempM, var_to_compute, ctx)) {
//       printf("ERROR: fmpz_mpoly_get_fmpz_poly returned non-zero error code.");
//       goto cleanup;
//     }

//     fmpz_poly_resultant(y + i, Temp, Temp2);
//   }

//   fmpz_poly_interpolate(Resultant, x, y, number_of_points);

//   printf("Real Degree: %ld\n", fmpz_poly_degree(Resultant));

//   if (!fmpz_mpoly_resultant(R, P, Q, var_to_compute, ctx)) {
//     printf("ERROR: fmpz_mpoly_resultant function call returned non-zero error code.\n");
//     goto cleanup;
//   }

//   if (!fmpz_mpoly_get_fmpz_poly(Resultant2, R, var_to_evoluate, ctx)) {
//     printf("ERROR: fmpz_mpoly_get_fmpz_poly returned non-zero error code.");
//     goto cleanup;
//   }

//   if (fmpz_poly_equal(Resultant, Resultant2))
//     printf("SUCCESS!!!!\n");
//   else
//     printf("Error :((\n");

//   // printf("Interpolation Method:\n");
//   // fmpz_poly_print_pretty(Resultant, "x");
//   // printf("\nNormal Method:\n");
//   // fmpz_poly_print_pretty(Resultant2, "x");
//   // printf("\n");

//   for (slong i = 0; i < number_of_points; i++) {
//     fmpz_clear(x + i);
//     fmpz_clear(y + i);
//   }

// cleanup:
//   flint_rand_clear(rand_state);
//   fmpz_mpoly_clear(P, ctx);
//   fmpz_mpoly_clear(Q, ctx);
//   fmpz_mpoly_clear(R, ctx);
//   fmpz_mpoly_clear(TempM, ctx);
//   fmpz_poly_clear(Temp);
//   fmpz_poly_clear(Temp2);
//   fmpz_poly_clear(Resultant);
//   fmpz_poly_clear(Resultant2);
//   fmpz_mpoly_ctx_clear(ctx);
// }

// int main() {
//   resultant_from_interpolation();
// }