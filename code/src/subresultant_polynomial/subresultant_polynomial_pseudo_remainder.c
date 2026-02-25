// gcc-15 -O2 subresultant_polynomial_pseudo_remainder.c $(pkg-config --cflags --libs flint) -oa

#include "subresultant_polynomial_pseudo_remainder.h"

int fmpz_poly_subresultant_polynomial_pseudo_remainder(
  fmpz_poly_t * subresultant_polynomials,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
) {
  fmpz_poly_t r0, r1, tmp_poly;
  fmpz_poly_init(r0);
  fmpz_poly_init(r1);
  fmpz_poly_init(tmp_poly);

  if (fmpz_poly_degree(Q) > fmpz_poly_degree(P)) {
    fmpz_poly_set(r0, Q);
    fmpz_poly_set(r1, P);
  } else {
    fmpz_poly_set(r0, P);
    fmpz_poly_set(r1, Q);
  }

  int i = 1;
  fmpz_t gama, beta, psi, d, tmp;
  fmpz_init(gama);
  fmpz_init(beta);
  fmpz_init(psi);
  fmpz_init(d);
  fmpz_init(tmp);

  while (fmpz_poly_degree(r1) > 0) {
    if (i == 1) {
      slong dValue = fmpz_poly_degree(r0) - fmpz_poly_degree(r1);
      fmpz_set_si(d, dValue);
      fmpz_set(gama, fmpz_poly_lead(r1));
      fmpz_set_si(beta, (dValue + 1) % 2 ? -1 : 1);
      fmpz_set_si(psi, -1);
    } else {
      fmpz_poly_set(subresultant_polynomials[i - 2], r1);

      fmpz_add_si(tmp, d, -1);
      fmpz_pow_fmpz(psi, psi, tmp);
      fmpz_set(tmp, gama);
      fmpz_mul_si(tmp, tmp, -1);
      fmpz_pow_fmpz(tmp, tmp, d);
      fmpz_divexact(psi, tmp, psi);

      fmpz_set_si(d, fmpz_poly_degree(r0) - fmpz_poly_degree(r1));

      fmpz_pow_fmpz(beta, psi, d);
      fmpz_mul(beta, beta, gama);
      fmpz_mul_si(beta, beta, -1);

      fmpz_set(gama, fmpz_poly_lead(r1));
    }

    fmpz_pow_fmpz(tmp, gama, d);
    fmpz_mul(tmp, tmp, gama);

    fmpz_poly_scalar_mul_fmpz(r0,  r0, tmp);

    fmpz_poly_rem(tmp_poly, r0, r1);

    fmpz_poly_scalar_divexact_fmpz(tmp_poly, tmp_poly, beta);
    fmpz_poly_set(r0, r1);
    fmpz_poly_set(r1, tmp_poly);

    i++;
  }

  fmpz_poly_set(subresultant_polynomials[i - 2], r1);

  fmpz_poly_clear(r0);
  fmpz_poly_clear(r1);
  fmpz_poly_clear(tmp_poly);

  fmpz_clear(gama);
  fmpz_clear(beta);
  fmpz_clear(psi);
  fmpz_clear(d);
  fmpz_clear(tmp);

  return 0;
}

// int main () {
//   fmpz_poly_t f, g;
//   fmpz_poly_init(f);
//   fmpz_poly_init(g);

//   // Example: f = 824x^5 - 65x^4 -814x^3 - 741x^2 - 979x- 764, g = 216x^4 + 663x^3 + 880x^2 + 916x + 617
//   // fmpz_poly_set_coeff_si(g, 0, -764);
//   // fmpz_poly_set_coeff_si(g, 1, -979);
//   // fmpz_poly_set_coeff_si(g, 2, -741);
//   // fmpz_poly_set_coeff_si(g, 3, -814);
//   // fmpz_poly_set_coeff_si(g, 4, -65);
//   // fmpz_poly_set_coeff_si(g, 5, 824);
//   fmpz_poly_set_coeff_si(g, 0, -5);
//   fmpz_poly_set_coeff_si(g, 1, 2);
//   fmpz_poly_set_coeff_si(g, 2, 8);
//   fmpz_poly_set_coeff_si(g, 3, -3);
//   fmpz_poly_set_coeff_si(g, 4, -3);
//   fmpz_poly_set_coeff_si(g, 6, 1);
//   fmpz_poly_set_coeff_si(g, 8, 1);

//   // fmpz_poly_set_coeff_si(f, 0, 617);
//   // fmpz_poly_set_coeff_si(f, 1, 916);
//   // fmpz_poly_set_coeff_si(f, 2, 880);
//   // fmpz_poly_set_coeff_si(f, 3, 663);
//   // fmpz_poly_set_coeff_si(f, 4, 216);
//   fmpz_poly_set_coeff_si(f, 0, 21);
//   fmpz_poly_set_coeff_si(f, 1, -9);
//   fmpz_poly_set_coeff_si(f, 2, -4);
//   fmpz_poly_set_coeff_si(f, 4, 5);
//   fmpz_poly_set_coeff_si(f, 6, 3);
  
//   fmpz_poly_t subresultant_polys[7];

//   for (int i = 0; i < 7; i++)
//     fmpz_poly_init(subresultant_polys[i]);

//   fmpz_subresultant_polynomial_pseudo_remainder(subresultant_polys,  g, f);
//   for (int i = 0; i < 7; i++) {
//     fmpz_poly_print_pretty(subresultant_polys[i], "x");
//     printf("\n");
//     fmpz_poly_clear(subresultant_polys[i]);
//   }

//   fmpz_poly_clear(f);
//   fmpz_poly_clear(g);

//   return 0;
// }