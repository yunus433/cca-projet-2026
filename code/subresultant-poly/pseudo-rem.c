// gcc-15 -O2 pseudo-rem.c $(pkg-config --cflags --libs flint) -oa
// QUESTION: Comparison with wikipedia?

#ifndef SUBRESPOLYPSEUDO_H
#define SUBRESPOLYPSEUDO_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_mat.h>

slong pow_slong(slong a, slong b) {
  slong ans = 1;

  if (b == 1) return a;
  if (b == 2) return a*a;

  while (b > 0) {
    if (b % 2) {
      ans *= a;
      b--;
    }

    if (b > 1)
      ans *= ans;
    b >>= 2;
  }

  return ans;
}

int fmpz_subresultant_polynomials_pseudo_remainder_sequence(
  fmpz_poly_t * subresultant_polynomials,
  const fmpz_poly_t P,
  const fmpz_poly_t Q
) {
  fmpz_poly_t r0, r1, tmp;
  fmpz_poly_init(r0);
  fmpz_poly_init(r1);
  fmpz_poly_init(tmp);

  if (fmpz_poly_degree(Q) > fmpz_poly_degree(P)) {
    fmpz_poly_set(r0, Q);
    fmpz_poly_set(r1, P);
  } else {
    fmpz_poly_set(r0, P);
    fmpz_poly_set(r1, Q);
  }

  int i = 1;
  slong gama, beta, psi, d;
  // Changer les types...

  while (fmpz_poly_degree(r1) > 0) {
    printf("%d\n", i);
    fmpz_poly_print_pretty(r0, "x");
    printf("\n");
    fmpz_poly_print_pretty(r1, "x");
    printf("\n");

    if (i == 1) {
      d = fmpz_poly_degree(r0) - fmpz_poly_degree(r1);
      gama = * fmpz_poly_lead(r1);
      beta = (d+1) % 2 ? -1 : 1;
      psi = -1;
    } else {
      fmpz_poly_set(subresultant_polynomials[i - 2], r1);

      psi = pow_slong(-1 * gama, d) / pow_slong(psi, d - 1);
      d = fmpz_poly_degree(r0) - fmpz_poly_degree(r1);
      beta = -1 * gama * pow_slong(psi, d);
      gama = * fmpz_poly_lead(r1);
    }

    printf("d: %lld, beta: %lld, gama: %lld\n", d, beta, gama);
    printf("pow %lld\n", pow_slong(gama, d+1));
    fmpz_poly_scalar_mul_si(r0,  r0, pow_slong(gama, d+1));
    printf("HERE\n");
    fmpz_poly_print_pretty(r0, "x");
    printf("\n");
    fmpz_poly_print_pretty(r1, "x");
    printf("\n");
    fmpz_poly_rem(tmp, r0, r1);
    fmpz_poly_print_pretty(tmp, "x");
    printf("\n");
    printf("HERE END\n");
    fmpz_poly_scalar_tdiv_si(tmp, tmp, beta);
    fmpz_poly_set(r0, r1);
    fmpz_poly_set(r1, tmp);

    i++;
  }

  fmpz_poly_set(subresultant_polynomials[i - 2], r1);

  return 1;
}

#endif

int main () {
  fmpz_poly_t f, g;
  fmpz_poly_init(f);
  fmpz_poly_init(g);

  // Example: f = 824x^5 - 65x^4 -814x^3 - 741x^2 - 979x- 764, g = 216x^4 + 663x^3 + 880x^2 + 916x + 617
  // fmpz_poly_set_coeff_si(g, 0, -764);
  // fmpz_poly_set_coeff_si(g, 1, -979);
  // fmpz_poly_set_coeff_si(g, 2, -741);
  // fmpz_poly_set_coeff_si(g, 3, -814);
  // fmpz_poly_set_coeff_si(g, 4, -65);
  // fmpz_poly_set_coeff_si(g, 5, 824);
  fmpz_poly_set_coeff_si(g, 0, -5);
  fmpz_poly_set_coeff_si(g, 1, 2);
  fmpz_poly_set_coeff_si(g, 2, 8);
  fmpz_poly_set_coeff_si(g, 3, -3);
  fmpz_poly_set_coeff_si(g, 4, -3);
  fmpz_poly_set_coeff_si(g, 6, 1);
  fmpz_poly_set_coeff_si(g, 8, 1);

  // fmpz_poly_set_coeff_si(f, 0, 617);
  // fmpz_poly_set_coeff_si(f, 1, 916);
  // fmpz_poly_set_coeff_si(f, 2, 880);
  // fmpz_poly_set_coeff_si(f, 3, 663);
  // fmpz_poly_set_coeff_si(f, 4, 216);
  fmpz_poly_set_coeff_si(f, 0, 21);
  fmpz_poly_set_coeff_si(f, 1, -9);
  fmpz_poly_set_coeff_si(f, 2, -4);
  fmpz_poly_set_coeff_si(f, 4, 5);
  fmpz_poly_set_coeff_si(f, 6, 3);
  
  fmpz_poly_t subresultant_polys[7];

  for (int i = 0; i < 7; i++)
    fmpz_poly_init(subresultant_polys[i]);

  fmpz_subresultant_polynomials_pseudo_remainder_sequence(subresultant_polys,  g, f);
  for (int i = 0; i < 7; i++) {
    fmpz_poly_print_pretty(subresultant_polys[i], "x");
    printf("\n");
  }

  return 0;
}