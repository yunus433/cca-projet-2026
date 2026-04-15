#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>


void fmpz_poly_normalize(fmpz_poly_t poly) {
    if (fmpz_poly_is_zero(poly))
        return;

    fmpz_t content;
    fmpz_init(content);

    // Compute GCD of coefficients
    fmpz_poly_content(content, poly);

    // Divide by content if not 1
    if (!fmpz_is_one(content))
        fmpz_poly_scalar_divexact_fmpz(poly, poly, content);

    // Make leading coefficient positive
    if (fmpz_sgn(fmpz_poly_lead(poly)) < 0)
        fmpz_poly_neg(poly, poly);

    fmpz_clear(content);
}