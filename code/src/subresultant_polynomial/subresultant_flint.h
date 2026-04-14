#ifndef SUBRESULTANT_FLINT_H
#define SUBRESULTANT_FLINT_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_poly.h>
#include <flint/fmpz_poly.h>
#include "../utils/utils.h"

void fmpz_poly_pseudo_rem_cohen_2(
    fmpz_poly_t R,
    const fmpz_poly_t A,
    const fmpz_poly_t B
);

// 🔥 suite des sous-résultants (PRS)
int fmpz_subresultant_flint(
    fmpz_poly_t *res,
    const fmpz_poly_t P,
    const fmpz_poly_t Q
);

#endif