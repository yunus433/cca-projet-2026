#ifndef FMPZ_MPOLY_SUBRESULTANT_FLINT_H
#define FMPZ_MPOLY_SUBRESULTANT_FLINT_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

/*
    Calcul d'une suite minimale de sous-résultants
    bivariés style FLINT.

    subresultants :
        tableau déjà initialisé de fmpz_mpoly_t

    len :
        nombre de sous-résultants calculés

    var :
        variable principale
*/

int
fmpz_mpoly_subresultant_sequence(
    fmpz_mpoly_t * subresultants,
    slong * len,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    slong var,
    const fmpz_mpoly_ctx_t ctx);

#endif