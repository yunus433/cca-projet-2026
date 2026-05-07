# cca-projet-2026

- fixer le naive sousresultant polynomials
- fixer les signes et dernier coeff de flint sous resultant poly
- ajoute methode interpolation pour le sousresultant des polybivariee SANS la vérification
        - comparer ça avec le FLINt sousresultant 
- trouver le methode FLINT bivariée pour les polynomials sousresultant 
- trouver le methode FLINT bivariée pour les sousresultant 
        1. trouver l'implementation
        2. generer des polynomes sur Z mais les traites comme Q (anneau générique)
- ajouter les benchmarks sur le rapport

Mr. Wang
- Ecrire comment faire l'interpolation, mieux expliquer resultant / sous-resultant

À Faire:
- Rapport: Ajouter l'algo, les explications, etc..
        Details de l'implementation
        Multimodulaire
        Sous resultant aussi, mais pais sous resultant poly
- Benchmarks: 
        length 32, inc 2
        exp bitsize 6-7 
        coeff: 50
        count: 5
- Expliquer les résultats
- FMPZ subresultant => implementation subresultant
- Essayer d'expliquer la condition d'élimination de sousresultants aussi

Presentation:
- un peu example (pseudo remainder aussi si on veut)
- pas d'autres algo, seulement methode avec interpolation
- les benchmarks

13 mercredi => slides etc..
21 jeudi => repetition

SLIDES: Après la vacance

BONUS:
- sur flint: comment faire les calculs sur un corps fini? On veut definir des polnomes sur un corps finis, pas Z
             on va faire des calcul multi-modulaires
        - pas besoin de faire, s'il y a une question on dit que tous les implem de FLINT sont multi-mod donc le notre aussi. 

Presentation:
15 mins
