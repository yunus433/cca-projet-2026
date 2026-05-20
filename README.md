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

À faire (15.05.2026):
- fixer tous les bugs (subresultant pseudo remainder, cas avec 0 aussi)
- Dire que subresultant algo c'est le même pas besoin de pseudo code
- Expliquer mieux la fluctuations avec une exemple de FLINT random
- Ajouter pourcentange de différence pour les benchmarks
- Faire resultant pour n=3 
- Pour le meilleur algo de resultant, comparer avec wolfram benchmark (temps)
- Rapport:
        - Intro change un peu, fixer les benchmarks etc.
        - Sylvester et Proposition 5 pas besoin, peut faire comme un rappel
        - chapitre 2 ajoute subresultant peut etre aussi trouve avec algo
        - conclusion: resultat theorique, resultats pratiques, etc..
- Presentation:
        - Avoir euclid coeff / subres coeff cote à cote
        - subresultant polynomial exemple
        - subresultant algo change un peu, le faire plus facile avec C1 C2 par ex
        - concataner les slides (6-7 slides) et avoir pause dans Latex
        - state of art - date etc..  Dire implementation est sur FLINT
        - partie theorique ajouter multimodulaire avec chinese remainder the
        - Dire que le GCD est le dernère polynome sousresultant non-bul
        - (D + 1)^(n-1)
        - Avoir les commentaires des benchmarks sur des boites sur des images
        - Changer images sur le meme diapos (le nb de page ne doit pas chagner!!!)
        - conclusion

13 mercredi => slides etc..
21 jeudi => repetition

SLIDES: Après la vacance

BONUS:
- sur flint: comment faire les calculs sur un corps fini? On veut definir des polnomes sur un corps finis, pas Z
             on va faire des calcul multi-modulaires
        - pas besoin de faire, s'il y a une question on dit que tous les implem de FLINT sont multi-mod donc le notre aussi. 

Presentation:
15 mins
