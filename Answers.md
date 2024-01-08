# Answers : Exercise 3

1) Pour pouvoir utiliser une matrice avec LAPACK/BLAS nous devons la déclarer comme étant un pointeur du type désiré : float, int, double, etc. Ensuite nous devons choisir ses dimensions m et n et initialiser les valeurs (m et n dépendent du format de stockage que l'on a choisit pour la matrice) puis allouer la mémoire nécessaire avec un malloc/calloc. Enfin il faut ensuite la remplir a l'aide d'une boucle.

2) Lapack_row_major définie la manière dont sont stockées les valeures de la matrice pour que les méthodes de lapack/blas puissent l'utiliser correctement.

3) La leading dimension sert à savoir comment son organisés les données en mémoire, c'est a dire pour un stockage en ligne par exemple, la leading dimension correspond à longueur d'une colonne (nombre d'éléments séparant 2 lignes). Cela sert a avoir un accès aux données correcte et efficace.

4) DGBMV effectue une multiplication matrice-vecteur sur une matrice stockée en format général band.

5) DGBTRF effectue une factorisation LU sur une matrice bande.

6) DGBTRS résous un systeme linéaire de la forme Ax = b ou A a été factorisé en sa "version" LU;

7) DGBSV résous un system d'equations lineaires AX = B ou X et B sont des matrices de dimensions N par NRHS. A sera factorisé par décomposition LU directement dans la méthode.

8) Avec BLAS pour calculer la norme du résidu relatif on peut, une fois que l'on a obtenu le résultat "b" de notre systeme (Ax = b), 
calculer le résidu à l'aide d'une daxpy puis calculer la norme du résidu en utilisant la fonction cblas_dnrm2.

## Exercice 5 :

2) Pour évaluer les performances en temps nous utilisons ici la fonction clock_gettime.