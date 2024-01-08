# Answers :

## Exercise 3

1) Pour pouvoir utiliser une matrice avec LAPACK/BLAS nous devons la déclarer comme étant un pointeur du type désiré : float, int, double, etc. Ensuite nous devons choisir ses dimensions m et n et initialiser les valeurs (m et n dépendent du format de stockage que l'on a choisit pour la matrice) puis allouer la mémoire nécessaire avec un malloc/calloc. Enfin il faut ensuite la remplir a l'aide d'une boucle.

2) Lapack_row_major définie la manière dont sont stockées les valeures de la matrice pour que les méthodes de lapack/blas puissent l'utiliser correctement.

3) La leading dimension sert à savoir comment son organisés les données en mémoire, c'est a dire pour un stockage en ligne par exemple, la leading dimension correspond à longueur d'une colonne (nombre d'éléments séparant 2 lignes). Cela sert a avoir un accès aux données correcte et efficace.

4) DGBMV effectue une multiplication matrice-vecteur sur une matrice stockée en format général band.

5) DGBTRF effectue une factorisation LU sur une matrice bande.

6) DGBTRS résous un systeme linéaire de la forme Ax = b ou A a été factorisé en sa "version" LU;

7) DGBSV résous un system d'equations lineaires AX = B ou X et B sont  des matrices de dimensions N par NRHS. A sera factorisé par décomposition LU directement dans la méthode.

8) Avec BLAS pour calculer la norme du résidu relatif on peut, une fois que l'on a obtenu le résultat "b" de notre systeme (Ax = b), 
calculer le résidu à l'aide d'une daxpy puis calculer la norme du résidu en utilisant la fonction cblas_dnrm2.

## Exercice 5 :

2) Pour évaluer les performances en temps nous utilisons ici la fonction clock_gettime :

    Des graphiques représentant les performances des différentes méthodes de résolution directes sont disponibles dans le
    dossier Benchmarks dans leurs emplacement respectifs.

    D'après les observations des graphiques, on peut observer que :

    - La complexité de DGBSV en temps semble être de l'ordre de O(N²) 

    - La complexité de DGBTRF (implémentation LAPACK) semble être
      de l'ordre de O(N), la croissance de la courbe est linéaire.

    - La complexité de DGBTRS semble aussi de par sa courbe avoir une
      complexité linéaire en O(N). Lors des benchmarks j'ai pu observer
      un pic pour une taille de matrice 300x300. Cela est surement du
      à une analyse peu précise. J'aurais pu faire des répétitions et
      des moyennes de temps d'execution en passant plus de temps sur
      les benchmarks.

## Exercice 7 :

3)  D'après l'analyse du vecteur de résidus de la méthode 
    Richardson Alpha nous observons que la méthode converge vers 0 à 
    partir d'environ 460 itérations. (Le graphique représentant ce 
    résultat est disponible dans le dossier Convergence).

## Exercice 8 :

3)  D'après l'analyse du vecteur de résidus de la méthode de Jacobi
    nous observons que la méthode converge vers 0 à partir 
    d'environ 130 itérations. (Le graphique représentant ce résultat
    est disponible dans le dossier Convergence).

## Exercice 8 :

3)  D'après l'analyse du vecteur de résidus de la méthode de 
    Gauss-Siedel nous observons que la méthode converge vers 0 à 
    partir d'environ 900 itérations. Hors nous devrions observer que 
    cette méthode converge environ deux fois plus vite que celle de 
    Jacobi, cela est probablement du à un bug dans le code dans la 
    manière de calculer une itération de Gauss-Siedel, de plus on 
    remarque que la courbe fluctue énormément ce qui montre bien un 
    comportement anormal. (Le graphique représentant ce 
    résultat est disponible dans le dossier Convergence).