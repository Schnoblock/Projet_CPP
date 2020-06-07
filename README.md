#Projet C++

Lucas Perrin

Projet de C++ effectué pour le second semestre du Master 1 Mathématiques appliquées.

Remarques sur l'ensemble du projet :

- Projet réalisé sur repl.it / Atom (sur mac).
- Aussi mis en ligne sur GitHub.

Etat du projet :

- Partie 1,2,3,4,5 entièrement codées sauf hosvd
- Partie 1,2 validées
- Partie 3 : qrpivot et svd (dépend de qrpivot) non validées
- Partie 4 : pmod non validée
- Partie 5 : non validée (tenseurtotal dépend de pmod et hosvd dépend de svd)

Points positifs :
- Sujet très intéréssant, apportant plus que la simple connaissance du langage avec des applications qui parlent à des élèves en parcours statistiques.
- Sujet complet vis à vis des particularités du langage, création et manipulation de classes, d'allocation dynamique, de passage par références, etc...

Difficultés rencontrés lors de la réalisation :

- Appréhension des notions d'allocation de mémoire. Cela a notamment amené à des "segmentation error" dues à des doubles libérations de pointeurs à cause de déclarations de variables mal formulées.
- Compréhension de la formule pour la récupération d'indices sur la classe 'Tenseur' (codage des fonctions 'phi()' et 'phi_inv()', ainsi que le constructeur de la classe 'Tenseur' passant par une Matrice, et la fonction membre '.mode' de la classe 'Matrice', et le produit pmod).
- Temps, dans le sens ou la donnée des exames à distances a conduit, personnellement, à une répartition du temps entre les différentes matières à étudier, réviser, ainsi que les différents projets qui a été en partie mal gérée pour un projet aussi long.
- Trouver les erreurs s'est avérée bcp plus technique que sur d'autres langages.
- Compréhension parfois difficile des relations entre classes notamment vis à vis de l'accès aux membres données privées (ou protégées) (problème contourné avec une fonction getdims() pour les matrices (utilisée dans pmod) et rencontré à nouveau sur hosvd.).
