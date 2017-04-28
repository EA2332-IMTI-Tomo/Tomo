#ifndef _VECTRA_PARSE_H_
#define _VECTRA_PARSE_H_



// analyse le fichier texte donné en entrée et détermine le nombre de lignes et de colonnes
// si le nombre de colonnes varie d'une ligne à l'autre, retourne faux
// sinon, écrit résultat dans nb_lines et nb_cols
// indiquer dans max_line_size un nombre de caractères suffisant pour contenir une ligne (typiquement 80, mettre 1000 par précaution)

// header_size: indique combien de lignes il faut omettre pour vérifier le critère de carrité (sic). Ces lignes ne sont pas omises pour le comptage de ligne


bool
retrieveFileDimensions(char* filename, size_t max_line_size, size_t &nb_lines, size_t &nb_cols, size_t header_size = 0);


// alloue une matrice et lui affecte les données du fichier spécifié
// matrice de taille nb_lines * nb_cols
// si le fichier contient plus de lignes ou colonnes, elles sont ignorées
// la lecture saute offset_lines linges au début et offset_cols colonnes si besoin (commentaires en première ligne par exemple) 

double**
readMatrixFromFile(char* filename, size_t max_line_size, size_t nb_lines, size_t nb_cols, size_t offset_lines, size_t offset_cols);





#endif // _VECTRA_PARSE_H_
