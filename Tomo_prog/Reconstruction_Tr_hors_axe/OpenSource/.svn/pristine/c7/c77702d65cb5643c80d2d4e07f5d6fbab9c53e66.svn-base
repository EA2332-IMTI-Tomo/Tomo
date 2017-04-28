#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <string.h>
#include <fstream>

#include "macros.h"
#include "vParse.h"

using namespace std;


// =============================================================================


// analyse le fichier texte donné en entrée et détermine le nombre de lignes et de colonnes
// si le nombre de colonnes varie d'une ligne à l'autre, retourne faux
// sinon, écrit résultat dans nb_lines et nb_cols
// indiquer dans max_line_size un nombre de caractères suffisant pour contenir une ligne (typiquement 80, mettre 1000 par précaution)

// header_size: indique combien de lignes il faut omettre pour vérifier le critère de carrité (sic). Ces lignes ne sont pas omises pour le comptage de ligne

bool
retrieveFileDimensions(char* filename, size_t max_line_size, size_t &nb_lines, size_t &nb_cols, size_t header_size)
{
  char* read_line;
  ARRAY_ALLOC(read_line, max_line_size, char);
  
  ifstream infile (filename);
  assert (infile.is_open());
    
  nb_lines = 0;
  nb_cols = 0;

  char* token;
  size_t nb_tokens = 0;

  while (! infile.eof() )
    {
      infile.getline(read_line, max_line_size, '\n');
      nb_lines++;

      // guess how many tokens are present in line
      nb_tokens = 0;

      token = strtok(read_line, " \t");
      while(token != NULL)
	{
	  nb_tokens++;
	  //printf("%s\n", token);
	  token = strtok(NULL, " \t");
	}
      
      // assuming all lines have same number of tokens
      if (nb_lines <= 1 + header_size)
	{nb_cols = nb_tokens;}
      else
	{
	  if (nb_tokens != nb_cols) 
	    if (nb_tokens == 0)
	      { nb_lines--; return true;} // fin de fichier
	    else
	      return false;
	}
      

    }
  
  return true;
}



// =============================================================================


// alloue une matrice et lui affecte les données du fichier spécifié
// matrice de taille nb_lines * nb_cols
// si le fichier contient plus de lignes ou colonnes, elles sont ignorées
// la lecture saute offset_lines linges au début et offset_cols colonnes si besoin (commentaires en première ligne par exemple) 

double**
readMatrixFromFile(char* filename, size_t max_line_size, size_t nb_lines, size_t nb_cols, size_t offset_lines, size_t offset_cols)
{
  char* read_line;
  double** matrix;
  ARRAY_ALLOC(read_line, max_line_size, char);
  MATRIX_ALLOC(matrix, nb_lines, nb_cols, double);

  ifstream infile (filename);
  assert (infile.is_open());

  char* token;
  size_t nb_tokens_read = 0;
  size_t nb_lines_read = 0;
  size_t i = 0, j = 0;


  while (! infile.eof()) 
    {
      infile.getline(read_line, max_line_size, '\n');
      nb_lines_read++;

      if (nb_lines_read <= offset_lines) continue;
      if (nb_lines_read - offset_lines > nb_lines) break;
      
      // guess how many tokens are present in line
      j = 0;

      token = strtok(read_line, " \t");
      nb_tokens_read = 1; 
      while(token != NULL)
	{
	  if (nb_tokens_read - offset_cols > nb_cols) break;
	  if (nb_tokens_read > offset_cols) 
	    {
	      //printf("%s\n", token);
	      matrix[i][j] = atof(token);
	      j++;
	    }
	  token = strtok(NULL, " \t");
	  nb_tokens_read++;
	}
      
    
      i++;
    }
  
  return matrix;
}
