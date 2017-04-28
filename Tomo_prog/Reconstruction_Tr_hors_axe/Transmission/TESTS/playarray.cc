#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

using namespace std;



/*
  Ce code permet de transférer une matrice (stockage par tableau 1D) dans une matrice plus grande
  recopie centrée, comblage des vides par répétitions des valeurs du centre
  Il est facile de modifier les valeurs mises pour combler les côtés
 */



#define ARRAY_MATRIXDISPLAY(array, lines, cols)		\
  {							\
    size_t i, j;					\
    for (i = 0; i < lines; i++)				\
      {							\
	printf("\n");					\
	for (j = 0; j < cols; j++)			\
	  printf("%d | ", array[i * lines + j]);	\
      }							\
  }




int 
main(void)
{
  int* TS;
  int* TD;
  int NS = 4; // côté du carré
  int SS = NS * NS;
  int ND = 8;
  int SD = ND * ND;
  int delta = (ND - NS) / 2; // doivent être pairs, nécessairement (delta uniforme en haut/bas/gauche/dte)
  cout << "delta:" << delta;

  ARRAY_ALLOC(TS, SS, int);
  ARRAY_ALLOC(TD, SD, int);

  for (int i = 0; i < SS; i++)
    TS[i] = i + 1;
  for (int i = 0; i < SD; i++)
    TD[i] = 0;
  
  cout << endl << "source:";
  ARRAY_MATRIXDISPLAY(TS, NS, NS);
  cout << endl << endl;
  ARRAY_MATRIXDISPLAY(TD, ND, ND);
  cout << endl << endl;
  
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------


  // TS est une matrice 4*4 contenant des nombres de 1 à 16, et représentée par un tableau 1D
  // TD fait 6 par 6: il faut le remplir avec TS au centre, et fixer des valeurs pour les gouttières latérales et haut/bas

  int i, j, run;
  int nsd = NS + delta;
  int nsdm1 = nsd - 1;
  int ns2d = nsd + delta;

  // curseurs posés sur le début du tableau source (petit) et du destination (plus grand)
  int cS = 0; int cD = 0;

  
  // on parcourt selon le tableau grand à remplir
  for (i = 0; i < ND; i++)
    {
      run = 0;

      // on saute les lignes du haut et du bas, on commence par le centre
      if (i < delta) 
	{ cD += ND; continue; }
      if (i > nsdm1)
	{ cD += ND; continue; }

      // on est sur une ligne du centre: il y a delta valeurs à interpoler avant et après
      // ----------------------------------------
      // avant:
      for (j = 0, run = delta -1; \
	   j < delta;		  \
	   j++, cD++, run--)
	{ TD[cD] = TS[cS + run]; }
	  
      // ----------------------------------------
      // centre de la ligne: on recopie les valeurs
      for (j = delta; j < nsd; \
	   j++, cD++, cS++)
	{ TD[cD] = TS[cS]; }

      cS--;
      
      // ----------------------------------------
      // fin de la ligne: 
      for (j = nsd, run = 0; \
	   j < ns2d;		    \
	   j++, cD++)
	{ TD[cD] = TS[cS - run]; }
     
      cS++;
    }


  // maintenant, il faut gérer les lignes du haut et du bas.
  // on se contente de répéter les lignes déjà écrites dans TD, TS n'est plus nécessaire
  cD = 0; 
  int cS2 = delta * ND; // curseur source sur la première ligne affectée

  for (i = 0; i < delta; i++)
    {
      for (j = 0; j < ND; j++, cD++, cS2++)
	{ TD[cD] = TD[cS2]; }
      cS2 -= ND; // on répète indéfiniment la première ligne affectée
    }	     

  cD = (delta + NS) * ND;
  int cS3 = (delta + NS - 1) * ND;
  for (i = 0; i < delta; i++)
    {
      for (j = 0; j < ND; j++, cD++, cS3++)
	{ TD[cD] = TD[cS3]; }
      cS3 -= ND;
    }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  
  cout << endl << "resultat: " << endl;
  cout << endl;
  ARRAY_MATRIXDISPLAY(TD, ND, ND);
  
  cout << endl << endl;
  
  return EXIT_SUCCESS;
}
