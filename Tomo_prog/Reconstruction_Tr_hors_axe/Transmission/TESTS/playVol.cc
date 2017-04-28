#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

using namespace std;



/*
  suite de playArray
  Cette fois-ci, on fait la même chose avec des matrices 3D stockées comme tableaux 1D
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


#define ARRAY_VOLUMEDISPLAY(array, lines, cols, stacks)		\
  {								\
  int* l_array = array;						\
  size_t k;							\
  								\
  for (k = 0; k < stacks; k++)					\
    {                                                   \ 
printf("\nstack: %d", k);				\
ARRAY_MATRIXDISPLAY(l_array, lines, cols);		\
l_array += lines * cols;				\
}							\
}							\




// ==============================================================================
//
// ==============================================================================


// works for src > dst, src and dst even, and of course src and dst isocubes (we give edge lengthes)
template <typename T>
void
copy_3Dmatrix_in_bigger(T* src, T* dst, size_t src_edge, size_t dst_edge);


// conversely,
template <typename T>
void
copy_3Dmatrix_in_smaller(T* src, T* dst, size_t src_edge, size_t dst_edge);


// ==============================================================================
//
// ==============================================================================



template <typename T>
void
copy_3Dmatrix_in_bigger(T* src, T* dst, size_t src_edge, size_t dst_edge)
{
  int delta = (dst_edge - src_edge) / 2;

  int i, j, k, run;
  int nsd = src_edge + delta;
  int nsdm1 = nsd - 1;
  int ns2d = nsd + delta;

  printf("\n delta: %d / nsd: %d / nsdm1: %d / ns2d: %d \n", delta, nsd, nsdm1, ns2d);

  // curseurs posés sur le début du tableau source (petit) et du destination (plus grand)
  int cS = 0; int cD = 0;
  // pour les sauvegardes des curseurs à chaque début de couche 2D 
  int clS, clD; 
  // pour les parcours locaux
  int cS2, cS3, cD2, cD3;


  
  // ***************************************************************************
  // ***************************************************************************
  for (k = delta; k < nsd; k++)
    {

      
      cD = k * dst_edge * dst_edge; // d'office, on place le curseur de destination en début de couche puisqu'on parcourt selon le volume de destination
      // sauvegardes locales des valeurs de cS et cD
      clS = cS; clD = cD;

      // MATRICE
      // on parcourt selon le tableau grand à remplir
      for (i = 0; i < dst_edge; i++)
	{
	  run = 0;

	  // on saute les lignes du haut et du bas, on commence par le centre
	  if (i < delta) 
	    { cD += dst_edge; continue; }
	  if (i > nsdm1)
	    { cD += dst_edge; continue; }

	  // on est sur une ligne du centre: il y a delta valeurs à interpoler avant et après
	  // ----------------------------------------
	  // avant:
	  for (j = 0, run = delta -1; \
	       j < delta;	      \
	       j++, cD++, run--)
	    { dst[cD] = src[cS + run]; }
	  
	  // ----------------------------------------
	  // centre de la ligne: on recopie les valeurs
	  for (j = delta; \ 
		 j < nsd;			\
	       j++, cD++, cS++)
	    { dst[cD] = src[cS]; }

	  cS--;
      
	  // ----------------------------------------
	  // fin de la ligne: 
	  for (j = nsd, run = 0;	    \
	       j < ns2d;		    \
	       j++, cD++)
	    { dst[cD] = src[cS - run]; }
     
	  cS++;
	}


      // maintenant, il faut gérer les lignes du haut et du bas.
      // on se contente de répéter les lignes déjà écrites dans TD, TS n'est plus nécessaire.

      // on repointe en début de matrice pour la coupe courante
      cD2 = clD; 
      cS2 = clD + delta * dst_edge;
      for (i = 0; i < delta; i++)
	{
	  for (j = 0; j < dst_edge; j++, cD2++, cS2++)
	    { dst[cD2] = dst[cS2]; }
	  cS2 -= dst_edge; // on répète indéfiniment la première ligne affectée
	} 


      cD3 = clD + (delta + src_edge) * dst_edge;
      cS3 = cD3 - dst_edge;
      for (i = 0; i < delta; i++)
	{
	  for (j = 0; j < dst_edge; j++, cD3++, cS3++)
	    { dst[cD3] = dst[cS3]; }
	  cS3 -= dst_edge;
	}

      
    }

  // ***************************************************************************
  // ***************************************************************************
  // couches du début

  // couches de fin

}





// ==============================================================================
//
// ==============================================================================



template <typename T>
void
copy_3Dmatrix_in_smaller(T* big, T* small, size_t big_edge, size_t small_edge)
{
  int delta = (big_edge - small_edge) / 2;
  int i, j, k, run;
  
  // curseurs posés sur le début du tableau source (petit) et du destination (plus grand)
  int cS = 0; int cB = 0;


  cB = delta * big_edge * big_edge; // on saute les premières couches du volume destination

  // ***************************************************************************
  for (k = 0; k < small_edge; k++)  
    {
      cB += delta * big_edge;

      // MATRICE
      // on parcourt selon le tableau grand 
      for (i = 0; i < small_edge; i++)
	{
	  cB += delta;
	  for (j = 0; j < small_edge; j++)
	    {
	      small[cS] = big[cB];
	      cB++; cS++;
	    }
	  cB += delta;
	}

      cB += delta * big_edge;
    }
  // ***************************************************************************

}








// ==============================================================================
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
// ==============================================================================



int 
main(void)
{
  int* TS;
  int* TD;
  int* TF;
  int NS = 4; // côté du carré
  int SS = NS * NS;
  int SMS = SS * NS;
  int ND = 8;
  int SD = ND * ND;
  int SMD = SD * ND;
  int delta = (ND - NS) / 2; // doivent être pairs, nécessairement (delta uniforme en haut/bas/gauche/dte)
  cout << "delta:" << delta;

  ARRAY_ALLOC(TS, SMS, int);
  ARRAY_ALLOC(TF, SMS, int);
  ARRAY_ALLOC(TD, SMD , int);

  for (int i = 0; i < SMS ; i++)
    TS[i] = i + 1;
  for (int i = 0; i < SMD; i++)
    TD[i] = 0;
  
  cout << endl << "source:";
  ARRAY_VOLUMEDISPLAY(TS, NS, NS, NS);
  cout << endl << endl;
  ARRAY_VOLUMEDISPLAY(TD, ND, ND, ND);
  cout << endl << endl;
  

  
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  
  copy_3Dmatrix_in_bigger(TS, TD, NS, ND);

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  copy_3Dmatrix_in_smaller(TD, TF, ND, NS);
  
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  
  cout << endl << "resultat: " << endl;
  cout << endl;
  ARRAY_VOLUMEDISPLAY(TD, ND, ND, ND);
  
  cout << endl << endl;

  cout << endl << "resultat fenêtrage: " << endl;
  cout << endl;
  ARRAY_VOLUMEDISPLAY(TF, NS, NS, NS);

  
  return EXIT_SUCCESS;
}

















/*
//exit(1);

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// TS est une matrice 4*4 contenant des nombres de 1 à 16, et représentée par un tableau 1D
// TD fait 6 par 6: il faut le remplir avec TS au centre, et fixer des valeurs pour les gouttières latérales et haut/bas

int i, j, k, run;
int nsd = NS + delta;
int nsdm1 = nsd - 1;
int ns2d = nsd + delta;

// curseurs posés sur le début du tableau source (petit) et du destination (plus grand)
int cS = 0; int cD = 0;
// pour les sauvegardes des curseurs à chaque début de couche 2D 
int clS, clD; 
// pour les parcours locaux
int cS2, cS3, cD2, cD3;



// ***************************************************************************
// ***************************************************************************
for (k = delta; k < nsd; k++)
{

      
cD = k * SD; // d'office, on place le curseur de destination en début de couche puisqu'on parcourt selon le volume de destination
// sauvegardes locales des valeurs de cS et cD
clS = cS; clD = cD;

// MATRICE
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
j < delta;	      \
j++, cD++, run--)
{ TD[cD] = TS[cS + run]; }
	  
// ----------------------------------------
// centre de la ligne: on recopie les valeurs
for (j = delta; \ 
j < nsd;			\
j++, cD++, cS++)
{ TD[cD] = TS[cS]; }

cS--;
      
// ----------------------------------------
// fin de la ligne: 
for (j = nsd, run = 0;	    \
j < ns2d;		    \
j++, cD++)
{ TD[cD] = TS[cS - run]; }
     
cS++;
}


// maintenant, il faut gérer les lignes du haut et du bas.
// on se contente de répéter les lignes déjà écrites dans TD, TS n'est plus nécessaire

// on repointe en début de matrice pour la coupe courante
cD2 = clD; 
cS2 = clD + delta * ND;

for (i = 0; i < delta; i++)
{
for (j = 0; j < ND; j++, cD2++, cS2++)
{ TD[cD2] = TD[cS2]; }
cS2 -= ND; // on répète indéfiniment la première ligne affectée
}	     

cD3 = clD + (delta + NS) * ND;
cS3 = cD3 - ND;
for (i = 0; i < delta; i++)
{
for (j = 0; j < ND; j++, cD3++, cS3++)
{ TD[cD3] = TD[cS3]; }
cS3 -= ND;
}

      
}

// ***************************************************************************
// ***************************************************************************
// couches du début

// couches de fin


  
// ***************************************************************************
// ***************************************************************************
*/
