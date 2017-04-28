#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "volumeCopy.h"

using namespace std;



// ==============================================================================
//
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
//
//
// ==============================================================================



void
copy_3Dmatrix_in_biggerC2C(double* srcR, double* srcI, cufftComplex* dst, size_t src_edge, size_t dst_edge)
{
  int delta = (dst_edge - src_edge) / 2;

  int i, j, k, run;
  int nsd = src_edge + delta;
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
	    { dst[cD].x = srcR[cS + run]; 
	      dst[cD].y = srcI[cS + run]; 
	    }
	  
	  // ----------------------------------------
	  // centre de la ligne: on recopie les valeurs
	  for (j = delta; \
		 j < nsd;			\
	       j++, cD++, cS++)
	    { dst[cD].x = srcR[cS]; 
	      dst[cD].y = srcI[cS]; 
	    }

	  cS--;
      
	  // ----------------------------------------
	  // fin de la ligne: 
	  for (j = nsd, run = 0;	    \
	       j < ns2d;		    \
	       j++, cD++)
	    { dst[cD].x = srcR[cS - run]; 
	      dst[cD].y = srcI[cS - run]; 
	    }
     
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
	    { dst[cD2].x = dst[cS2].x; 
	      dst[cD2].y = dst[cS2].y; 
	    }
	  cS2 -= dst_edge; // on répète indéfiniment la première ligne affectée
	} 


      cD3 = clD + (delta + src_edge) * dst_edge;
      cS3 = cD3 - dst_edge;
      for (i = 0; i < delta; i++)
	{
	  for (j = 0; j < dst_edge; j++, cD3++, cS3++)
	    { dst[cD3].x = dst[cS3].x; 
	      dst[cD3].y = dst[cS3].y;
	    }
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
//
// ==============================================================================





void
copy_3DmatrixC2C_in_smaller(cufftComplex* big, double* smallR, double* smallI, size_t big_edge, size_t small_edge)
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
	      smallR[cS] = big[cB].x;
	      smallI[cS] = big[cB].y;
	      cB++; cS++;
	    }
	  cB += delta;
	}

      cB += delta * big_edge;
    }
  // ***************************************************************************

}







