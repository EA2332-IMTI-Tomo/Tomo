#include <omp.h>
#include <iostream>
#include <cstdlib>

#include "vChrono.h" // requiert boost

#include "macros.h"

using namespace std;

#define TYPE double
#define TABSIZE 134217728
// 512 cube

// --------------------------------------------------

inline TYPE
nozero(TYPE f)
{
  return (f) ? f : 1; 
}

// --------------------------------------------------

int 
main()
{
  size_t i = 0;
  size_t tabsize = 512 * 512 * 512;
  TYPE supRedSafe;

  
  unsigned short int *sup_redon;
  TYPE *reel_arc, *imag_arc;
  ARRAY_ALLOC(sup_redon, tabsize, unsigned short int);
  ARRAY_ALLOC(reel_arc, tabsize, TYPE);
  ARRAY_ALLOC(imag_arc, tabsize, TYPE);
  

  // remplissage aléatoire de sup_redon
  for (i = 0; i < tabsize; i++)
    { 
      sup_redon[i] = i % 256;
      reel_arc[i] = sup_redon[i] / 3;
      imag_arc[i] = reel_arc[i] * sup_redon[i];
    }


  //   ==================================================

  //double time0 = omp_get_wtime();
  
  vCHRONO_SET(beginOMP, "avec omp * "); 
  vCHRONO_START(beginOMP);
#pragma omp parallel for num_threads(3)
  for(int i = 0; i < TABSIZE; i++)
    {
      //#pragma omp parallel private(i)
      TYPE n = sup_redon[i];
      TYPE factor = n ? 1.0f / n : 1.0f;
      // accès non protégé à VC_Volume / V_RealPart / V_ImagPart
      reel_arc[i] *= n;
      imag_arc[i] *= n;
    }
  vCHRONO_STOP(beginOMP);
  // 1.59

  //cout << endl << "omp on: elapsed:" << omp_get_wtime() - time0;



  vCHRONO_SET(beginRAW, "sans omp *"); 
  vCHRONO_START(beginRAW);

  for(i = 0; i < TABSIZE; i++)
    {
      //supRedSafe = nozero(sup_redon[i]);
      TYPE n = sup_redon[i];
      TYPE factor = n ? 1.0f / n : 1.0f;

      // accès non protégé à VC_Volume / V_RealPart / V_ImagPart
      reel_arc[i] *= factor;
      imag_arc[i] *= factor; 
    }
  
  vCHRONO_STOP(beginRAW);
  // 1.82


  vCHRONO_SET(beginRAW2, "avec omp lent "); 
  vCHRONO_START(beginRAW2);

  #pragma omp parallel for 
  for(int i = 0; i < TABSIZE; i++)
    {
      supRedSafe = sup_redon[i];
      supRedSafe = (supRedSafe)? 1.0f/supRedSafe : 1.0;
      // accès non protégé à VC_Volume / V_RealPart / V_ImagPart
      reel_arc[i] *= supRedSafe;
      imag_arc[i] *= supRedSafe;
    }
  vCHRONO_STOP(beginRAW2);



    vCHRONO_SET(beginOMP2, "avec omp /"); 
    vCHRONO_START(beginOMP2);
#pragma omp parallel for //num_threads(3) //(dynamic, CHUNKSIZE)
  for(i = 0; i < TABSIZE; i++)
    {
      //#pragma omp parallel private(i)
      if (sup_redon[i])
	{
	  reel_arc[i] /= sup_redon[i];
	  imag_arc[i] /= sup_redon[i];
	}
    }
  vCHRONO_STOP(beginOMP2);
  // 1.31 


  return EXIT_SUCCESS;
}


// descendent tous les 3 à 1s avec -O3 (4)
// apparemment, omp ne marche pas
