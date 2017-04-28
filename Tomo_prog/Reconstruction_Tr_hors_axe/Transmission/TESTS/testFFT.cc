#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <Magick++.h>
#include <fftw3.h>
#include <cstring>
#include <fstream>
#include <assert.h>
#include <exception>

using namespace std;

#include "main.h"
#include "util.h"
#include "readloop.h"
#include "fourier.h"
#include "memory.h"
#include "IO.h"

#ifdef CUDA
#include "cuFourier.h"
#endif


/* -------------------------------------------------------------------------- */
// 
/* -------------------------------------------------------------------------- */

// nombre de threads employés sous fftw en TF3D
int g_fftw_threads = 4;






char* g_MASKPATH;
char* g_OUTPUT_DIR;
char* g_INPUT_DIR;
char* g_INFILE_RADIX; 



/* -------------------------------------------------------------------------- */
// Usage 
/* -------------------------------------------------------------------------- */



int main(int argc, char *argv[])
{
  int Nxmax_Rf;
  Nxmax_Rf = atoi(argv[1]);
  
  int N_tab=64*Nxmax_Rf*Nxmax_Rf*Nxmax_Rf;


  
  // ==========================================================
  // allocation des volumes de calcul
  // ==========================================================



  


  
  // initialize multithreaded fftw 
  if (g_fftw_threads > 1)
    assert(fftw_init_threads());
  

  

  double *reel_arc_shift, *imag_arc_shift;      
  /////////////////////creation des tableaux qui recoivent le circshift
  //partie reelle du papillon dans l'espace reciproque  //
  ARRAY_NEW(reel_arc_shift, N_tab, double);
  
  // partie imaginaire du papillon dans l'espace reciproque //
  ARRAY_NEW(imag_arc_shift, N_tab, double);


  // chargement des données
  // *******************************************************
  // circshift apres TF3D

  
  double *final_reel_shift, *final_imag_shift;
  try
    {
      //partie reelle du resultat final shifté (centré)
      final_reel_shift=new double[N_tab];
      //partie imaginaire du resultat final shifté (centré)
      final_imag_shift=new double[N_tab];
    }
  catch(std::bad_alloc& ex)
    {cerr << endl << "memory allocation failure" << ex.what(); return EXIT_FAILURE; } 
  
  
  TIMES("Circshift retour",
	{
	  circshift3D_memcpy(final_reel, final_reel_shift, dv0rf, dv0rf, dv0rf);
	  circshift3D_memcpy(final_imag, final_imag_shift, dv0rf, dv0rf, dv0rf);
	  delete final_reel;
	  delete final_imag;
	});


 
  
  // *******************************************************
  // TF3D
  
  printf("*******************************************\n");
  printf("TF3D\n");
  printf("*******************************************\n");
  fingers();


  // à ce stade, le volume ne fait plus que 508^3, soit 4*127 de large. 
  // Il y a une ligne de rognée en 128*128 qui se répercute
  

  TIMES("TF3D fftw", {    
      TF3D(N_tab, 4 * Nxmax_Rf, reel_arc_shift, imag_arc_shift);      
    });  
  TIMES("TF3D cuda", {    
      TF3D_cuda(4 * Nxmax_Rf, 4 * Nxmax_Rf, 4 * Nxmax_Rf, reel_arc_shift, imag_arc_shift);
    });







  //libération memoire allouée pour les threads
  void fftw_cleanup_threads(void);
  


  return EXIT_SUCCESS;
}


