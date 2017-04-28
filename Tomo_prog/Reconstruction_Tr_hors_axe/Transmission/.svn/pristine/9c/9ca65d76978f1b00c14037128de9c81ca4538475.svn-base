#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>

using namespace std;



#include "main.h"
//#include "fourier.h"
#include "memory.h"
#include "util.h"


#ifdef CUDA
#include "cuFourier.h"
#endif




// bientot complètement obsolet, en cours de remplacement complet


// cette fonction modifie reel_arc et imag_arc de façon à les remplacer par leur transformée inverse de fourier.
// l'idée sous-jacente est d'économiser de la mémoire vive compte-tenu de la place que prend un seul volume.
// un seul volume temporaire est employé comme tampon pour les permutations de circshift.
void
fourier3DInverse_cshift(double* reel_arc, double* imag_arc, int tab_size, int cube_edge, bool doCuda)
{
  // *******************************************************
  // circshift avant invocation

  double *reel_arc_shift, *imag_arc_shift;  
  double *temp, *temp2;


  ARRAY_NEW(temp, tab_size, double);
  circshift3D_memcpy(reel_arc, temp, cube_edge, cube_edge, cube_edge);
  reel_arc_shift = temp;

  temp2 = reel_arc; 
  circshift3D_memcpy(imag_arc, temp2, cube_edge, cube_edge, cube_edge);
  imag_arc_shift = temp2;

  memcpy(imag_arc, temp, tab_size * sizeof(double)); // dest, src, size
  reel_arc_shift = imag_arc; delete temp;


  // on n'a pas besoin de conserver les versions shiftées et non-shiftées en même temps.
  // de ce fait, on écrase reel_arc et imag_arc, mais pour limiter les permutations, reel_arc_shift se retrouve dans imag_arc et réciproquement


  // *******************************************************
  // TF3D inverse

  
  printf("*******************************************\n");
  fprintf(stdout, "TF3D%s \n", (doCuda) ? "avec Cuda" : "avec fftw");
  printf("*******************************************\n");
  fingers();
  
  
TIMES("TF3D Inverse", {                
    
#ifdef CUDA
    if (doCuda)    
      TF3D_cuda(cube_edge, cube_edge, cube_edge, reel_arc_shift, imag_arc_shift, false);
      //TF3D_cuda_cubepower(cube_edge, 512, reel_arc_shift, imag_arc_shift);
    else
#endif
      // TF3D(cube_edge, reel_arc_shift, imag_arc_shift, FFTW_BACKWARD);      
      // obsolet
  });
 
 
  // *******************************************************
  // circshift après invocation


  ARRAY_NEW(temp, tab_size, double);
  memcpy(temp, reel_arc_shift, tab_size * sizeof(double)); // dest, src, size
  circshift3D_memcpy(imag_arc_shift, imag_arc, cube_edge, cube_edge, cube_edge);
  circshift3D_memcpy(temp, reel_arc, cube_edge, cube_edge, cube_edge);
  delete temp;
}
