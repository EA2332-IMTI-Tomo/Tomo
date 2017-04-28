#include <omp.h>
#include <iostream>
#include <cstdlib>

#include <cuda.h>
#include <cuda_runtime.h>

#include "vChrono.h"

#include "macros.h"

using namespace std;


// --------------------------------------------------
// --------------------------------------------------

inline float
nozero(float f)
{
  return (f) ? f : 1; 
}

// --------------------------------------------------
// --------------------------------------------------

__global__ void 
MOYENNE_d(unsigned short int* counts, float *tr, float *ti, size_t tabsize, size_t *calls)
{
  // on infère du bloc (sous-division du kernel) et du thread 
  // (atome d'exécution) l'indice unique pour laquelle cette fonction
  // est invoquée par le GPU

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int val;
  if (i < tabsize) 
    {
      val = counts[i];
      val = (val) ?  val : 1;

      tr[i] /= val;
      ti[i] /= val;
      *calls = (*calls)+1;
    }
}


__global__ void 
MOYENNE2_d(unsigned short int* counts, float *tr, size_t tabsize, size_t *calls)
{
  // on infère du bloc (sous-division du kernel) et du thread 
  // (atome d'exécution) l'indice unique pour laquelle cette fonction
  // est invoquée par le GPU

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int val;
  if ((i < tabsize) && counts[i])
    {
      tr[i] /= counts[i];
      //*calls = (*calls)+1;
    }
}
// --------------------------------------------------
// --------------------------------------------------

int 
main(int argc, char **argv)
{
  int GRIDSIZE = 50;
  if (argc - 1)
    GRIDSIZE = atoi(argv[1]);
  
  cout << endl << GRIDSIZE << endl;


  size_t i = 0;
  size_t tabsize = 256 * 256 * 256;
  float supRedSafe;

  size_t calls = 0;
  
  unsigned short int *sup_redon;
  float *reel_arc, *imag_arc;
  ARRAY_ALLOC(sup_redon, tabsize, unsigned short int);
  ARRAY_ALLOC(reel_arc, tabsize, float);
  ARRAY_ALLOC(imag_arc, tabsize, float);

  // pour test
  float *reel_arc_comp, *reel_arc_err;
  ARRAY_ALLOC(reel_arc_comp, tabsize, float);
  ARRAY_ALLOC(reel_arc_err, tabsize, float);

  unsigned short int *sup_redon_d;
  float *reel_arc_d, *imag_arc_d;

  
  cudaMalloc( (void**) &sup_redon_d, tabsize * sizeof( unsigned short int ));
  cudaMalloc( (void**) &reel_arc_d, tabsize * sizeof( float ));
  cudaMalloc( (void**) &imag_arc_d, tabsize * sizeof( float ));


  for (i = 0; i < tabsize; i++)
    { 
      sup_redon[i] = 5 + (i % 256) + (2 * i) % 84;
      reel_arc[i] = 0.32f + sup_redon[i] / 1.33f ;
      imag_arc[i] = reel_arc[i] * 0.89f * float(sup_redon[i]);
    }



  vCHRONO_SET(timer1, "toto1 sur cpu");
  vCHRONO_START(timer1);

  for(int i = 0; i < tabsize; i++)
    {
      if (sup_redon[i])
	{
	  reel_arc[i] /= sup_redon[i];
	  imag_arc[i] /= sup_redon[i];
	}
    }

  /*
    supRedSafe = sup_redon[i];
    supRedSafe = (supRedSafe)? supRedSafe : 1;
    // accès non protégé à VC_Volume / V_RealPart / V_ImagPart
    reel_arc[i] /= supRedSafe;
    imag_arc[i] /= supRedSafe;
    }
    // 1.59
    /*

    vCHRONO_STOP(timer1);

  */

  cudaMemcpy(sup_redon_d, sup_redon, tabsize * sizeof(unsigned short int), cudaMemcpyHostToDevice);
  cudaMemcpy(reel_arc_d, reel_arc, tabsize * sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(imag_arc_d, imag_arc, tabsize * sizeof(float), cudaMemcpyHostToDevice);


  /* dimensionnement de la "grille" de calcul */
  int blocksize = 512; // nombre de threads par bloc
  dim3 dimBlock( blocksize ); 
  dim3 dimGrid( ceil(float (GRIDSIZE) / float(dimBlock.x)));

  
  vCHRONO_SET(timer2, "toto2 sur gpu");
  vCHRONO_START(timer2);

  //   MOYENNE_d<<<dimGrid, dimBlock>>>(sup_redon_d, reel_arc_d, imag_arc_d, tabsize, &calls);
  MOYENNE2_d<<<16, 64>>>(sup_redon_d, reel_arc_d, tabsize, &calls);
  cudaThreadSynchronize();
  vCHRONO_STOP(timer2);


  cudaMemcpy(reel_arc_comp, reel_arc_d, tabsize * sizeof(float), cudaMemcpyDeviceToHost);
  cudaThreadSynchronize();
  
  bool notify = true;
  
  for (size_t i = 0; i < tabsize; i++)
    {
      reel_arc_err[i] = reel_arc_comp[i] - reel_arc[i];

      if (reel_arc_err[i] && notify)
	{
	  cout << endl << "error at: " << i << " of "
	       << reel_arc_err[i] << ": " << reel_arc_comp[i] << " - " << reel_arc[i];
	}
    }
  cout << endl << calls << "calls";

  return EXIT_SUCCESS;
}


// descendent tous les 3 à 1s avec -O3 (4)
// apparemment, omp ne marche pas
