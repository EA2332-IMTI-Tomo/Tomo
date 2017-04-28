#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>


using namespace std;



#include "main.h"
#include "cuFourier.h"


#include "volumeCopy.h"

// ----------------------------------------------------------------------
// 
// ----------------------------------------------------------------------


void 
TF2D_cuda(double entree_reelle[], double entree_imag[], double fft_reel[], double fft_imag[], int taille_x, int taille_y)
{
  int i, N = taille_x * taille_y;

  cufftHandle plan;
  cufftComplex *h_data, *d_data;

  // allocate array in CPU space
  assert(h_data = (cufftComplex *) calloc(N, sizeof(cufftComplex)));

  // allocate array in GPU memory
  cudaMalloc((void**) &d_data, N * sizeof(cufftComplex));

  
  // set host array to values to process
  for(i = 0 ; i < N ; i++)
    {
      h_data[i].x = entree_reelle[i];
      h_data[i].y = entree_imag[i];
    }
  // and clone it to its GPU twin
  cudaMemcpy(d_data, h_data, N * sizeof(cufftComplex), cudaMemcpyHostToDevice);
  

  // compute GPU FFT (Z2Z: doubleComplex to doubleComplex)
  //   assert( cufftPlan2d(&plan, taille_x, taille_y, CUFFT_C2C) == CUFFT_SUCCESS );
  HANDLE_FFT_ERROR( cufftPlan2d(&plan, taille_x, taille_y, CUFFT_C2C) );
  HANDLE_FFT_ERROR( cufftExecC2C(plan, d_data, d_data, CUFFT_FORWARD) );

  // retrieve GPU results
  cudaMemcpy(h_data, d_data, N * sizeof(cufftComplex), cudaMemcpyDeviceToHost); 

  // transfer host results to parameter arrays
  for(i = 0; i < N; i++)
    {
      fft_reel[i] = h_data[i].x;
      fft_imag[i] = h_data[i].y;
    }
  

  cufftDestroy(plan);
  free(h_data);
  cudaFree(d_data); 

}


// ----------------------------------------------------------------------
// 3D
// ----------------------------------------------------------------------


// *****************************************************************************
// version "standard": TF3D inverse, C2C


void
TF3D_cuda(int NX, int NY, int NZ, double *reel_arc_shift, double *imag_arc_shift, bool forward_or_backward)
{
  int i; 
  int N3D = NX * NY * NZ; 
  
  cufftHandle plan;
  cufftComplex *h_data, *d_data;
  cufftResult fft_result;

  // allocate array in CPU space
  //assert(h_data = (cufftComplex *) calloc(N3D, sizeof(cufftComplex)));
  assert ( cudaHostAlloc(&h_data, N3D * sizeof(cufftComplex), cudaHostAllocDefault) == cudaSuccess );

  // allocate array in GPU memory
  assert ( cudaMalloc((void**) &d_data, N3D * sizeof(cufftComplex)) == cudaSuccess );
  
  // set host array to values to process
  for(i = 0 ; i < N3D ; i++)
    {
      h_data[i].x = reel_arc_shift[i];
      h_data[i].y = imag_arc_shift[i];
    }
  // and clone it to its GPU twin
  assert ( cudaMemcpy(d_data, h_data, N3D * sizeof(cufftComplex), cudaMemcpyHostToDevice) == cudaSuccess );


  // --------------------------------------------------
  // CUFFT core start

  // compute GPU FFT (Z2Z: doubleComplex to doubleComplex)
  assert( cufftPlan3d(&plan, NX, NY, NZ, CUFFT_C2C) == CUFFT_SUCCESS );
  ////assert( cufftExecC2C(plan, d_data, d_data, CUFFT_INVERSE) == CUFFT_SUCCESS );
  
  fft_result = cufftExecC2C(plan, d_data, d_data,   \
			    (forward_or_backward) ? CUFFT_FORWARD : CUFFT_INVERSE);
  if (fft_result != CUFFT_SUCCESS) {
    cufftResult_display(fft_result);
    exit (EXIT_FAILURE);
  }
  
  // CUFFT core end
  // --------------------------------------------------


  // retrieve GPU results
  assert( cudaMemcpy(h_data, d_data, N3D * sizeof(cufftComplex), cudaMemcpyDeviceToHost) == cudaSuccess );

  // transfer host results to parameter arrays
  for(i = 0; i < N3D; i++)
    {
      reel_arc_shift[i] = h_data[i].x;
      imag_arc_shift[i] = h_data[i].y;
    }

  cufftDestroy(plan);
  // free(h_data);
  cudaFreeHost(h_data);
  cudaFree(d_data); 
}




// *****************************************************************************
// version cube: on prend un volume cubique en entrée, on le
// transforme en taille supérieure cubique n^2 pour effectuer la TF
// dans des conditions optimisées, puis on récupère en sortie la même
// taille qu'en entrée

// params: arête du cube du volume d'entrée, arête du cube sur lequel on calcule la TF, puis volumes d'entrée (et sortie) correspondant à N_src
// pour l'instant, on ne gère que des arêtes de taille paire.

void
TF3D_cuda_cubepower(int N_src, int N_tf, double *reel_arc_shift, double *imag_arc_shift, bool forward_or_backward)
{
  int i; 
  int N3D = N_src * N_src * N_src;
  int N3D_tf = N_tf * N_tf * N_tf;
  assert ( N_src < N_tf );
  assert ( N_src % 2 == 0 );
  assert ( N_tf % 2 == 0 );
  int delta = (N_tf - N_src) / 2; 

  cufftHandle plan;
  cufftComplex *h_data, *d_data;
  cufftResult fft_result;

  // allocate array in CPU space
  //assert(h_data = (cufftComplex *) calloc(N3D, sizeof(cufftComplex)));
  assert ( cudaHostAlloc(&h_data, N3D_tf * sizeof(cufftComplex), cudaHostAllocDefault) == cudaSuccess );

  // allocate array in GPU memory
  assert ( cudaMalloc((void**) &d_data, N3D_tf * sizeof(cufftComplex)) == cudaSuccess );
  

  /*
  // set host array to values to process
  // step 1 : fill center window with source values
  for(i = 0 ; i < N3D ; i++)
    {
      h_data[i].x = reel_arc_shift[i];
      h_data[i].y = imag_arc_shift[i];
    }
  */
  copy_3Dmatrix_in_biggerC2C(reel_arc_shift, imag_arc_shift, h_data, N_src, N_tf);

  // step 3: and clone it to its GPU twin
  assert ( cudaMemcpy(d_data, h_data, N3D_tf * sizeof(cufftComplex), cudaMemcpyHostToDevice) == cudaSuccess );


  // --------------------------------------------------
  TIME_START(cufft, "CUFFT core start");

  // compute GPU FFT (Z2Z: doubleComplex to doubleComplex)
  assert( cufftPlan3d(&plan, N_tf, N_tf, N_tf, CUFFT_C2C) == CUFFT_SUCCESS );
  ////assert( cufftExecC2C(plan, d_data, d_data, CUFFT_INVERSE) == CUFFT_SUCCESS );
  
  fft_result = cufftExecC2C(plan, d_data, d_data, \
			    (forward_or_backward) ? CUFFT_FORWARD : CUFFT_INVERSE);
  if (fft_result != CUFFT_SUCCESS) {
    cufftResult_display(fft_result);
    exit (EXIT_FAILURE);
  }
  
  TIME_END(cufft, "CUFFT core end  ");
  // --------------------------------------------------


  // retrieve GPU results
  assert( cudaMemcpy(h_data, d_data, N3D_tf * sizeof(cufftComplex), cudaMemcpyDeviceToHost) == cudaSuccess );
  /*
  // transfer host results to parameter arrays
  for(i = 0; i < N3D; i++)
    {
      reel_arc_shift[i] = h_data[i].x;
      imag_arc_shift[i] = h_data[i].y;
    }
  */
  copy_3DmatrixC2C_in_smaller(h_data, reel_arc_shift, imag_arc_shift, N_tf, N_src);

  cufftDestroy(plan);
  // free(h_data);
  cudaFreeHost(h_data);
  cudaFree(d_data); 
}








// reel_arc_shift et imag_arc_shift réécrits 
void
TF3D_cuda_svg(int NX, int NY, int NZ, double *reel_arc_shift, double *imag_arc_shift)
{
  int i, N3D = NX * NY * NZ;
  
  cufftHandle plan;
  cufftComplex *h_data, *d_data;
  
  // allocate array in CPU space
  assert(h_data = (cufftComplex *) calloc(N3D, sizeof(cufftComplex)));
  
  // allocate array in GPU memory
  cudaMalloc((void**) &d_data, N3D * sizeof(cufftComplex));
  
  // set host array to values to process
  for(i = 0 ; i < N3D ; i++)
    {
      h_data[i].x = reel_arc_shift[i];
      h_data[i].y = imag_arc_shift[i];
    }
  // and clone it to its GPU twin
  cudaMemcpy(d_data, h_data, N3D * sizeof(cufftComplex), cudaMemcpyHostToDevice);


  // compute GPU FFT (Z2Z: doubleComplex to doubleComplex)
  assert( cufftPlan3d(&plan, NX, NY, NZ, CUFFT_C2C) == CUFFT_SUCCESS );
  assert( cufftExecC2C(plan, d_data, d_data, CUFFT_FORWARD) == CUFFT_SUCCESS );


  
  // retrieve GPU results
  cudaMemcpy(h_data, d_data, N3D * sizeof(cufftComplex), cudaMemcpyDeviceToHost); 

  // transfer host results to parameter arrays
  for(i = 0; i < N3D; i++)
    {
      reel_arc_shift[i] = h_data[i].x;
      imag_arc_shift[i] = h_data[i].y;
    }

  cufftDestroy(plan);
  free(h_data);
  cudaFree(d_data); 
}







void
cufftResult_display(cufftResult result)
{
  switch (result) {
  case CUFFT_SUCCESS:
    cout << endl << "CUFFT_SUCCESS";
    break;
  case CUFFT_EXEC_FAILED:
    cout << endl << "CUFFT_FAILURE";
    break;
  case CUFFT_INVALID_VALUE:
    cout << endl << "CUFFT_INVALID_VALUE";
    break;
  case CUFFT_INVALID_PLAN:
    cout << endl << "CUFFT_INVALID_PLAN";
    break;
  case CUFFT_SETUP_FAILED:
    cout << endl << "CUFFT_SETUP_FAILED";
    break;
  default:
    cout << endl << "unknown code";
  }
}




