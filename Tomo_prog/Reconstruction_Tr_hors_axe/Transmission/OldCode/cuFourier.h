#ifndef __CU_FOURIER__
#define __CU_FOURIER__



#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>


/* ====================================================================== */
// 
/* ====================================================================== */


/* #include "../CUDA/cuda_by_example/common/book.h" */
// pris du bouquin


static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


static char *mkString[9]=
{
  (char *) "CUFFT_SUCCESS",
  (char *) "CUFFT_INVALID_PLAN",
  (char *) "CUFFT_ALLOC_FAILED",
  (char *) "CUFFT_INVALID_TYPE",
  (char *) "CUFFT_INVALID_VALUE",
  (char *) "CUFFT_INTERNAL_ERROR",
  (char *) "CUFFT_EXEC_FAILED",
  (char *) "CUFFT_SETUP_FAILED",
  (char *) "CUFFT_INVALID_SIZE"
};


static void HandleFFTError( cufftResult err,
                         const char *file,
                         int line ) {
    if (err != CUFFT_SUCCESS) {
      printf( "\n CuFFT Failure: %s in %s at line %d\n", mkString[(int)err],
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_FFT_ERROR( err ) (HandleFFTError( err, __FILE__, __LINE__ ))



/* ====================================================================== */
// 
/* ====================================================================== */


// transformée de fourier, Complex to Complex, single precision
// forward_or_backward: true = forward, false = backward
// nécessite que les volumes d'entrée aient été circshiftés, ne circshifte pas ceux de sortie
void
TF3D_cuda(int NX, int NY, int NZ, double *reel_arc_shift, double *imag_arc_shift, bool forward_or_backward);


void
TF3D_cuda_cubepower(int N_src, int N_tf, double *reel_arc_shift, double *imag_arc_shift, bool forward_or_backward);


void 
TF2D_cuda(double entree_reelle[], double entree_imag[], double fft_reel[], double fft_imag[], int taille_x, int taille_y);


void
cufftResult_display(cufftResult result);

#endif // __CU_FOURIER__
