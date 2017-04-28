#ifndef __FOURIER_2D__
#define __FOURIER_2D__


#ifdef CUDA
#include "cuFourier.h"
#endif

#include <fftw3.h>

/* ---------------------------------------------------------------------- */

extern int g_fftw_threads;

/* ---------------------------------------------------------------------- */


void
compute_wisdom_filename_f(char* result, size_t dim_x, size_t dim_y, bool forw_or_back, size_t nb_threads);

void
compute_wisdom_filename_d(char* result, size_t dim_x, size_t dim_y, bool forw_or_back, size_t nb_threads);



bool 
importWisdom_f(const char* filename);

bool 
importWisdom_d(const char* filename);



void 
saveWisdom_f(const char* filename);

void 
saveWisdom_d(const char* filename);


void 
TF2D_compute_wisdom_f(int taille_x, int taille_y, bool forw_or_back, size_t nb_threads);

void 
TF2D_compute_wisdom_d(int taille_x, int taille_y, bool forw_or_back, size_t nb_threads);


fftwf_complex*
fftwf_volume_alloc(size_t dim_x, size_t dim_y);

fftw_complex*
fftw_volume_alloc(size_t dim_x, size_t dim_y);


void 
TF2D_f(float* entree_reelle, float* entree_imag, float* fft_reel, float* fft_imag, size_t taille_x, size_t taille_y, fftwf_complex* temp_1, fftwf_complex* temp_2, const char* wisdom_filename, size_t nb_threads );

void 
TF2D_d(double*  entree_reelle, double* entree_imag, double* fft_reel, double* fft_imag, size_t taille_x, size_t taille_y, fftw_complex* temp_1, fftw_complex* temp_2, const char* wisdom_filename, size_t nb_threads );
 












#endif // __FOURIER__



/* ---------------------------------------------------------------------- */
/*
void 
TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y);


// version économe de TF2D, requiert qu'on lui fournisse deux tableaux temp_1 et temp_2 alloués à taille_x * taill
e_y éléments
void 
TF2D_eco(double entree_reelle[], double entree_imag[], double fft_reel[], double fft_imag[], \
         int taille_x, int taille_y,                                    \
         fftw_complex* temp_1, fftw_complex* temp_2);

// version trafiquée pour créer du wisdom. sale
void 
TF2D_eco_dummy(double entree_reelle[], double entree_imag[], double fft_reel[], double fft_imag[], \
         int taille_x, int taille_y,                                    \
         fftw_complex* temp_1, fftw_complex* temp_2);

// présuppose que les volumes d'entrée aient été circshiftés au préalable
void
TF3D(size_t cube_edge, double *inout_reel, double *inout_imag, int fftw_constant_BACKorFORW);
//FFTW_BACKWARD or FFRW_FORWARD


// s'utilise directement avec des volumes non-shiftés
void
fourier3DInverse(double* reel_arc, double* imag_arc, int tab_size, int cube_edge, bool doCuda);
*/
