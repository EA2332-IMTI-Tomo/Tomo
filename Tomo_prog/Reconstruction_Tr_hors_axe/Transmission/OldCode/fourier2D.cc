#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>
#include <assert.h>
#include "macros.h"

/*
#ifdef CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#endif
*/

using namespace std;



//#include "main.h"
#include "fourier2D.h"




// =============================================================================
// gestion wisdom
// =============================================================================


// ----------------------------------------------------------------------
// calcul du nom du fichier wisdom


void
compute_wisdom_filename_f(char* result, size_t dim_x, size_t dim_y, bool forw_or_back, size_t nb_threads)
{
  assert(result);
  sprintf(result, "fftwf_%dx%d_%s_%dT", dim_x, dim_y, \
	  (forw_or_back) ? "forw" : "back", \
	  nb_threads);
}


void
compute_wisdom_filename_d(char* result, size_t dim_x, size_t dim_y, bool forw_or_back, size_t nb_threads)
{
  assert(result);
  sprintf(result, "fftw_%dx%d_%s_%dT", dim_x, dim_y, \
	  (forw_or_back) ? "forw" : "back", \
	  nb_threads);
}


// ----------------------------------------------------------------------
// import: renvoie faux si fichier non trouvé


bool 
importWisdom_f(const char* filename)
{
  FILE* wisdom_file;
  wisdom_file = fopen(filename, "r");

  if (wisdom_file) 
    {
      fftwf_import_wisdom_from_file(wisdom_file);
      fclose(wisdom_file);
      return true;
    }
  
  return false;
}


bool 
importWisdom_d(const char* filename)
{
  FILE* wisdom_file;
  wisdom_file = fopen(filename, "r");

  if (wisdom_file) 
    {
      fftw_import_wisdom_from_file(wisdom_file);
      fclose(wisdom_file);
      return true;
    }
  
  return false;
}


// ----------------------------------------------------------------------
// save


void 
saveWisdom_f(const char* filename)
{
  FILE* wisFile;
  assert(wisFile = fopen(filename, "w"));
  fftwf_export_wisdom_to_file(wisFile);
  fclose(wisFile);
}


void 
saveWisdom_d(const char* filename)
{
  FILE* wisFile;
  assert(wisFile = fopen(filename, "w"));
  fftw_export_wisdom_to_file(wisFile);
  fclose(wisFile);
}



// ----------------------------------------------------------------------
// calcul wisdom


void 
TF2D_compute_wisdom_f(int taille_x, int taille_y, bool forw_or_back, size_t nb_threads)
{
  
  int N = taille_x * taille_y;

  fftwf_complex* in = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * N);
  fftwf_complex* out = (fftwf_complex*) fftw_malloc(sizeof(fftwf_complex) * N);
  
  fftwf_plan_with_nthreads(nb_threads);
  fftwf_plan p;


  ARRAY_DEC_ALLOC(wisdom_file, 256, char);
  compute_wisdom_filename_f(wisdom_file, taille_x, taille_y, forw_or_back, nb_threads);
  
  cerr << endl << "computing wisdom, plz wait";
  p = fftwf_plan_dft_2d( taille_x, taille_y, in, out, forw_or_back ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_EXHAUSTIVE );
    
  assert(p);
  fftwf_execute(p); /* repeat as needed */
  
  saveWisdom_f(wisdom_file);
  cerr << endl << "saved: " << wisdom_file;


  fftwf_destroy_plan(p);
}


void 
TF2D_compute_wisdom_d(int taille_x, int taille_y, bool forw_or_back, size_t nb_threads)
{
  
  int N = taille_x * taille_y;

  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  
  fftw_plan_with_nthreads(nb_threads);
  fftw_plan p;


  ARRAY_DEC_ALLOC(wisdom_file, 256, char);
  compute_wisdom_filename_d(wisdom_file, taille_x, taille_y, forw_or_back, nb_threads);
  
  cerr << endl << "computing wisdom, plz wait";
  p = fftw_plan_dft_2d( taille_x, taille_y, in, out, forw_or_back ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_EXHAUSTIVE );
    
  assert(p);
  fftw_execute(p); /* repeat as needed */
  
  saveWisdom_d(wisdom_file);
  cerr << endl << "saved: " << wisdom_file;


  fftw_destroy_plan(p);
}



// =============================================================================
// calcul FFT 2D avec wisdom si disponible
// =============================================================================


fftwf_complex*
fftwf_volume_alloc(size_t dim_x, size_t dim_y)
{
  fftwf_complex *vol;
  vol = (fftwf_complex *) fftwf_malloc(sizeof (fftwf_complex) * dim_x * dim_y);
  return vol;
}


fftw_complex*
fftw_volume_alloc(size_t dim_x, size_t dim_y)
{
  fftw_complex *vol;
  vol = (fftw_complex *) fftw_malloc(sizeof (fftw_complex) * dim_x * dim_y);
  return vol;
}


// =============================================================================
// 
// temp_1 et temp_2 sont les volumes de type interne fftw_complex 
// (afin éviter l'allocation/désallocation)
// wisdom_filename: si fourni, utilisé afin de calculer plus vite. Sinon, mode ESTIMATE
// !! c'est de votre responsabilité de vérifier que le nombre de thread correspond aussi au wisdom


void 
TF2D_f(float* entree_reelle, float* entree_imag, float* fft_reel, float* fft_imag, size_t taille_x, size_t taille_y, fftwf_complex* temp_1, fftwf_complex* temp_2, const char* wisdom_filename, size_t nb_threads ) 
{
  static int calls = 0;
  int N = taille_x * taille_y;


  fftwf_plan_with_nthreads(nb_threads);
  fftwf_plan p;
  
  for(int cpt = 0 ; cpt < N ; cpt++)
    {
      temp_1[cpt][0] = entree_reelle[cpt];
      temp_1[cpt][1] = entree_imag[cpt];
    }


  bool have_wisdom = (wisdom_filename != NULL);
  if (have_wisdom)
    {
      importWisdom_f(wisdom_filename);
      p = fftwf_plan_dft_2d( taille_x, taille_y, temp_1, temp_2, FFTW_FORWARD, FFTW_EXHAUSTIVE | FFTW_WISDOM_ONLY );
    }
  else
    {
      p = fftwf_plan_dft_2d( taille_x, taille_y, temp_1, temp_2, FFTW_FORWARD, FFTW_ESTIMATE );
    }


  assert(p);
  fftwf_execute(p); /* repeat as needed */
  
  for(int cpt = 0; cpt < (N); cpt++)
    {
      fft_reel[cpt] = temp_2[cpt][0];
      fft_imag[cpt] = temp_2[cpt][1];
    }

  fftwf_destroy_plan(p);
  calls++;
  //   if (g_fftw_threads > 1)  void fftw_cleanup_threads(void);
}



void 
TF2D_d(double*  entree_reelle, double* entree_imag, double* fft_reel, double* fft_imag, size_t taille_x, size_t taille_y, fftw_complex* temp_1, fftw_complex* temp_2, const char* wisdom_filename, size_t nb_threads ) 
{
  static int calls = 0;
  int N = taille_x * taille_y;
  size_t cpt;

  fftw_plan_with_nthreads(nb_threads);
  fftw_plan p;

#pragma omp parallel for num_threads( 3 ) private(cpt)
  for(cpt = 0 ; cpt < N ; cpt++)
    {
      temp_1[cpt][0] = entree_reelle[cpt];
      temp_1[cpt][1] = entree_imag[cpt];
    }


  bool have_wisdom = (wisdom_filename != NULL);
  if (have_wisdom)
    {
      //if (! calls ) 
      importWisdom_d(wisdom_filename);
      p = fftw_plan_dft_2d( taille_x, taille_y, temp_1, temp_2, FFTW_FORWARD, FFTW_EXHAUSTIVE | FFTW_WISDOM_ONLY );
    }
  else
    {
      p = fftw_plan_dft_2d( taille_x, taille_y, temp_1, temp_2, FFTW_FORWARD, FFTW_ESTIMATE );
    }


  assert(p);
  fftw_execute(p); /* repeat as needed */
  
#pragma omp parallel for num_threads( 3 ) private(cpt)
  for(cpt = 0; cpt < (N); cpt++)
    {
      fft_reel[cpt] = temp_2[cpt][0];
      fft_imag[cpt] = temp_2[cpt][1];
    }

  fftw_destroy_plan(p);
  calls++;
  //   if (g_fftw_threads > 1)  void fftw_cleanup_threads(void);
}
 




















// ----------------------------------------------------------------------
// oldcode
// ----------------------------------------------------------------------

/*




void 
TF2D_eco(double entree_reelle[], double entree_imag[], double fft_reel[], double fft_imag[], int taille_x, int taille_y, fftw_complex* temp_1, fftw_complex* temp_2)
{
  static int calls = 0;
  int N = taille_x * taille_y;

  
  //fftw_plan_with_nthreads( 3 );

  fftw_plan p;
  
  for(int cpt = 0 ; cpt < N ; cpt++)
    {
      temp_1[cpt][0] = entree_reelle[cpt];
      temp_1[cpt][1] = entree_imag[cpt];
    }

  //calcul du plan, parametre servant a calculer et optimiser le FFT
  //  if (! calls) 
  importWisdom_d("./wisdom2D_double_512x512");

  p = fftw_plan_dft_2d( taille_x, taille_y, temp_1, temp_2, FFTW_FORWARD, FFTW_EXHAUSTIVE | FFTW_WISDOM_ONLY );

  assert(p);
  fftw_execute(p); // repeat as needed
  
  for(int cpt = 0; cpt < (N); cpt++)
    {
      fft_reel[cpt] = temp_2[cpt][0];
      fft_imag[cpt] = temp_2[cpt][1];
    }

  fftw_destroy_plan(p);
  calls++;
  //   if (g_fftw_threads > 1)  void fftw_cleanup_threads(void);
}


 */


// ----------------------------------------------------------------------
// 
// ----------------------------------------------------------------------

/*

// uniquement pour calculer le fichier wisdom
void 
TF2D_eco_dummy(double entree_reelle[], double entree_imag[], double fft_reel[], double fft_imag[], int taille_x, int taille_y, fftw_complex* temp_1, fftw_complex* temp_2)
{
int N = taille_x * taille_y;

//   if (g_fftw_threads > 1) fftw_plan_with_nthreads( g_fftw_threads );

fftw_plan p2D;
  
for(int cpt = 0 ; cpt < N ; cpt++)
{
temp_1[cpt][0] = entree_reelle[cpt];
temp_1[cpt][1] = entree_imag[cpt];
}

//calcul du plan, parametre servant a calculer et optimiser le FFT
p2D = fftw_plan_dft_2d( taille_x, taille_y, temp_1, temp_2, FFTW_FORWARD, FFTW_EXHAUSTIVE);
assert(p2D);

fftw_execute(p2D); // repeat as needed 
  
saveWisdom_d("./wisdom2D_double");
  
fftw_destroy_plan(p2D);
cout << endl << "wisdom plan computed";
exit(EXIT_SUCCESS);
  
//   if (g_fftw_threads > 1)  void fftw_cleanup_threads(void);
}
 
*/

 /*


   void 
   TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y)
   {
   int N=taille_x*taille_y;
   //Déclaration des variables pour la FFT : entree,sortie et "fftplan"
   fftw_complex *in, *out;
   fftw_plan p;
   //Réservation memoire
   in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);



   //Récupération de l'image dans la partie reelle de l'entree
   for(int cpt=0;cpt<N;cpt++)
   {
   in[cpt][0]=entree_reelle[cpt];
   in[cpt][1]=entree_imag[cpt];
   }


   //calcul du plan, parametre servant a calculer et optimiser le FFT
   p=fftw_plan_dft_2d( taille_x,  taille_y, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
   fftw_execute(p); // repeat as needed 

   for(int cpt=0;cpt<(N);cpt++)
   {

   fft_reel[cpt]=out[cpt][0];
   fft_imag[cpt]=out[cpt][1];
   }


   fftw_destroy_plan(p);
   fftw_free(in); fftw_free(out);
   }

 */

 // reel_arc_shift et imag_arc_shift réécrits 
 // FFTW_BACKWARD or FFTW_FORWARD

 /*

   void
   TF3D(size_t cube_edge, double *reel_arc_shift, double *imag_arc_shift, int fftw_constant_BACKorFORW)
   {
   int N3D = cube_edge * cube_edge * cube_edge;

   //Déclaration des variables pour la FFT : entre,sortie et "fftplan"
   fftw_complex *in3D, *out3D;

   if (g_fftw_threads > 1)
   fftw_plan_with_nthreads( g_fftw_threads );   // si parallele 

   fftw_plan p3D;
   //Réservation memoire
   in3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3D);
   //out3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3D);
   // p = fftw_plan_dft_2d(nx,ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

   //Récupération de l'image dans la partie reelle de l'entree
   for(int cpt=0;cpt<(N3D);cpt++)
   {
   in3D[cpt][0]=reel_arc_shift[cpt];
   in3D[cpt][1]=imag_arc_shift[cpt];
   }
   // destruction "virtuelle" des tableaux
   //   delete[] reel_arc_shift;
   //   delete[] imag_arc_shift;

   // reel_arc_shift et imag_arc_shift non nécessaires à ce stade
   out3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3D);
 
   //calcul du plan, parametre servant a calculer et optimiser le FFT
  

   p3D = fftw_plan_dft_3d(cube_edge, cube_edge, cube_edge, in3D, out3D, fftw_constant_BACKorFORW, FFTW_ESTIMATE);
   //                     int n0  n1     n2   fftw_complex*    int sign
   fftw_execute(p3D); // repeat as needed
    

   fftw_destroy_plan(p3D);
   fftw_free(in3D);
  
   for(int cpt=0;cpt<(N3D);cpt++)
   {
   reel_arc_shift[cpt]=out3D[cpt][0];
   imag_arc_shift[cpt]=out3D[cpt][1];
   }
   fftw_free(out3D);


   if (g_fftw_threads > 1)
   void fftw_cleanup_threads(void); 
   }



 */
