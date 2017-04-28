#ifndef _FOCUS_CORRECTOR_
#define _FOCUS_CORRECTOR_



#ifndef OMP_THREADS
#define OMP_THREADS 8
#endif


#include <iostream>
using namespace std;

extern size_t g_fftw_threads;


#include "cv.h"
#include "highgui.h"

#include "vImg.hpp"

#include <AIR_Volume.h>
#include "FFTW_Image.hpp"



#include "cvDisplayArray.h"






// classe pour l'autofocus des hologrammes
// nécessite beaucoup d'effets de bord et d'allocations initiales
//
// idée: chaque hologramme (fréquentiel, avec Re et Im) peut être rétropropagé( d )
// trouver la valeur de d => hologramme soit le plus net et focus dans l'espace image


class Focus_Corrector{



  // =========================================================================== 
  // Data
  // =========================================================================== 

  
 private:


  // =========================================================================== 
  // parameters
  

  // spectre à refocaliser
  double *p_input_spectrum_r, *p_input_spectrum_i;
  // taille du spectre
  size_t p_dimX, p_dimY; 
  size_t p_nbPixels; // déduit

  size_t p_rayon;
  size_t p_nb_planes_retroprop; // nombre de rétropropagations effectuées pour chaque holograme

  // varie à chaque angle
  size_t p_spec_posX, p_spec_posY; // position du spéculaire 
  size_t p_current_angle;

  
  // pour produire des fichiers de debug (repropagations essentiellement)
  string debug_out_dir;
  size_t debug_nb_angles;
  bool i_debug;
  
  // montrer valeurs minimisées 
  bool i_show_vals;
  
  // montrer les images openCV 
  bool i_show_pics;

  

  // =========================================================================== 
  // local data

  
  // objet de traitement FFTW 2D
  FFTW_Image<double, fftw_complex> *i_Fourier2D;
  
  bool i_allocated, i_spectrum_set;
  
  
  // ----------------------------------------
  // images & volumes

  double *spectrePropag_normRe; double *spectrePropag_normIm;
  double *spectrePropag_normRe_shift; double *spectrePropag_normIm_shift;

  double *planObjetRe_shift; double *planObjetIm_shift;

  // volume stockant en mémoire l'ensemble des rétropropagations effectuées pour l'angle courant
  //double *objet3DRe; double *objet3DIm; // tué
  AIR_Volume<double> *V_Repropag_R, *V_Repropag_I;

  // volume constitué des nb_angle meilleurs hologrammes repropagés (debug)
  AIR_Volume<double> *V_RepropagBest_R, *V_RepropagBest_I;

  // objet refocalisé
  double *planObjetRe_refocal; double *planObjetIm_refocal;   
  // spectre correspondant
  double *spectreRe_refocal; double *spectreIm_refocal; 
  // idem avec franges remises = résultat final
  double *spectreRe_refocal_final; double *spectreIm_refocal_final; 

  // valeurs des gradients pour chaque rétropropagation
  double *stack_gradients;

  // images intermédiaires de rétropropagation
  double* l_retrop_kz; 
  double* l_retrop_rephase_r; double *l_retrop_rephase_i;
  double* l_retrop_propag_r; double *l_retrop_propag_i;

  // gradient: images intermédiaires pour calcul (la plupart ne sont plus utilisées)
  double* l_grad_ampli;
  double* l_grad_ampliGrad;
  double* l_grad_phase;
  double* l_grad_ampli1; double* l_grad_ampli2;
  double* l_grad_GradImage;

  // gradient: images opencv permanentes 
  Mat* l_cv_img, *l_cv_img_result;
  Mat* l_cv_grad_x, *l_cv_grad_y;
  
  // fenêtre dans gradient (opt)
  bool i_window_set;
  size_t i_window_ul_x, i_window_ul_y, i_window_lr_x, i_window_lr_y;   

  //

  // =========================================================================== 
  // Constructors & Destructors
  // =========================================================================== 
  

 public:


  Focus_Corrector() {
    
    p_dimX = p_dimY = p_nbPixels = 0;
    p_rayon = p_nb_planes_retroprop = p_spec_posX = p_spec_posY = 0;
    p_current_angle = 0;

    i_allocated = i_debug = i_spectrum_set = i_window_set = false;
    i_window_ul_x = i_window_ul_y = i_window_lr_x = i_window_lr_y = 0;
    debug_nb_angles = 1;
    V_RepropagBest_R = V_RepropagBest_I = NULL;

    i_show_pics = i_show_vals = false;
  }


  ~Focus_Corrector() { 
    if (i_allocated) 
      this -> unallocate(); 
  }
      

  NREADER(bool, i_allocated, get_allocated);


  // =========================================================================== 
  // Methods
  // =========================================================================== 


 private:


  // performs retropropagation of input spectrum (instance var) for given delta value
  // saves result in p_propag_spectrum
  void retropropagation(int delta, double* dst_real, double* dst_imag);

  
  // computes gradient of a complex image (not spectrum) of size (p_dimX, p_dimY)
  double gradient_complex( double* img_real, double* img_imag);

  double gradient_complex_window( double* img_real, double* img_imag, size_t ul_x, size_t ul_y);


  // =========================================================================== 
  

 public:

  
  // 1x au début
  // performs all allocations given the size of input spectrum
  void set_dimensions( size_t dimX, size_t dimY, size_t rayon )
  {
    p_rayon = rayon;
    p_dimX = dimX; p_dimY = dimY;
    p_nbPixels = p_dimX * p_dimY;
  }
  
  
  // 1x au début, non répétable
  // defines retropropagation to be done half_steps times before and after
  void set_retropropagation_steps(size_t half_steps) { 
    p_nb_planes_retroprop = 2 * half_steps; 
  }
  
  
  // 1x après les deux précédentes
  void allocate();
  void unallocate();


  // à chaque nouvel angle
  void set_input_spectrum_infos( size_t spec_posX, size_t spec_posY, size_t current_angle)
  { 
    p_spec_posX = spec_posX; p_spec_posY = spec_posY;
    p_current_angle = current_angle;
  }


  // à chaque nouvel angle
  // reads from outside the data hologram to refocus
  void get_input_spectrum( double* src_real, double* src_imag );


  // performs retropropagations, computes sum of gradients for each candidate, selects best
  void process();
  
  // writes to outside the refocused version of input hologram
  void set_output_spectrum( double* dst_real, double* dst_imag );


  // optionnal: enable the use of windowed functions
  void set_window(size_t ul_x, size_t ul_y, size_t lr_x, size_t lr_y)
  {
    i_window_ul_x = ul_x; i_window_ul_y = ul_y;
    i_window_lr_x = lr_x; i_window_lr_y = lr_y;
    i_window_set = true;
  }


  // --------------------------------------------------
  // montrer les valeurs de gradient de chaque repropagation 
  NREADWRITER(bool, i_show_vals, show_vals);

  // précédent + sauver les piles de repropagations pour les angles 1 à nb_angle
  void debug_mode(string output_dir, size_t nb_angles); 

  // affichage opencv 
  void set_show_pics()
  {
    i_show_pics = true;
    namedWindow("ydiff", CV_WINDOW_AUTOSIZE); // penser à détruire
  }


  // =========================================================================== 
};


#endif // _FOCUS_CORRECTOR_
