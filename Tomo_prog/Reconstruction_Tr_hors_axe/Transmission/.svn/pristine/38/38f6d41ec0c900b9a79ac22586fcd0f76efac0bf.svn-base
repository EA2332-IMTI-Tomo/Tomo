#ifndef __FOCUS_CORRECTOR_CC__
#define __FOCUS_CORRECTOR_CC__


#include "Focus_Corrector.h"

#include "FFTW_Image.hpp"
#include "vArray.h"
#include "util_Image.h"


typedef struct {
  double Re,Im;
}nb_complexe;


typedef struct {
  int x,y;
}Var2D;



// =============================================================================
// =============================================================================


void
Focus_Corrector::allocate()
{
  ASSERT( ! i_allocated );
 
  // a priori inutile
  ARRAY_ALLOC(p_input_spectrum_r, p_nbPixels, double);
  ARRAY_ALLOC(p_input_spectrum_i, p_nbPixels, double);
  // pas de free pour le moment

  ARRAY_ALLOC(spectrePropag_normRe, p_nbPixels, double);
  ARRAY_ALLOC(spectrePropag_normIm, p_nbPixels, double);
  ARRAY_ALLOC(spectrePropag_normRe_shift, p_nbPixels, double);
  ARRAY_ALLOC(spectrePropag_normIm_shift, p_nbPixels, double);
  ARRAY_ALLOC(planObjetRe_shift, p_nbPixels, double);
  ARRAY_ALLOC(planObjetIm_shift, p_nbPixels, double);
  

  // volume stockant en mémoire l'ensemble des rétropropagations effectuées pour l'angle courant
  //ARRAY_ALLOC(objet3DRe, p_nbPixels * p_nb_planes_retroprop, double);
  //ARRAY_ALLOC(objet3DIm, p_nbPixels * p_nb_planes_retroprop, double);
  
  V_Repropag_R = new AIR_Volume<double> ( p_dimX, p_dimY, p_nb_planes_retroprop );
  V_Repropag_I = new AIR_Volume<double> ( p_dimX, p_dimY, p_nb_planes_retroprop );

  V_Repropag_R -> allocate(); 
  V_Repropag_I -> allocate();

  //V_RepropagBest_R/I: dans debug uniquement

  ARRAY_ALLOC(planObjetRe_refocal, p_nbPixels, double);
  ARRAY_ALLOC(planObjetIm_refocal, p_nbPixels, double);
  ARRAY_ALLOC(spectreRe_refocal, p_nbPixels, double);
  ARRAY_ALLOC(spectreIm_refocal, p_nbPixels, double);
  ARRAY_ALLOC(spectreRe_refocal_final, p_nbPixels, double);
  ARRAY_ALLOC(spectreIm_refocal_final, p_nbPixels, double);

  // valeurs des gradients pour chaque rétropropagation
  ARRAY_ALLOC(stack_gradients, p_nb_planes_retroprop, double);


  // specifique retropropagation
  ARRAY_ALLOC(l_retrop_kz, p_nbPixels, double);
  ARRAY_ALLOC(l_retrop_rephase_r, p_nbPixels, double); 
  ARRAY_ALLOC(l_retrop_rephase_i, p_nbPixels, double);
  ARRAY_ALLOC(l_retrop_propag_r, p_nbPixels, double); 
  ARRAY_ALLOC(l_retrop_propag_i, p_nbPixels, double); 

  //spécifique calcul gradient
  ARRAY_ALLOC(l_grad_ampli, p_nbPixels, double);
  ARRAY_ALLOC(l_grad_ampliGrad, p_nbPixels, double);
  ARRAY_ALLOC(l_grad_phase, p_nbPixels, double);
  ARRAY_ALLOC(l_grad_ampli1, p_nbPixels, double);
  ARRAY_ALLOC(l_grad_ampli2, p_nbPixels, double);
  ARRAY_ALLOC(l_grad_GradImage, p_nbPixels, double);


  i_Fourier2D = new FFTW_Image<double, fftw_complex>(p_dimX, p_dimY, 3); //g_fftw_threads);
  i_Fourier2D -> allocate();
  i_Fourier2D -> set_wisdom_dir("/usr/local/phd/fftw_wisdom");


  // to free, someday!
  l_cv_img = new cv::Mat(p_dimX, p_dimY, CV_64F); 
  l_cv_grad_x = new cv::Mat(p_dimX, p_dimY, CV_64F); 
  l_cv_grad_y = new cv::Mat(p_dimX, p_dimY, CV_64F); 
  l_cv_img_result = new cv::Mat(p_dimX, p_dimY, CV_64F); 
  

  // alloué dans la fonction
  ///cv_grad_window = NULL;  // A  VOIR
  

  i_allocated = true;
}


// -----------------------------------------------------------------------------


void 
Focus_Corrector::unallocate()
{ 
  ASSERT( i_allocated );
  // attention, il s'avère que ça plante aux premiers free
  
  free(spectrePropag_normRe); free(spectrePropag_normIm);
  free(spectrePropag_normRe_shift); free(spectrePropag_normIm_shift);
  free(planObjetRe_shift); free(planObjetIm_shift);
  
  free(planObjetRe_refocal); free(planObjetIm_refocal);
  free(spectreRe_refocal); free(spectreIm_refocal);
  free(spectreRe_refocal_final); free(spectreIm_refocal_final);

  free(stack_gradients);

  free(l_retrop_kz);
  free(l_retrop_rephase_r); free(l_retrop_rephase_i);
  free(l_retrop_propag_r); free(l_retrop_propag_i);

  free(l_grad_ampli); free(l_grad_ampliGrad);
  free(l_grad_phase); 
  free(l_grad_ampli1); free(l_grad_ampli2);
  free(l_grad_GradImage);

  i_Fourier2D -> unallocate(); 

  
 
  V_Repropag_R -> unallocate(); 
  V_Repropag_I -> unallocate(); 

  if (i_debug) 
    {
      ASSERT( V_RepropagBest_R -> get_allocated() && V_RepropagBest_I -> get_allocated() );
      char stack_file[256]; 
      
      sprintf(stack_file, "./Repropagated_Best_R");
      V_RepropagBest_R -> change_files(stack_file);
      V_RepropagBest_R -> write_files();

      sprintf(stack_file, "./Repropagated_Best_I");
      V_RepropagBest_I -> change_files(stack_file);
      V_RepropagBest_I -> write_files();
      
      V_RepropagBest_R -> unallocate(); 
      V_RepropagBest_I -> unallocate(); 
    }


  i_allocated = false;
}


// -----------------------------------------------------------------------------


void
Focus_Corrector::debug_mode(string output_dir, size_t nb_angles) 
{
  MSG_ASSERT( i_allocated, "trop tôt");
  MSG_ASSERT( ! i_debug, "déjà fait");

  debug_out_dir = output_dir;
  debug_nb_angles = nb_angles;
  i_debug = i_show_vals = true; // show_pics trop intrusif

  V_RepropagBest_R = new AIR_Volume<double> ( p_dimX, p_dimY, debug_nb_angles );
  V_RepropagBest_I = new AIR_Volume<double> ( p_dimX, p_dimY, debug_nb_angles );

  V_RepropagBest_R -> allocate(); 
  V_RepropagBest_I -> allocate();
}


// =============================================================================
// =============================================================================


void 
Focus_Corrector::process()
{ 
  ASSERT( i_allocated );
  

  int xm0 = p_spec_posX - ( p_dimX / 2 );
  int ym0 = p_spec_posY - ( p_dimY / 2 );

  
  int decalage_spec_X = -xm0 + p_dimX / 2;
  int decalage_spec_Y = -ym0 + p_dimX / 2;  
  int decalage_specInv_X = xm0 - p_dimX / 2;
  int decalage_specInv_Y = ym0 - p_dimX / 2;

  
  //  cout << endl << "Focus corrector: peak "


  // pile des repropagations d'un même hologramme effectuées pour l'angle courant
  double *stack_R = V_Repropag_R -> get_data(); // objet3DRe
  double *stack_I = V_Repropag_I -> get_data(); // objet3DIm


  // volume constitué des nb_angle meilleurs hologrammes repropagés
  double *best_R, *best_I;
  if (i_debug) {
    best_R = V_RepropagBest_R -> get_data();
    best_I = V_RepropagBest_I -> get_data();
  }
  

  // ----------------------------------------
  // step 1: on calcule les rétropropagations du spectre courant pour p_nb_planes_retroprop valeurs de deltaZ
  // on les stocke dans un volume (R & I)
  // on stocke une valeur de gradient pour chaque image

  int deltaZ;
  for (int steps = 0; steps < p_nb_planes_retroprop; steps++)
    {
      deltaZ = steps - (p_nb_planes_retroprop / 2); 
      
      // prend p_input_spectrum_r/i et calcule la rétropropag courante (deltaZ) dans spectrePropag_norm
      this -> retropropagation(deltaZ, spectrePropag_normRe, spectrePropag_normIm);
      
      // Eliminer les franges d'inclinaision
      circshift2(spectrePropag_normRe, spectrePropag_normRe_shift, p_dimX, p_dimY, decalage_spec_X, decalage_spec_Y);
      circshift2(spectrePropag_normIm, spectrePropag_normIm_shift, p_dimX, p_dimY, decalage_spec_X, decalage_spec_Y);

      // notre spectre a son "0" à la position du spéculaire: on effectue un décalage non-centré
      
      i_Fourier2D -> import_from(spectrePropag_normRe_shift, spectrePropag_normIm_shift, false);
      i_Fourier2D -> set_fourier_backward();
      // circshift classique: on veut l'image résultante
      i_Fourier2D -> export_to(planObjetRe_shift, planObjetIm_shift, true);
      
      
      // Accumulation de l'image générée dans un volume de nb_retropag images
      memcpy( stack_R + steps * p_nbPixels, planObjetRe_shift, p_nbPixels * sizeof(double) );
      memcpy( stack_I + steps * p_nbPixels, planObjetIm_shift, p_nbPixels * sizeof(double) );
      
      // calcul du gradient "complexe" des deux images courantes et accumulation dans un tableau
      stack_gradients[steps] = this -> gradient_complex(planObjetRe_shift, planObjetIm_shift);
    }

  
  // sauvegarde de la pile de repropagations pour débogage
  if (i_debug && (p_current_angle < 100)) 
    {
      char stack_file[256]; 
      ASSERT( strlen( debug_out_dir.c_str() ) < 220 );
      sprintf(stack_file, "%s/stack_repro_a_%03d_I", debug_out_dir.c_str(), (int)p_current_angle);
      V_Repropag_I -> change_files(stack_file);
      V_Repropag_I -> write_files();
    }
  
  
  // on confie le tableau à vArray pour trouver le minimum
  vArray<double> Stack_gradients(stack_gradients, p_nb_planes_retroprop);
  size_t pos_minGrad = Stack_gradients.position_min();


  // affichages valeurs
  if (i_show_vals)
      cerr << endl << ">>> angle:" << p_current_angle << "/ min is: " << pos_minGrad << "     || " << Stack_gradients;
  if (pos_minGrad != (p_nb_planes_retroprop / 2)) // si delta != 0: ca a fonctionné
    cout << endl << "correction effective pour hologramme: " << p_current_angle << "( " \
	 << pos_minGrad << " : " << p_nb_planes_retroprop << " )";
  

  // on récupère le plan objet choisi
  memcpy(planObjetRe_refocal, stack_R + pos_minGrad * p_nbPixels, p_nbPixels * sizeof(double));
  memcpy(planObjetIm_refocal, stack_I + pos_minGrad * p_nbPixels, p_nbPixels * sizeof(double));
 
  
  // DEBUG: sauvegarde de l'image sélectionnée dans un volume de nb_angles couches
  if (i_debug)
    { // bug pas compris, revoir les historiques
      memcpy(best_R + (p_current_angle - 1) * p_nbPixels, planObjetRe_refocal, p_nbPixels * sizeof(double)); 
      memcpy(best_I + (p_current_angle - 1) * p_nbPixels, planObjetIm_refocal, p_nbPixels * sizeof(double)); 
    }



  ////VECTRA: l'objet est pas circshifté en sortie
  //// conclusion: il faut le circshifter là, et donc dans l'espace image, avant le repassage freq

  i_Fourier2D -> import_from(planObjetRe_refocal, planObjetIm_refocal, true);
  i_Fourier2D -> set_fourier_forward();
  i_Fourier2D -> set_fourier_normalize();
  i_Fourier2D -> export_to(spectreRe_refocal, spectreIm_refocal, false);

  // uncentered: donne que les franges!! seulement 2 marche
  ///Ajouter les franges d'inclinaision
  circshift2(spectreRe_refocal, spectreRe_refocal_final, p_dimX, p_dimY, decalage_specInv_X, decalage_specInv_Y); 
  circshift2(spectreIm_refocal, spectreIm_refocal_final, p_dimX, p_dimY, decalage_specInv_X, decalage_specInv_Y); 
}


// =============================================================================
// =============================================================================

// retroPropagSA(deltaZ, fft_shift_normRe, fft_shift_normIm, spectrePropag_normRe, spectrePropag_normIm, dimInit, rayon, CoordSpecI );

// void  retroPropagSA(float deltaZ, double *fft_shift_normRe, double *fft_shift_normIm, double * spectrePropag_normRe, double * spectrePropag_normIm, Var2D dimSpctHolo, double rayon, Var2D CoordSpecI)


void
Focus_Corrector::retropropagation( int delta, double* dst_real, double* dst_imag )
{
  ASSERT( i_allocated );
  MSG_ASSERT( i_spectrum_set, "no spectrum to propagate");

  float deltaZ= (float) delta / 20.0f;

  Var2D KMAX = { round(p_dimX / 2), round(p_dimY / 2)}; 
  Var2D CoordSpecC = {p_spec_posX -KMAX.x, p_spec_posY - KMAX.y};
  
  double r2 = p_rayon * p_rayon;
  double arg_kz_arc = 0, kz_arc =0;


  // kz, rephase, spectrePropag
  memset(l_retrop_kz, 0, p_nbPixels * sizeof(double));
  //memset(l_retrop_rephase_r, 0, p_nbPixels * sizeof(double));
  //memset(l_retrop_rephase_i, 0, p_nbPixels * sizeof(double));
  memset(l_retrop_propag_r, 0, p_nbPixels * sizeof(double));
  memset(l_retrop_propag_i, 0, p_nbPixels * sizeof(double));


  double kz0 = sqrt(r2 - CoordSpecC.x * CoordSpecC.x - CoordSpecC.y * CoordSpecC.y);
  int NXMAX_CARRE = KMAX.x * KMAX.x; 

  ////VECTRA: Ce calcul semble constant, le faire une fois pour toutes
  ///CALCULER LES KZ-----------------------------------
  for (int ky = -KMAX.y; ky < KMAX.y; ky++) 
    { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
      int ky_carre=ky*ky;

      for (int kx = -KMAX.x; kx < KMAX.x; kx++) 
	{ //on balaye l'image 2D en y, centre au milieu
	  //kyi*dim+kxi ...calcul du cpt du tableau 1D de l'image 2D
	  int cpt = (ky + KMAX.y ) * p_dimX + kx +KMAX.x; // dimSpctHolo.x
	  if(kx * kx + ky_carre < NXMAX_CARRE)//ne pas depasser l'ouverture numérique pour 1 hologramme
	    {
	      double kz_carre = r2 -kx * kx - ky_carre; //altitude au carré des données
	      l_retrop_kz[cpt] = sqrt(kz_carre) - kz0;
	    }
	  else  
	    { l_retrop_kz[cpt] = 0; }
	} 
    }
  //SAV(kz, NbPixHolo, "/home/hui/maniptomo/IDP/SA/kz.bin", FLOAT,"a+b");


  ///---------------------CALCULER LE SPECTRE REPROPAGE-----------------------------------
  for(int pix = 0; pix < p_nbPixels ; pix++)  // NbPixHolo
    {
      double kz_delta = l_retrop_kz[pix] * deltaZ;
      double rephase_r = cos( kz_delta );
      double rephase_i = sin( kz_delta );


      l_retrop_propag_r[pix] = p_input_spectrum_r[pix] * rephase_r - \
	p_input_spectrum_i[pix] * rephase_i;

      l_retrop_propag_i[pix] = p_input_spectrum_r[pix] * rephase_i + \
	p_input_spectrum_i[pix] * rephase_r;
    }


  int cptSpec = p_dimX * p_spec_posY + p_spec_posX;
  
  double max_module = pow( l_retrop_propag_r[cptSpec], 2) + pow( l_retrop_propag_i[cptSpec], 2);
  double max_part_reel = l_retrop_propag_r[cptSpec];
  double max_part_imag = l_retrop_propag_i[cptSpec];

  for(int cpt = 0; cpt < p_nbPixels; cpt++) 
    {
      dst_real[cpt] = (l_retrop_propag_r[cpt] * max_part_reel + l_retrop_propag_i[cpt] * max_part_imag) / max_module;
      dst_imag[cpt] = (l_retrop_propag_i[cpt] * max_part_reel - l_retrop_propag_r[cpt] * max_part_imag) / max_module;
    }

  // SAV_Re(spectrePropag_norm, NbPixHolo, "/home/mat/tomo_test/spectrePropag_Norm.bin", FLOAT,"a+b");
}


// =============================================================================
// =============================================================================


void
Focus_Corrector::set_output_spectrum( double* img_real, double* img_imag)
{
  ASSERT( i_allocated );
  memcpy(img_real, spectreRe_refocal_final, p_nbPixels * sizeof(double));
  memcpy(img_imag, spectreIm_refocal_final, p_nbPixels * sizeof(double));
}


// =============================================================================
// =============================================================================


void 
Focus_Corrector::get_input_spectrum( double* src_real, double* src_imag )
{
  ASSERT( i_allocated );
  //p_input_spectrum_r = src_real; p_input_spectrum_i = src_imag;
  memcpy( p_input_spectrum_r, src_real, p_nbPixels * sizeof(double) );
  memcpy( p_input_spectrum_i, src_imag, p_nbPixels * sizeof(double) );
  i_spectrum_set = true;
}

// =============================================================================
// =============================================================================


inline double
norm(double x, double y)
{
  return (sqrt ( x * x + y * y ) );
} 


void 
translateImg(Mat &img_src, Mat &img_dst, int offsetx, int offsety)
{
  Mat trans_mat = (Mat_<double>(2,3) << 1, 0, offsetx, 0, 1, offsety);
  warpAffine(img_src,img_dst,trans_mat,img_src.size());
}


// calcule la valeur de gradient de l'image d'entrée (supposée minimale au focus)
double
Focus_Corrector::gradient_complex( double* img_real, double* img_imag)
{

  /// calcul amplitude de l'image d'entrée -> l_cv_img (on se sert d'images préallouées)
  double* cv_data = (double *) l_cv_img -> data; // cast du pointeur de données 
  for(int pixel = 0; pixel < p_nbPixels; pixel++) 
    cv_data[pixel] = norm( img_real[pixel], img_imag[pixel] );

  // on décale l'image norme d'un pixel vertical (puisque la fibre est horizontale)
  Mat* cv_shifty = l_cv_grad_x;
  translateImg( *l_cv_img, *cv_shifty, 0, 1);

  // on soustrait à l'image non décalée => gradient y
  Mat* cv_diff = l_cv_img_result;
  cv::subtract( *l_cv_img, *cv_shifty, *cv_diff );
  
  
  // --------------------------------------------------
  // sélection d'une fenêtre si précisée
  size_t x_min = 0, x_max = p_dimX;
  size_t y_min = 0, y_max = p_dimY;
  if (i_window_set)
    {
      x_min = i_window_ul_x; x_max = i_window_lr_x;
      y_min = i_window_ul_y; y_max = i_window_lr_y;
    }

  // --------------------------------------------------
  // parcours des valeurs 
  size_t x, y;
  double sommegrad = 0;
  double* diff_data = (double *)cv_diff -> data;


  for (x = x_min; x < x_max; x++)
    for (y = y_min; y < y_max; y++)
      {
	size_t i = y * p_dimY + x;
	sommegrad += abs(diff_data[i]);
      }

  // --------------------------------------------------
  // displays
  if (i_show_vals)
    cout << sommegrad << " / ";

  
  if (i_show_pics) {
    imshow("ydiff", *cv_diff);
    cvWaitKey(0);
  }

  
  return sommegrad;
}


// *****************************************************************************
// 
// *****************************************************************************

double
Focus_Corrector::gradient_complex_window( double* img_real, double* img_imag, size_t ul_x, size_t ul_y)
{
  ASSERT(false);

  // marche sans doute mais pas encore appellé
  /*

  MSG_ASSERT(l_window_tmp1 -> get_allocated_p(), "call set_window first");

  // on découpe une fenêtre dans img_real / imag positionnée par ul (upper left corner)
  l_window_tmp1 -> import_from_bigger(img_real, p_dimX, p_dimY, ul_x, ul_y);
  l_window_tmp2 -> import_from_bigger(img_imag, p_dimX, p_dimY, ul_x, ul_y);

  // on calcule l'image module (amplitude)
  l_window_tmp3 -> module(*l_window_tmp1, *l_window_tmp2);
  double* grad_ampli = l_window_tmp3 -> get_data();


  // load double data into an openCV image
  double* cv_data = (double *) l_cv_img -> data;
  memcpy( cv_data, grad_ampli, p_nbPixels * sizeof(double) );


  /// sobel
  int scale = 1;
  int delta = 0;
  int ddepth = CV_64F;
  //Mat grad_x, grad_y, imgResult;

  /// Gradient X
  //Scharr( img1[0], grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
  cv::Sobel( *l_cv_img, *cv_grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
  /// Gradient Y
  //Scharr( img1[0], grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
  cv::Sobel( *l_cv_img, *cv_grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
  /// Total Gradient (approximate)
  cv::addWeighted( *cv_grad_x, 0.0, *cv_grad_y, 1.0, 0, *l_cv_img_result );


  double* cv_result_data = (double *) l_cv_img_result -> data;
  memcpy( l_grad_ampliGrad, cv_result_data, p_nbPixels * sizeof(double) );

  circshift2(l_grad_ampliGrad, l_grad_ampli1, p_dimX, p_dimY, 0, -1);
  circshift2(l_grad_ampliGrad, l_grad_ampli2, p_dimX, p_dimY, -1, 0);


  double sommegrad = 0;
  for(int cpt = 0; cpt < p_nbPixels; cpt++)
    {
      double ampliGrad_i = l_grad_ampliGrad[cpt];
      double gradImage_i =  norm( ampliGrad_i - l_grad_ampli1[cpt],	\
				  ampliGrad_i - l_grad_ampli2[cpt]);

      l_grad_GradImage[cpt] = gradImage_i; // pour débogage
      sommegrad += gradImage_i;
    }

    return sommegrad;
  */

  return 0;
}



// =============================================================================
// fonctions locales annexes
// =============================================================================



#endif //__FOCUS_CORRECTOR_CC__
