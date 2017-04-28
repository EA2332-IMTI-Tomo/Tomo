#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>


// #include "vChrono.h"

#include "avPoint.h"

#include "recon_includes.h"
#include "util.h"
// #include "util_Magick.h"
#include "util_Image.h"
//#include "readloop.h"
#include "Compute.h"
// #include "memory.h"
#include "vectra.h"

// #include "fourier2D.h"
#include "vChronos.h"

#include <boost/thread.hpp>  
#include <boost/date_time.hpp>



#include "FFTW_Image.hpp"


#ifdef SSE4
#include "FFTW_Image_SSE4.h"
#endif

#include "Holo_Process.h"





// jeu de fonctions hui/matthieu pour Autofocus
// #include "Refocus.h"


void AutofocusGrad_jojo(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe, double *sortieIm, size_t dimX, size_t dimY, double rayon, size_t SpecPosX, size_t SpecPosY, int cpt_angle);

 
// =============================================================================
// =============================================================================


void
square_crop(double *src, double *dst, size_t src_size, size_t dst_size);


void 
cpu_kernel_mapping(int l_xm0, int l_ym0, size_t p_Nxmax, size_t p_Nymax, \
		   float p_rayon, float p_delta_zmax, double zm0, \
		   size_t c_dv0s2rf, size_t c_dv1s2rf, size_t c_dv0rf, size_t c_dv2s2rf, size_t c_dv1xdv0rf, \
		   unsigned short int *i_sup_redon, RECON_TYPE *i_reel_arc, RECON_TYPE *i_imag_arc, \
		   double *fft_reel_shift_norm, double *fft_imag_shift_norm);

 
// =============================================================================
// =============================================================================


Holo_Process::Holo_Process()
{
  i_prefs_file = i_prefs_size_image = i_prefs_output_vols = i_prefs_rayon = i_prefs_path = false;
  i_compute_ready = i_batch_ready = i_proper_names_p = false;
  show_hdp_p = false;
  
  d_save_stack_spots_p = false; d_VF_StackSpots_R = NULL;
  d_save_stack_windows_p = false; d_VF_StackWindows_R = NULL;
  d_save_stack_holograms_p = false; d_VF_StackHolograms_R = NULL;
  d_show_peaks_p = false;

  // on traite séparément la gestion du hors-axe, car elle est facultative
  // les autres valeurs sont obligatoires, et la classe vérifie que les 
  // set_ respectifs ont été apellés
  p_off_axis_mode = p_psh_mode = false; 
  p_off_axis_circle_cx = p_off_axis_circle_cy = p_off_axis_circle_r = 0;
  p_focus_correction_mode = p_modref_correction_mode = p_read_from_mem_mode = false;

  holo_filename = NULL;

  i_angle_last_done = 0;
  i_Refocalizer = NULL;
}

 
// =============================================================================
// =============================================================================


Holo_Process::~Holo_Process()
{
  set_free();
}

 
// =============================================================================
// =============================================================================


void
Holo_Process::set_psh()
{
  p_off_axis_mode = false;
  p_psh_mode = true;
}


void
Holo_Process::set_off_axis_single(size_t off_axis_circle_cx, size_t off_axis_circle_cy, size_t off_axis_r)
{
  p_off_axis_mode = true;
  p_psh_mode = false;
  p_off_axis_circle_cx = off_axis_circle_cx;
  p_off_axis_circle_cy = off_axis_circle_cy;

  // correction par défaut d'une erreur constatée en HA+HDP 
  // si le radius du cercle de découpe est omis, on laisse p_Nxmax à défaut (nb: en HDP seul, il est calculé)
  p_off_axis_circle_r = off_axis_r ? off_axis_r : p_Nxmax;

  
  // p_rayon; Le rayon est celui de la sphère sur laquelle on prélève les calottes.
  // p_Nxmax est le rayon de la calotte prélevée
  cout << endl << "OFF AXIS MODE: p_Nxmax & p_Nymax = " << p_Nxmax;
  p_Nxmax = p_Nymax = p_off_axis_circle_r; // rayon de la calotte 
  cout << endl << "OFF AXIS MODE: p_Nxmax & p_Nymax set to: " << p_Nxmax;
  
}


void
Holo_Process::set_off_axis_and_psh(size_t off_axis_circle_cx, size_t off_axis_circle_cy, size_t off_axis_r)
{
  p_off_axis_mode = true;
  p_psh_mode = true;
  p_off_axis_circle_cx = off_axis_circle_cx;
  p_off_axis_circle_cy = off_axis_circle_cy;
  p_off_axis_circle_r = off_axis_r ? off_axis_r : p_Nxmax;
}

 
// =============================================================================
// =============================================================================


void
Holo_Process::set_autofocus_gradient(size_t nb_planes_eachside)
{
  ASSERT( ! p_focus_correction_mode );

  const size_t c_Nxmax2 = 2 * p_Nxmax, c_Nymax2 = 2 * p_Nymax;

  MSG_ASSERT(i_prefs_size_image && i_prefs_rayon, "Unable to set autofocus before set_image_settings() and set_rayon_settings() are called");

  i_Refocalizer = new Focus_Corrector();
  i_Refocalizer -> set_dimensions(c_Nxmax2, c_Nymax2, p_rayon);
  i_Refocalizer -> set_retropropagation_steps( nb_planes_eachside );
  i_Refocalizer -> allocate();


  cerr << endl << "AUTOFOCUS MODE accounted" << endl;
  p_focus_correction_mode = true;
}


// =============================================================================
// =============================================================================


bool
Holo_Process::set_image_settings(size_t window_edge_x, size_t window_edge_y, size_t window_size, bool autocenter_window)
{
  // =============================================================================
  // A) on tâche de récupérer et vérifier la taille des images acquises

  MSG_ASSERT( i_prefs_file, \
	      "Holo_Process::set_image_settings: set_file_settings or set_memory_settings must be called first");
  MSG_ASSERT( ! p_read_from_mem_mode, \ 
	      "inutile en mode lecture mémoire");

  // si on lit de la mémoire, la taille des images est déjà connue. 
  // sinon, on déduit la taille des images à partir de la première
  
  bool file_found_p = false;
  size_t max_tries = 25;
  boost::posix_time::seconds sleepTime(3);

  while (! file_found_p)
    {
      // on teste un certain nombre du fichier de peur d'erreurs
      for (size_t i = 0; i < max_tries; i++)
	{
	  compute_hologram_filename(p_images_radix, holo_filename, i, 1);  
	  if (vectra::file_exists_p(holo_filename))
	    {file_found_p = true; break;}
	}

      // si aucune image n'a été trouvée, c'est aussi peut-être parce que l'acquisition n'a pas commencé   
      if (! file_found_p)
	{
	  cerr << endl << "holo_process / image_settings : waiting for acquisition to start";
	  boost::this_thread::sleep(sleepTime);     
	}
    }

  MSG_ASSERT( file_found_p, "Holo_Process::set_image_settings : ne trouve aucune image à analyser, extraction taille hologrammes  impossible" );

  i_vImgHoloReader = new vImg<unsigned char>(holo_filename);
  p_image_dim_x = i_vImgHoloReader -> get_x_dim();
  p_image_dim_y = i_vImgHoloReader -> get_y_dim();
  ASSERT(p_image_dim_x * p_image_dim_y);

  
  // =============================================================================
  // B) Gestion de le fanêtre de découpe si besoin
  // si la taille d'image est lue, la fenêtre de découpe est un paramètre fixe de la reconstruction

  p_window_edge_x = window_edge_x;
  p_window_edge_y = window_edge_y;
  p_window_size = window_size;

  if (autocenter_window)
    {
      p_window_edge_x = (p_image_dim_x - window_size) / 2;
      p_window_edge_y = (p_image_dim_y - window_size) / 2;
    }

  MSG_ASSERT( p_window_size && ( p_window_edge_x + p_window_size <= p_image_dim_x ) \
	      && ( p_window_edge_y + p_window_size <= p_image_dim_y ), \
	      "Holo_Process::set_image_settings : fenêtre de découpe sort de l'image");


  c_N = window_size * window_size;
  i_holo_resize_needed_p = (p_window_size != p_image_dim_x || p_window_size != p_image_dim_y);
  cout << endl;
  if (i_holo_resize_needed_p)
    cout <<"besoin de découper les acquisition";
  else 
    cout << "acquisitions à même taille que fenêtre"; cout.flush();
  
  
  // =============================================================================
  // C) Bilan 
  
  fprintf(stdout, "\ Taille des interférogrammes: ( %d x %d ) \n", p_image_dim_x, p_image_dim_y); 
  if (i_holo_resize_needed_p)
    fprintf(stdout, "\ Taille de la fenêtre de découpe: ( %d x %d ) \n", p_window_size, p_window_size);


  return i_prefs_size_image = true;
}

 
// =============================================================================
// =============================================================================


bool 
Holo_Process::set_output_volumes(AIR_Volume<RECON_TYPE>* V_fourier_R, AIR_Volume<RECON_TYPE>* V_fourier_I, AIR_Volume<unsigned short int>* V_supredon)
{
  this -> V_out_fourier_R = V_fourier_R;
  this -> V_out_fourier_I = V_fourier_I;
  this -> V_out_supredon = V_supredon;

  // check qui va bien
  
  ASSERT(V_out_fourier_R -> get_allocated() && V_out_fourier_I -> get_allocated() \
	 && V_out_supredon -> get_allocated()); 
  ASSERT(V_out_fourier_R -> same_dim_p(*V_out_fourier_I));
  ASSERT(V_out_fourier_R -> same_dim_p(*V_out_supredon));
  ASSERT(V_out_fourier_R -> same_order_p(*V_out_fourier_I));
  ASSERT(V_out_fourier_R -> same_order_p(*V_out_supredon));
  
  i_reel_arc = V_out_fourier_R -> get_data();
  i_imag_arc = V_out_fourier_I -> get_data();
  i_sup_redon = V_out_supredon -> get_data();
  ASSERT( i_reel_arc != NULL && i_imag_arc != NULL && i_sup_redon != NULL ); 


  return i_prefs_output_vols = true;
}


// =============================================================================
// =============================================================================


bool
Holo_Process::set_rayon_settings(size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, size_t xm0_limite, size_t ym0_limite, float rayon, float delta_zmax)
{
  p_Nxmax = Nxmax;
  p_Nymax = Nymax;
  p_Nxmax_Rf = Nxmax_Rf;
  p_rayon = rayon;
  p_delta_zmax = delta_zmax;

  ASSERT( p_Nxmax * p_Nymax * p_Nxmax_Rf * p_rayon * p_delta_zmax);

  c_rayon_inf = xm0_limite * xm0_limite + ym0_limite * ym0_limite;
  c_dv0rf = 4 * Nxmax_Rf; 
  c_dv0s2rf = 2 * Nxmax_Rf;
  c_dv1s2rf = c_dv0s2rf;
  c_dv1xdv0rf = 16 * Nxmax_Rf * Nxmax_Rf;
  c_dv2s2rf = c_dv0s2rf;
  
  return i_prefs_rayon = true;
}

 
// =============================================================================
// =============================================================================


bool
Holo_Process::set_file_settings(const char* images_radix, bool proper_names)
{
  p_images_radix = images_radix;
  i_proper_names_p = proper_names;
  ASSERT(holo_filename = (char *) calloc(strlen(p_images_radix) + 256 + GENERATED_FILENAME_LEN, sizeof(char)));

  p_read_from_mem_mode = false;
  return i_prefs_file = true;
}
 
// =============================================================================
// =============================================================================

/*
bool
Holo_Process::set_memory_settings(unsigned char** images_array, size_t image_dim_x, size_t image_dim_y)
{
  p_image_dim_x = image_dim_x;
  p_image_dim_y = image_dim_y;
  ASSERT(p_image_dim_x * p_image_dim_y);

  i_images_array = images_array;

  p_read_from_mem_mode = true;
  return i_prefs_file = true;
}
*/
  

bool
Holo_Process::set_memory_settings(const AIR_Volume<unsigned char> *V_stack )
{
  p_image_dim_x = V_stack -> get_x_dim();
  p_image_dim_y = V_stack -> get_y_dim();
  ASSERT(p_image_dim_x * p_image_dim_y);
  
  i_V_hologramStack = ( AIR_Volume<unsigned char> *) V_stack;
  
  if (holo_filename == NULL) 
    ASSERT(holo_filename = (char *) calloc(strlen(p_images_radix) + 256 + GENERATED_FILENAME_LEN, sizeof(char)));

  p_window_edge_x = p_window_edge_y = 0;
  p_window_size = p_image_dim_x;
  ASSERT(p_window_size == p_image_dim_y); 
  c_N = p_window_size * p_window_size;
  i_holo_resize_needed_p = false;

  i_vImgHoloReader = new vImg<unsigned char>(p_image_dim_x, p_image_dim_y);

  p_read_from_mem_mode = true;
  i_prefs_size_image = true;
  return i_prefs_file = true;
}

// =============================================================================
// =============================================================================


bool
Holo_Process::set_modref_correction(const char* mod_ref_file)
{
  p_modref_correction_mode = vectra::file_exists_p(mod_ref_file);
  if (p_modref_correction_mode) 
    {
      p_modref_filename = mod_ref_file;    
      cout << endl << "mod ref correction mode with file:" << mod_ref_file << endl;
    }
  return p_modref_correction_mode;
}


// =============================================================================
// =============================================================================


bool
Holo_Process::set_path_settings(std::map<string, string> &map_paths)
{
  l_MASKPATH = map_paths["MASKPATH"].c_str();
  return i_prefs_path = true;
}


// =============================================================================
// =============================================================================


bool 
Holo_Process::set_prepare()
{
  ASSERT( i_prefs_file && i_prefs_rayon && i_prefs_size_image &&i_prefs_path  );
  //  i_prefs_output_vols

  // --------------------------------------------------
  // allocations groupées
  //tampon_image_cv = cvCreateImage( cvSize(p_window_size, p_window_size), IPL_DEPTH_8U, 1); 


  ARRAY_SSE_ALLOC(masque_hamming, c_N, 32, unsigned char);  
  ARRAY_SSE_ALLOC(phase1, c_N, 32, unsigned char); 
  ARRAY_SSE_ALLOC(phase2, c_N, 32, unsigned char); 
  ARRAY_SSE_ALLOC(phase3, c_N, 32, unsigned char); 
  ARRAY_SSE_ALLOC(phase4, c_N, 32, unsigned char);  

  //vImg *vIphase1
  vImgPhase1 = new vImg<unsigned char>(phase1, p_window_size, p_window_size);
  vImgPhase2 = new vImg<unsigned char>(phase2, p_window_size, p_window_size);
  vImgPhase3 = new vImg<unsigned char>(phase3, p_window_size, p_window_size);
  vImgPhase4 = new vImg<unsigned char>(phase4, p_window_size, p_window_size);


  // expérimental, à désallouer
  ARRAY_SSE_ALLOC(phase1n, c_N, 32, double); 
  ARRAY_SSE_ALLOC(phase2n, c_N, 32, double); 
  ARRAY_SSE_ALLOC(phase3n, c_N, 32, double); 
  ARRAY_SSE_ALLOC(phase4n, c_N, 32, double);  

  
  vImg<unsigned char>* vImgModRef;
  if (p_modref_correction_mode)
    { 
      ARRAY_SSE_ALLOC(mod_ref, c_N, 32, unsigned char); 
      ARRAY_SSE_ALLOC(mod_refn, c_N, 32, double); 
      vImgModRef = new vImg<unsigned char>(mod_ref, p_window_size, p_window_size);
    }

  ARRAY_SSE_ALLOC(plan_reel, c_N, 32, double);  
  ARRAY_SSE_ALLOC(plan_imag, c_N, 32, double);

  // nb de pixels de l'image fourier 2D (2 * Nxmax, 2 * Nymax)
  const size_t l_NxNy4 = 4 * p_Nxmax * p_Nymax;
  ARRAY_ALLOC(i_centre, l_NxNy4, int);
  ARRAY_SSE_ALLOC(fft_reel_shift, l_NxNy4, 32, double); 
  ARRAY_SSE_ALLOC(fft_imag_shift, l_NxNy4, 32 ,double); 
  ARRAY_SSE_ALLOC(fft_module_shift, l_NxNy4, 32, double);  
  ARRAY_SSE_ALLOC(fft_reel_shift_norm, l_NxNy4, 32, double);  
  ARRAY_SSE_ALLOC(fft_imag_shift_norm, l_NxNy4, 32, double); 
  if (p_focus_correction_mode) {
    ARRAY_SSE_ALLOC(fft_reel_shift_norm_focus, l_NxNy4, 32, double);  
    ARRAY_SSE_ALLOC(fft_imag_shift_norm_focus, l_NxNy4, 32, double); 
  }
  

  d_VF_StackSpots_R = NULL; // on ne peut l'allouer que dans le batch

  // TF d'images caméra découpées 
  //
  // wisdom, double, t. cumulé. 1T -> 2.70, 2T-> 1.70 3T -> 1.30s. 4T -> 1.30s! 
  #ifdef SSE4
  VF_Fourier2D = new FFTW_Image_SSE4<double, fftw_complex>(p_window_size, p_window_size, 3);
  #else
  VF_Fourier2D = new FFTW_Image<double, fftw_complex>(p_window_size, p_window_size, 3); //g_fftw_threads);
  #endif
 
  VF_Fourier2D -> allocate();
  VF_Fourier2D -> set_wisdom_dir("/usr/local/phd/fftw_wisdom");
  //   VF_Fourier2D -> create_wisdom(true); VF_Fourier2D -> create_wisdom(false);


  // TF dans le champ de fourier
  #ifdef SSE4
  I_Fourier_Holo = new FFTW_Image_SSE4<double, fftw_complex>(2 * p_Nxmax, 2 * p_Nymax, 3);
  #else
  I_Fourier_Holo = new FFTW_Image<double, fftw_complex>(2 * p_Nxmax, 2 * p_Nymax, 3);
  #endif
  I_Fourier_Holo -> allocate();
  I_Fourier_Holo -> set_wisdom_dir("/usr/local/phd/fftw_wisdom");
  

  // --------------------------------------------------


  i_nb_proj = i_points_faux = i_jumeau_elimine = i_centres_exclus = 0;
  //i_cpt_max = 0; i_max_val = 0;
  

  // --------------------------------------------------
  // gestion masques

  // peut être généré aussi, à voir
  c_mask_size = 30;
  char* l_chaine = str_concat(l_MASKPATH, "/k_30x30.pgm");
  ARRAY_SSE_ALLOC(cache_jumeau, c_mask_size * c_mask_size, 32, unsigned char);
  remplit_tableau_cv(cache_jumeau, l_chaine, c_mask_size, c_mask_size);


  i_hamming = new vImg<unsigned char>(masque_hamming, p_window_size, p_window_size);
  i_hamming -> tuckey(0.15f);
    
  
  // --------------------------------------------------
  // gestion mod_ref (correction du fond par référence)

  if (p_modref_correction_mode)
    {
      i_vImgHoloReader -> reload_pgm(p_modref_filename);
      vImgModRef -> copy_from_bigger(*i_vImgHoloReader, p_window_edge_x, p_window_edge_y); // copié dans mod_ref
      
      // élimination des zéros en prévision de division
#pragma omp parallel for num_threads( 4 ) 
      for(size_t pixel = 0; pixel < c_N; pixel++)
	{
	  if (mod_ref[pixel] == 0) {mod_ref[pixel] = 1;}
	}
#pragma omp barrier


      normalize_vectra1(mod_ref, mod_refn, c_N);
    }

  // --------------------------------------------------


  cout << endl << "ready for launch";
  return i_compute_ready = true;
}


// =============================================================================
// =============================================================================


void
Holo_Process::set_free()
{
  if (i_compute_ready)
    {
      //cvReleaseImage( &tampon_image_cv );

      free(masque_hamming);
      _mm_free(phase1); _mm_free(phase2); _mm_free(phase3); _mm_free(phase4);
      free(phase1n); free(phase2n); free(phase3n); free(phase4n);
      free(plan_reel); free(plan_imag); 

      if (p_modref_correction_mode) 
	{ free(mod_ref); free(mod_refn); }


      free(fft_reel_shift); free(fft_imag_shift); free(fft_module_shift); 
      free(fft_reel_shift_norm); free(fft_imag_shift_norm);

      // if (p_focus_correction_mode)
      //      {free(fft_reel_shift_norm); free(fft_imag_shift_norm);}
      // pas possible pour cause affectation pointeurs
      
      free(cache_jumeau);
      free(i_centre);

      VF_Fourier2D -> unallocate();
      I_Fourier_Holo -> unallocate();

      //releaser les vImg///

      /* plante :/
      if (p_focus_correction_mode && i_Refocalizer -> get_allocated())
	i_Refocalizer -> unallocate();
      */
      
    }


  // on vide les infos de débug
  if (d_keep_centers)
    {
      ASSERT(v_angles.size() == v_centers.size());
      size_t x, y;
      double dist;
      const avPoint<int> Centre(p_Nxmax, p_Nymax, 0);
      avPoint<int> Point(p_Nxmax, p_Nymax, 0);
      
      size_t angle_count = v_angles.size();
      

      std::ofstream file_centers(d_centers_log_file.c_str());
      if (file_centers.is_open())
	{
	  
	  for (size_t i = 0; i < angle_count; i++)
	    { 
	      y = v_centers[i] / (2 * p_Nxmax);
	      x = v_centers[i] % (2 * p_Nxmax);
	      Point.set(x, y, 0);
	      dist = Centre.distance_eucl(Point);
	  
	      file_centers << "angle: " << v_angles[i] << " peak: [" << x << " " << y << "] dist: " << dist << endl;
	    }
	  
	  file_centers.flush();
	  file_centers.close();
	}
      else cerr << endl << "can't write to file" << d_centers_log_file;

      // si je fais cela, ça détruit a posteriori le fichier (sans déconner putain)
      //v_angles.clear(); v_centers.clear(); 
    }
  
  i_compute_ready = false;
}


// =============================================================================
// =============================================================================


bool
Holo_Process::batch_prepare(size_t start_angle, size_t batch_size, size_t max_error)
{
  ASSERT(i_compute_ready);

  v_batch_angles.reserve(batch_size);
  
  size_t error_count = 0;
  size_t i = 0;
    

  if (! p_read_from_mem_mode)
    {
      // CAS 1: lecture des interférogrammes par fichiers
      // la synchronization est effectuée en vérifiant la présence sur disque des fichiers

      while (v_batch_angles.size() < batch_size)
	{
	  size_t angle_num = start_angle + i;
	  i++;
	  compute_hologram_filename(p_images_radix, holo_filename, angle_num, 1);
      
	  // on contrôle la validité d'un angle en vérifiant que la dernière 
	  // image de cet angle existe
	  if (vectra::file_exists_p(holo_filename))
	    {
	      v_batch_angles.push_back(angle_num); 
	      error_count = 0; 
	    }
	  else
	    { 
	      error_count++; 

	      if (error_count < 2)
		cerr << endl << "wrong: num=" << angle_num << " " << holo_filename << " batch_size:" << batch_size << " vector size:" <<v_batch_angles.size(); 
	      else if (error_count == 2)
		cerr << endl << "skipping errors for this batch";
	      
	    }
	        
	  if (error_count >= max_error) 
	    break;
	}

    }
  else
    {
      // CAS 2: lecture des interférogrammes en mémoire
      // Là, tous les angles déclarés existent (alors que sur disque, on en efface)
      // devrait au moins marcher avec numéro du dernier angle acquis partagé
      // Pour le moment, ça ne fonctionne QUE si tous les interféros sont déjà présents

      while (v_batch_angles.size() < batch_size)
	v_batch_angles.push_back( ++i ); 


      // test mémoire si dispo
      unsigned char dummy;
      size_t interf_count = batch_size;
      if (p_psh_mode) interf_count *= 4;

      for (size_t i = 0; i < batch_size; i++)
	{
	  unsigned char* myslice = (*i_V_hologramStack)[i];
	  dummy = myslice[0];
	}

    }

  //cout << endl << endl << "------------" << v_batch_angles.size() << " " << batch_size;


  if (v_batch_angles.size())
    i_batch_ready = true;
  else
    i_batch_ready = false; 
  
  return i_batch_ready;
}


// =============================================================================
// =============================================================================


bool
Holo_Process::batch_prepare_realtime(size_t max_batch_size, size_t min_batch_size_start, size_t mysleep_ms)
{
  ASSERT(i_compute_ready);
  //vChronos timer_proc("realtime batch making"); timer_proc.start();

  boost::posix_time::seconds sleepTime(mysleep_ms);

  //const size_t max_batch_size = 200;
  //const size_t min_batch_size_start = 15;
  const size_t max_read_errors = 2;

  v_batch_angles.reserve(max_batch_size);
  
  size_t error_count = 0;
  size_t i = 0;
  size_t last_inserted = 0 + i_angle_last_done; // le premier angle s'apelle 1
  size_t target_angle;
  bool had_nap = false;

  MSG_ASSERT(! p_read_from_mem_mode, "pas encore programmé ce jour");

  while (v_batch_angles.size() < max_batch_size)
    {

      // on n'arrive *vraiment* pas à lire l'angle suivant (max_read_errors essais)
      if (error_count >= max_read_errors) 
	{
	  // si on a un minimum d'angles OK, on lance le batch en l'état
	  if (v_batch_angles.size() > min_batch_size_start) 
	    { 
	      i_batch_ready = (v_batch_angles.size() > 0);
	      //timer_proc.stop();
	      return i_batch_ready;
	    }
	  else // pas assez d'angles
	    {
	      // première solution: on dort 500ms 
	      if (! had_nap)
		{
		  cerr << endl << "holo_process:: sleeping for more angles";
		  cout << endl << "expects: " << holo_filename; 
		  cout << "HDP: " << p_psh_mode << "HA: "  << p_off_axis_mode;
		  boost::this_thread::sleep(sleepTime);     
		  had_nap = true;
		}
	      else // si on a déjà dormi 500ms
		{ 
		  // si on est au début de l'acquis, on est tenu de dormir indéfiniment jusqu'au lancement de l'instrument
		  if (i_angle_last_done < 1)
		    {
		      cerr << endl << "holo_process:: waiting for acquisition to start";
		      // pas besoin de faire du zèle tant que l'acquis ne démarre pas
		      boost::this_thread::sleep(sleepTime);     
		      boost::this_thread::sleep(sleepTime);     
		      boost::this_thread::sleep(sleepTime);     
		      boost::this_thread::sleep(sleepTime);     
		      had_nap = false;
		    }
		  else
		    {
		      // si on est plus loin que le début, on peut supposer qu'on est à la fin 
		      // c'est pas grave si détection abusive de fin, car aucun hologramme n'est effacé du disque
		      cerr << endl << "slept too much, no more angles to read from microscope";
		      break;
		    }
		}
	    }
	}
      
      target_angle = last_inserted + 1;
      // en HDP, on attend de voir le 4è interférogramme avant de déclarer l'angle valide.
      size_t holo_num = ( p_psh_mode ? 4 : 1 ); 
      compute_hologram_filename(p_images_radix, holo_filename, target_angle, holo_num);
      
      if (vectra::file_exists_p(holo_filename))
	{
	  v_batch_angles.push_back(target_angle); 
	  error_count = 0; 
	  last_inserted++;
	}
      else
	error_count++;	
    }

  
  if (v_batch_angles.size())
    i_batch_ready = true;
  else
    i_batch_ready = false;
  
  //timer_proc.stop();
  return i_batch_ready;
}


// =============================================================================
// =============================================================================


void
Holo_Process::batch_launch()
{
  ASSERT(i_batch_ready && i_prefs_output_vols);

  MSG_ASSERT( (! p_off_axis_mode || p_off_axis_circle_r), "hors axe configuré mais pas de rayon de découpe" );


  // constantes locales

  size_t batch_size = v_batch_angles.size();

  const size_t c_window_dim_x = p_window_size;
  const size_t c_window_dim_y = p_window_size;
  const size_t c_Nxmax2 = 2 * p_Nxmax, c_Nymax2 = 2 * p_Nymax;
  const size_t c_NxNy4 = 4 * p_Nxmax * p_Nymax;
  const size_t c_mask_halfsize = c_mask_size / 2;
  const double c_inv_255 = 1.0f / 255.0f;
  const float c_rayon_sqr = sqr(p_rayon);
  const float c_delta_zmax_sqr = sqr(p_delta_zmax);
    
  

  // pour svg pile images debug
  double* spot_data = NULL;
  double* window_data = NULL;
  double* holo_data = NULL;
  unsigned int spot_offset = 0;
  unsigned int window_offset = 0;
  unsigned int holo_offset = 0;
  //FFTW_Image<double, fftw_complex> *I_Fourier_Holo; 


  // --------------------------------------------------
  if (d_save_stack_spots_p)
    {
      ASSERT( !d_VF_StackSpots_R );
      d_VF_StackSpots_R = new AIR_Volume<double> ( 2 * p_Nxmax, 2 * p_Nymax, batch_size );
      d_VF_StackSpots_R -> allocate();
  
      spot_data = (d_VF_StackSpots_R -> get_allocated()) ? d_VF_StackSpots_R -> get_data() : NULL;
    }  
  if (d_save_stack_windows_p)
    {
      ASSERT( !d_VF_StackWindows_R );
      d_VF_StackWindows_R = new AIR_Volume<double> ( 2 * p_Nxmax, 2 * p_Nymax, batch_size );
      d_VF_StackWindows_R -> allocate();
  
      window_data = (d_VF_StackWindows_R -> get_allocated()) ? d_VF_StackWindows_R -> get_data() : NULL;
    }
  if (d_save_stack_holograms_p)
    {
      ASSERT( !d_VF_StackHolograms_R );
      d_VF_StackHolograms_R = new AIR_Volume<double> ( 2 * p_Nxmax, 2 * p_Nymax, batch_size );
      d_VF_StackHolograms_R-> allocate();
  
      holo_data = (d_VF_StackHolograms_R -> get_allocated()) ? d_VF_StackHolograms_R -> get_data() : NULL;
    }
  // --------------------------------------------------

  //   cout << endl << "BATCH READY, SIZE IS:" << batch_size;
  //   cout << endl << "STACK IS:" << 2 * p_Nxmax << "x" <<  2 * p_Nymax << "x" << batch_size;


  // ---------------------------------------------------------------------------
  // iteration sur les angles
  // ---------------------------------------------------------------------------

  /*
  // section de code pour afficher uniquement les hologrammes recalés
  double* holorecal_tmp;
  ARRAY_ALLOC(holorecal_tmp, 2 * p_Nxmax * 2 * p_Nymax, double);
  cvDisplayArray<double>* showRecal_tmp;
  showRecal_tmp = new cvDisplayArray<double>("HolorecalTMP", holorecal_tmp, 2 * p_Nxmax, 2 * p_Nymax);
  */


  for (size_t i = 0; i < batch_size; i++) 
    {
      size_t current_angle = v_batch_angles[i];

      // affichage périodique progression
      if (! (current_angle % 50))
	cout << endl << "processing: " << current_angle;

      /*
      this -> hologram_show_debug(current_angle, holorecal_tmp);
      showRecal_tmp -> updateImage();
      showRecal_tmp -> showImage();
      cvWaitKey(30);
      continue;
      */

      // lit les interférogrammes, TF, crop dans fréquences instrumentales
      // résultat écrit dans fft_reel_shift & fft_imag_shift
      this -> hologram_compute( current_angle );

      
      // ----------------------------------------------------------------------
      // Recherche du maximum du module dans ref non centré + 
      // Normalisation par le pic central
      // ----------------------------------------------------------------------


      const size_t peak_1D_pos = peak_normalize(c_NxNy4, fft_reel_shift, fft_imag_shift, \
						fft_reel_shift_norm, fft_imag_shift_norm);

      // Coordonnées dans l'espace 2D à partir de l'indice 1D: (xc,yc)=(0,0)=en haut à gauche de l'image
      // mi: coordonnées objet (spéculaire?). l_xmi avant. c_Nymax2 = 2 * p_Nxmax
      const size_t peak_x = peak_1D_pos % (2 * p_Nxmax);
      const size_t peak_y = peak_1D_pos / (2 * p_Nymax);


      if (d_show_peaks_p) {
	cout << endl << "peak n°" << current_angle << "=" << peak_1D_pos << " x=" <<  peak_x << " y=" << peak_y << " val=(" << fft_reel_shift[peak_1D_pos] << ";" << fft_imag_shift[peak_1D_pos] << ")";      
	//cout << "HACK OLIVE:"; fft_reel_shift_norm[peak_1D_pos] = 0; fft_imag_shift_norm[peak_1D_pos] = 0;     
      }


      // renvoie la position du centre, qui peut servir
      // i_centre[cpt_max]++;  ne sert plus
      if (d_keep_centers) // _log_file) 
	{ 
	  v_angles.push_back( current_angle ); 
	  v_centers.push_back( peak_1D_pos ); 
	}
      
      
      // ----------------------------------------------------------------------
      // debug: sauvegarde de l'hologramme recalé dans un volume
      // ----------------------------------------------------------------------


      if(holo_data)
	{
	  // le circshift centré sur la position du pic...
	  size_t circshift_decal_x = 2 * p_Nxmax - peak_x;
	  size_t circshift_decal_y = 2 * p_Nymax - peak_y;
	  I_Fourier_Holo -> import_excentric(fft_reel_shift_norm, fft_imag_shift_norm, circshift_decal_x, circshift_decal_y);
	  // va nous donner un hologramme spatial sans franges d'inclinaison
	  I_Fourier_Holo -> set_fourier_backward();
	  // après re-circshift centré
	  I_Fourier_Holo -> export_R_to(holo_data + holo_offset, true);

	  holo_offset += 2 * p_Nxmax * 2 * p_Nymax; 
	}


      // ----------------------------------------------------------------------
      // autofocus ici
      // ----------------------------------------------------------------------


      if (p_focus_correction_mode)
	{	  
	  /*
	    AutofocusGrad_jojo( fft_reel_shift_norm, fft_imag_shift_norm, \
	    fft_reel_shift_norm_focus, fft_imag_shift_norm_focus, \
	    c_Nxmax2, c_Nymax2, p_rayon, l_xmi, l_ymi, current_angle );
	  */
	 
	  i_Refocalizer -> set_input_spectrum_infos( peak_x, peak_y, current_angle );
	  i_Refocalizer -> get_input_spectrum( fft_reel_shift_norm, fft_imag_shift_norm );
	  i_Refocalizer -> process(); /// ICI
	  i_Refocalizer -> set_output_spectrum( fft_reel_shift_norm_focus, fft_imag_shift_norm_focus );
	}
      else // le programme est conçu en assumant qu'on refocuse par défaut, par commodité
	{
	  fft_reel_shift_norm_focus = fft_reel_shift_norm;
	  fft_imag_shift_norm_focus = fft_imag_shift_norm;
	}

      

      // ----------------------------------------------------------------------
      // Elimination du jumeau (ssi en HDP *pure* )
      // ----------------------------------------------------------------------


      if ( ! p_off_axis_mode )
	{
	  twin_eliminate_redux(peak_x, peak_y, fft_reel_shift_norm_focus, fft_imag_shift_norm_focus, cache_jumeau);
	  //twin_eliminate(l_xmi, l_ymi, fft_reel_shift_norm, fft_imag_shift_norm);
	} 

      

      // ----------------------------------------------------------------------
      // projection de l'holo (-> sphere cap) dans le volume 3D
      // ----------------------------------------------------------------------
      
      angle_mapping(current_angle, peak_x, peak_y, fft_reel_shift_norm_focus, fft_imag_shift_norm_focus);
      // fft_shift_norm: hologramme

      // ----------------------------------------------------------------------
      //if (show_hdp_p) cvWaitKey(1);
      

      
      // debug: sauvegarde du spot de fourier dans une pile
      if (spot_data)
	{ 
	  memcpy(spot_data + spot_offset, fft_reel_shift_norm_focus, 2 * p_Nxmax * 2 * p_Nymax * sizeof(double));
	  spot_offset += 2 * p_Nxmax * 2 * p_Nymax; 
	}
      
      // debug: sauvegarde de la découpe dans fourier dans une pile
      if (window_data)
	{ 
	  memcpy(window_data + window_offset, fft_reel_shift, 2 * p_Nxmax * 2 * p_Nymax * sizeof(double));
	  window_offset += 2 * p_Nxmax * 2 * p_Nymax; 
	}


      i_angle_last_done = current_angle;
    }
  // ---------------------------------------------------------------------------
  // FIN iteration sur les angles
  // ---------------------------------------------------------------------------

  

  // sauvegarde données contrôle (stack des spots de fourier 2D)
  if ( d_save_stack_spots_p && d_VF_StackSpots_R -> get_allocated() ) 
    {
      cout << endl << "fourier debug:" << d_stack_spots_file;
      d_VF_StackSpots_R -> change_files(d_stack_spots_file); //("svg_fourier_stack_spots");
      d_VF_StackSpots_R -> write_files();
      d_VF_StackSpots_R -> unallocate();
    }
  if ( d_save_stack_windows_p && d_VF_StackWindows_R -> get_allocated() ) 
    {
      cout << endl << "fourier debug:" << d_stack_windows_file;
      d_VF_StackWindows_R -> change_files(d_stack_windows_file); //("svg_fourier_stack_spots");
      d_VF_StackWindows_R -> write_files();
      d_VF_StackWindows_R -> unallocate();
    }
  if ( d_save_stack_holograms_p && d_VF_StackHolograms_R -> get_allocated() ) 
    {
      d_VF_StackHolograms_R -> change_files(d_stack_holograms_file); //("svg_fourier_stack_spots");
      d_VF_StackHolograms_R -> write_files();
      d_VF_StackHolograms_R -> unallocate();
    }

  v_batch_angles.clear();
  ASSERT( ! v_batch_angles.size() );

  //timer2.stop();

  i_batch_ready = false;
}


// =============================================================================
// =============================================================================


void
Holo_Process::display()
{
  printf("\ncampics: [%d x %d] | window size: %d  at (%d, %d)", p_image_dim_x, p_image_dim_y, p_window_size, p_window_edge_x, p_window_edge_y);
  printf("\nrayon_inf (norme xmo, ym0): %d", c_rayon_inf);
  printf("\nMax: Nx: %d, Ny: %d, NxRf: %d", p_Nxmax, p_Nymax, p_Nxmax_Rf);
  printf("\nr: %f , dmax: %f", p_rayon, p_delta_zmax);
} 


// =============================================================================
// =============================================================================


void
Holo_Process::d_record_centers(size_t angle_max, string file)
{
  d_centers_log_file = file;
  d_keep_centers = true;

  // à vider
  v_angles.reserve(angle_max);
  v_centers.reserve(angle_max);
}



// =============================================================================
// =============================================================================

void
Holo_Process::d_record_fourier_spots(string filename)
{
  d_save_stack_spots_p = true;
  d_stack_spots_file = filename;
  // rien de plus n'est possible sans être dans le batch
}


void
Holo_Process::d_record_fourier_windows(string file)
{
  d_save_stack_windows_p = true;
  d_stack_windows_file = file;
}


void
Holo_Process::d_record_final_holograms(string file)
{
  d_save_stack_holograms_p = true;
  d_stack_holograms_file = file;
}


// =============================================================================
// =============================================================================


void
Holo_Process::show_hdp()
{
  /*
  if (show_hdp_p) {
    showHDPr -> updateImage(); showHDPr -> showImage();
    showHDPi -> updateImage(); showHDPi -> showImage();
  }
  */
}


// =============================================================================
// =============================================================================


void
Holo_Process::hologram_show_debug(size_t current_angle, double* holo_recale_output)
{
  ASSERT(i_compute_ready);

  hologram_compute(current_angle);

  // ----------------------------------------------------------------------
  // Recherche du maximum du module dans ref non centré + 
  // Normalisation par le pic central
  // ----------------------------------------------------------------------
  
  const size_t c_NxNy4 = 4 * p_Nxmax * p_Nymax;
  const size_t peak_1D_pos = peak_normalize(c_NxNy4, fft_reel_shift, fft_imag_shift, \
					    fft_reel_shift_norm, fft_imag_shift_norm);
  
  // Coordonnées dans l'espace 2D à partir de l'indice 1D: (xc,yc)=(0,0)=en haut à gauche de l'image
  // mi: coordonnées objet (spéculaire?). l_xmi avant. c_Nymax2 = 2 * p_Nxmax
  const size_t peak_x = peak_1D_pos % (2 * p_Nxmax);
  const size_t peak_y = peak_1D_pos / (2 * p_Nymax);


  // le circshift centré sur la position du pic...
  size_t circshift_decal_x = 2 * p_Nxmax - peak_x;
  size_t circshift_decal_y = 2 * p_Nymax - peak_y;
  I_Fourier_Holo -> import_excentric(fft_reel_shift_norm, fft_imag_shift_norm, circshift_decal_x, circshift_decal_y);
  // va nous donner un hologramme spatial sans franges d'inclinaison
  I_Fourier_Holo -> set_fourier_backward();
  // après re-circshift centré
  I_Fourier_Holo -> export_R_to(holo_recale_output, true);
}


// =============================================================================
// =============================================================================


#include "Holo_Process_tools.cc"
