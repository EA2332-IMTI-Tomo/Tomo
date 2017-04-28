#ifndef _HOLO_PROCESS_
#define _HOLO_PROCESS_


#include <iostream>
using namespace std;

#include "cv.h"
#include "highgui.h"


#ifndef OMP_THREADS
#define OMP_THREADS 8
#endif


#include <AIR_Volume.h>
#include "FFTW_Image.hpp"


#include "Focus_Corrector.h"

#include "vImg.hpp"
#include "cvDisplayArray.h"

#include <iostream>
#include <vector>


extern size_t g_fftw_threads;


#ifndef RECON_TYPE
#define RECON_TYPE float
#endif



// classe pour le traitement progressif d'hologrammes
// nécessite beaucoup d'effets de bord et d'allocations initiales

class Holo_Process{


  // ===========================================================================
  // ===========================================================================
 protected:


  // ===========================================================================
  // input_data
  // ===========================================================================

  // volumes de sortie, alloués dans l'appellant 
  //
  // contiennent les coefficients de fourier non normalisés calculés ici
  AIR_Volume<float>* V_out_fourier_R;
  AIR_Volume<float>* V_out_fourier_I;
  // de même que l'histogramme des recouvrements (normalisation)
  AIR_Volume<unsigned short int> *V_out_supredon;

  // dimensions internes
  size_t p_Nxmax, p_Nymax, p_Nxmax_Rf;
  float p_rayon, p_delta_zmax;

  // taille images hologramme (caméra)
  size_t p_image_dim_x, p_image_dim_y;
  // taille et position fenêtre de découpe dans susnommées (fenêtre carrée)
  size_t p_window_edge_x, p_window_edge_y, p_window_size; 
  // si la fenêtre est plus petite que l'image, on doit cropper
  bool i_holo_resize_needed_p;
  // dans ce cas, on utilisera cet objet pour la lecture
  vImg<unsigned char>* i_vImgHoloReader; 

  const char* p_images_radix;
  const char* p_modref_filename;

  // gestion du hors-axe (cumulable au phase shifting )
  // (va juste modifier le crop Nxmax et désactiver la chasse au jumeau)
  bool p_psh_mode;
  bool p_off_axis_mode;
  size_t p_off_axis_circle_cx, p_off_axis_circle_cy, p_off_axis_circle_r;
  
  // gestion de la correction de focus
  bool p_focus_correction_mode;
  Focus_Corrector *i_Refocalizer;
  
  // doit-on exécuter correction par mod_ref (si demandé ET image modref lisible)
  bool p_modref_correction_mode;

  // Lecture des interférogrammes à partir de la mémoire (char**) ou du disque?
  bool p_read_from_mem_mode;


  // --------------------------------------------------
  // garde en mémoire le n° du dernier angle traîté par batch_launch
  size_t i_angle_last_done;


  
  // ===========================================================================
  // local_data
  // ===========================================================================

  // chemins
  const char* l_MASKPATH;


  // variables internes classe: a-t-on apellé tous les initialiseurs?
  bool i_prefs_file, i_prefs_size_image, i_prefs_output_vols, i_prefs_rayon, i_prefs_path;
  bool i_compute_ready, i_batch_ready;

 public: // temporaire
  // vecteur contenant les numéros d'angles valides concernés par le prochain calcul de batch
  vector<size_t> v_batch_angles;
 protected:

  // ===========================================================================
  
  RECON_TYPE *i_reel_arc, *i_imag_arc;
  unsigned short int* i_sup_redon;

  //IplImage *tampon_image_cv;
  unsigned char *masque_hamming, *cache_jumeau;
  // reservation de 4 tableaux (image 2D)  pour le phase shifting
  unsigned char *phase1, *phase2, *phase3, *phase4, *mod_ref;
  // de même, leurs versions normalisées (expérimental)
  double *phase1n, *phase2n, *phase3n, *phase4n, *mod_refn;

  vImg<unsigned char>* i_hamming; 

  // objets vImg utilisants phase1/4
  vImg<unsigned char> *vImgPhase1, *vImgPhase2, *vImgPhase3, *vImgPhase4;

  // si on lit à partir de la mémoire:
  unsigned char **i_images_array; // deprecated
  AIR_Volume<unsigned char> *i_V_hologramStack;

  // pour mettre la position des centres translatés, on crée une variable 2D de la taille d'un plan apres tomo
  int *i_centre;

  // classe de gestion des transformées de fourier 2D d'images lues de caméra
  FFTW_Image<double, fftw_complex> *VF_Fourier2D;

  // idem pour images de champ de fourier
  FFTW_Image<double, fftw_complex> *I_Fourier_Holo; 


  double *plan_reel, *plan_imag;
  // fft 2d: à remplacer par classe
  double *fft_reel_shift_norm, *fft_imag_shift_norm;
  double *fft_reel_shift, *fft_imag_shift, *fft_module_shift;
  double *fft_reel_shift_norm_focus, *fft_imag_shift_norm_focus;
  // fft_reel_tmp: avant Crop à NxMax   |  fft_reel: après

  // nom du fichier image (hologramme) en cours
  char *holo_filename;
  // noms au propre (%03d x2) ou pas (par défaut)
  bool i_proper_names_p; 


  
  // constantes après initialisation
  size_t c_rayon_inf;
  //nombre de pixel dans l'image de destination
  size_t c_N;
  //dimx_arc | dimension totale apres correction par Rf
  size_t c_dv0rf; 
  //dimx_arc / 2
  size_t c_dv0s2rf;
  //
  size_t c_dv1s2rf;
  //dimx_arc*dimy_arc
  size_t c_dv1xdv0rf;
  //dimz_arc
  size_t c_dv2s2rf;

  // taille totale du masque en pixel pour elimination objet jumeau;
  size_t c_mask_size;
  

  // variables de calcul
  size_t i_nb_proj, i_points_faux, i_jumeau_elimine, i_centres_exclus;

  // déclaration des visualiseurs de données par étapes (désactivés par défaut)
  //cvDisplayArray<double> *showHDPr, *showHDPi; 
  bool show_hdp_p;


  // ===========================================================================
  // pour debug
  
  bool d_keep_centers;
  string d_centers_log_file;

  // vecteur contenant les numéros d'angles traités
  std::vector<size_t> v_angles;
  // vecteur contenant les position des pics fouriers associés (0 à 4 * nxmax * nymax)
  std::vector<size_t> v_centers;

  
  // pour le débug: stack des spots de fourier 2D (autant de couches que d'angles traités)
  AIR_Volume<double> *d_VF_StackSpots_R;
  bool d_save_stack_spots_p;
  //const char* d_stack_spots_file;
  string d_stack_spots_file;
  // la même avec les découpes dans fourier
  AIR_Volume<double> *d_VF_StackWindows_R;
  bool d_save_stack_windows_p;
  string d_stack_windows_file;
  // la même avec les holos reconstruits centrés autour de position du pix
  AIR_Volume<double> *d_VF_StackHolograms_R;
  bool d_save_stack_holograms_p;
  string d_stack_holograms_file;
  bool d_show_peaks_p;


 public:
  void d_record_fourier_spots(string d_stack_file); // uniquement si un seul batch
  void d_record_fourier_windows(string d_stack_file); // uniquement si un seul batch

  void d_record_centers(size_t angle_max, string file);
  void d_record_final_holograms(string d_stack_file); // uniquement si un seul batch
  void d_show_peaks(bool b) { d_show_peaks_p = b; }

  // ===========================================================================
  // METHODES
  // ===========================================================================
 public:

  
  Holo_Process();
  ~Holo_Process();


  // --------------------------------------------------
  // série de fonctions à apeller pour configurer le calcul global
  // on peut les appeler dans n'importe quel ordre, sauf set_image_settings


  // --------------------
  // réglages de taille des images d'entrée (interférogrammes) et fenêtre de découpe
  // si autocenter_window, edge_x et edge_y sont calculés et écrasés
  bool set_image_settings(size_t window_edge_x, size_t window_edge_y, size_t window_size, bool autocenter_window);
  // la taille des images est déduite de la première image: 
  // DONC set_file_settings doit avoir été appelé, ou bien set_memory_settings (seul cas où l'ordre compte)
  

  // --------------------  
  // les volumes doivent avoir été alloués 
  bool set_output_volumes(AIR_Volume<RECON_TYPE>* V_out_fourier_R, AIR_Volume<RECON_TYPE>* V_out_fourier_I, AIR_Volume<unsigned short int>* V_out_supredon);


  // --------------------
  // réglage interférogrammes entrée: disque OU mémoire
  //
  // racine du nom de fichier complet d'une image (chemin relatif ou absolu) (e.g /tmp/session123-record)
  // proper: par défaut, %d %03d (false), %03d %03d (true)
  // si on demande une correction mod_ref, le fichier image correspondant doit être fourni et son existence sera vérifiée
  // si mode lecture interférogrammes sur disque
  bool set_file_settings(const char* images_radix, bool proper_names = false);

  // *OU* si mode lecture interférogrammes en mémoire
  //bool set_memory_settings(unsigned char** images_list, size_t image_dim_x, size_t image_dim_y);
  bool set_memory_settings(const AIR_Volume<unsigned char> *V_stack );
  // le nombre d'interférogrammes est déduit du nombre d'angles passé dans le paramétrage batch et de psh/pas psh.
  // mais la taille des images ne peut plus être déduite des fichiers, elle est donc passée


  // devrait être bientôt évincé: indique où sont les masques, dont tous sauf un sont générés maintenant
  bool set_path_settings(std::map<string, string> &map_paths);


  bool set_rayon_settings(size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, size_t xm0_limite, size_t ym0_limite, float rayon, float delta_zmax);

  
  // --------------------
  // opt: mod_ref correction
  bool set_modref_correction(const char* mod_ref_file);


  // fin série
  // --------------------------------------------------


  // --------------------------------------------------
  // accesseurs
  NREADER(size_t, p_Nxmax, get_nxmax);
  NREADER(size_t, p_Nymax, get_nymax);
  NREADER(size_t, p_image_dim_x, get_camera_dimx);
  NREADER(size_t, p_image_dim_y, get_camera_dimy);
  


  // --------------------------------------------------
  // fonctions optionnelles

  void set_psh(); // par défaut, valeurs déjà ok pour
  void set_off_axis_single(size_t off_axis_circle_cx, size_t off_axis_circle_cy, size_t off_axis_r);
  // ajouter le hors-axe au phase-shifting: préciser le centre du cercle des fréquences objet
  void set_off_axis_and_psh(size_t off_axis_circle_cx, size_t off_axis_circle_cy, size_t off_axis_r);
  
  // correction de focus par gradient
  // param: nb plans rétropropagés up / down (on doit calculer le double)
  void set_autofocus_gradient(size_t half_nb_planes);
    

  // fin série
  // --------------------------------------------------



  // fonction à appeller pour terminer la configuration du calcul global
  bool set_prepare();

  // lorsque tous les batchs sont passés (appel automatique si oublié)
  void set_free();


  void focus_debug_mode(string output_dir, size_t nb_angles, bool show_vals, bool show_pics) 
  {
    ASSERT(i_Refocalizer -> get_allocated());

    if (nb_angles)
      i_Refocalizer -> debug_mode(output_dir, nb_angles);
    
    // gestion des sorties de débog.
    if (show_vals)
      i_Refocalizer -> set_show_vals(true);  
    if (show_pics)
    i_Refocalizer -> set_show_pics();
    //  i_Refocalizer -> debug_mode("./", 20);
  }

  // --------------------------------------------------
  // une fois le calcul configuré ci-dessus, on peut alors lancer une série de batches sur tout ou partie des hologrammes

  // les fonction batch_prepare vont vérifier la présence des fichiers sur le disque (ou en mémoire)
  // et définir les n° des hologrammes disponibles correspondants
  // à disposition de batch_launch
  
  
  // définit un batch à partir d'un angle de départ (cf angle_last) et d'une taille escomptée de batch.
  // les candidats sont lus à la suite (pas de jump). max_error consécutives de lecture -> batch achevé 
  bool batch_prepare(size_t start_angle, size_t batch_size, size_t max_error);
  // on peut apeller cette fonction plusieurs fois sur des batchs consécutifs

  // tentative de batch temps réel
  // max_batch_size: batch traité automatiquement si angles dispo dépassent ce seuil
  // min_batch_size_start: au départ seulement, attend ce seuil d'angles dispo avant de traiter le batch
  bool batch_prepare_realtime(size_t max_batch_size, size_t min_batch_size_start, size_t sleep_ms = 100);
  // les appels de cette fonction sont successifs, car la classe garde en mémoire le n° du dernier angle traité
  
  // indique la taille du dernier batch configuré (zéro si aucun)
  size_t batch_query_size() const
  { return v_batch_angles.size(); }

  // indique le numéro du dernier angle traîté à l'instant
  size_t batch_query_angle_last_done() const
  { return i_angle_last_done; }

  // exécute le batch précédemment configuré
  void batch_launch();
  
 protected:
  // fonction-outil de batch_launch
  void hologram_compute(size_t current_angle);
 public:

  // fonction spéciale calculant uniquement l'hologramme recalé
  void hologram_show_debug(size_t current_angle, double* holo_recale_output);

  void display();

  
  void set_show_hdp() { show_hdp_p = true; }
  

  int* get_centres() const { return i_centre; }

  // --------------------------------------------------
  // fonctions-outil
 protected:
    
  template <class U>
    U sqr(U x) { return x * x; }

  // 1) phase-shifting: 4 images caméra -> 1 interférogramme complexe 512*512.
  void phase_shifting(size_t num_angle, double* plan_reel_out, double* plan_imag_out);
  // ALT: no-phase-shifting (HA pur):  (grande) image caméra -> la même.
  void no_phase_shifting(size_t num_angle, double* plan_reel_out, double* plan_imag_out);
  // NB: on recourt à phase_shifting quand on cumule PSH et HA
  
  // 2) recherche du pic non centré et normalisation
  size_t peak_normalize(size_t crop_size, double* in_fft_r, double* in_fft_i, \
			double* out_fft_r, double* out_fft_i);

  // 3) élimination du jumeau
  void twin_eliminate(int xmi, int ymi, double* io_fft_r, double* io_fft_i);
  void twin_eliminate_redux(int xmi, int ymi, double* io_fft_r, double* io_fft_i, unsigned char* cache_jumeau);
  
  // 4) mapping 3D
  void angle_mapping(size_t angle_num, int xmi, int ymi, double* in_fft_r, double* in_fft_i);

  
  void show_hdp();

  
  // fonctions outil génériques
  // la plupart devraient désormais être déportées dans vImg.hpp
  template <typename T>
    void extract_subImage(T* src, T* dst, size_t src_dimx, size_t src_dimy, size_t edge_x, size_t edge_y, size_t dst_dimx, size_t dst_dimy);
  void square_crop(double *src, double *dst, size_t src_size, size_t dst_size);
  void disk_crop(double *src, double *dst, size_t src_size, size_t c_x, size_t c_y, size_t rayon);
  // calcule le nom du fichier hologramme. Mode correct: tous numéros à %03d (pas supporté par défaut)
  void compute_hologram_filename(const char* images_radix, char* filename, size_t cpt_angle, size_t cpt_dp);

};


#endif
