/***************************************************************************
 *   Copyright (C) 2007 by mat (C) 2010 by vectra                          *
 *   matthieu.debailleul@uha.fr / bailleul@vectraproject.com               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


// p_ paramètre d'entrée du domaine de la physique ou de l'optique
// g_ variable globale visible dans tous les fichiers (qui "l'invitent" par extern)


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "vChrono.h" // requiert boost
#include "vChronos.h"


#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <exception>
#include <map>




#include "recon_includes.h"
#include "util.h"
#include "Compute.h"

#include "recon_parser.h"

// #include "FFTW_Image.h"
// #include "FFTW_Image.cc"
#include "FFTW_Image.hpp"


#ifdef CUDA
#include "cuda_Compute.h"
#endif

using namespace std;


/* -------------------------------------------------------------------------- */
// variables globales
/* -------------------------------------------------------------------------- */

ofstream g_log_file;


/* -------------------------------------------------------------------------- */
// Usage 
/* -------------------------------------------------------------------------- */


static void
usage(int argc, char **argv) 
{

#ifdef CUDA
  fprintf(stderr, "\n CUDA-enabled executable: %s", argv[0]);
#endif

  if ((argc - 1) == 0)
    {
      printf("Processes acquired planes images to compute a single result volume\n");
      printf("Values for I/O directories and processing parameters are read from mandatory ./config.txt file and can be redefined as follow:\n");
      printf("mandatory: <file> configuration file \n");
      printf("opt: -input_dir <dir> directory holding input images \n");
      printf("opt: -output_dir <dir> directory where results are written \n");
      printf("opt: -output_radix <dir> results volumes named radix_R/I.hdr/img \n");
      printf("opt: -masks <dir> directory holding masks \n");
      printf("opt: -in_radix <name> input images filename radix (program adds image number + extension) \n");
      printf("opt: -fftw_threads <int> number of threads to activate when computings ffts with fftw software library\n");      
      printf("opt: -image_dim_x(y) <int> to define the size of the input images (*ONLY* required if reading from memory stack) \n");      
      printf("opt: -window_dimx(y) <int> to define the window size on aforementioned images. This will define output size for reconstructed image \n");      
      printf("opt: -edge_x(y) <int> to define the upper left corner position of aforementioned frame \n");      
      printf("opt: -x(y)m0_limit <int> to set exclusion for some centers \n");      
      printf("opt: -first_angle <int> \n");      
      printf("opt: -final_angle <int> \n");      
      printf("opt: -inc_angle <int> \n");      
      printf("opt: -proper_numbers <0/1> proper numbers designates current %30d format, unproper old format is %3d\n");      
      printf("opt: -oil_index <float> \n");      
      printf("opt: -aperture <float> \n");  
      printf("opt: -lambda <float> to set wavelength \n");  
      printf("opt: -factor <int> pour définir grossissement du telan+objectif \n");  
      printf("opt: -pixel_size <float> \n");  
      printf("opt: -rf_val <float> \n");  
      printf("opt: -realtime_mode <dummy> real-time reconstruction and visualization during acquisition \n");
      printf("NB: on peut avoir PSH, PSH+HA, HA\n");
      printf("opt: -psh <bool> pour activer la modulation de phase (phase-shift).\n");
      printf("opt: -off_axis <bool> pour passer en mode hors-axe. Cette option vous engage à définir -circle*.\n");
      printf("opt: -circle_cx <c_x> \n");
      printf("opt: -circle_cy <c_y> : si off-axis, position du centre du cercle de fréquences, à mesurer lors du décalage du mirroir du faisceau de référence \n");
      printf("opt: -circle_r <r> : si off-axis, rayon du cercle, à mesurer lors du décalage du mirroir du faisceau de référence \n");
      printf("opt: -autofocus <dummy> pour effectuer une correction de focus par gradient (pour l'instant, seulement en mode cpu, offline, classique. Résultat inopérant ce jour) \n");
      printf("opt: -movie_on_slice <int>: active le mode film, enregistre uniquement la coupe donnée (def: 0) ) \n");
      printf("opt: -movie_every_angle <int>: active le mode film, saute le nombre précisé d'angles entre chaque enregistrement de coupe  (def: 1) ) \n");
      printf("opt: -movie_cut_plane <int>: active le mode film, précise le plan de coupe. 0=xy, 1=xz, 2=yz (def: xy) ) \n");
      printf("opt: -nv <dummy> pour désactiver la visualisation automatique à la fin pendant l'écriture volumes.\n");
      printf("opt: -mod_ref <mod_ref_file.pgm> correction du fond \n");  
      printf("opt: -fourier_normalize <dummy> correction de normalisation fftw, expérimentale \n");  

      exit(EXIT_FAILURE);
    }
}



/* -------------------------------------------------------------------------- */
// 
/* -------------------------------------------------------------------------- */


int main(int argc, char *argv[])
{
  usage(argc, argv);
  
  
  std::map<string, bool> map_bools;
  std::map<string, size_t> map_dims;
  std::map<string, string> map_paths;
  std::map<string, float> map_physics;
  // données expérimentales calculées
  std::map<string, float> map_experiment;


  init_vals(map_dims, map_bools, map_paths, map_physics);


  // ======================================================================
  // mode film

  //bool movie_mode = false;
  movieParams movie_settings{0, 1, _xy};
  

  // by default, we read all images at once after acquisition for blind reconstruction
  enum ReconType ReconMode = _offline; 
  // we can go online for reconstruction during acquisition (software
  // will sleep until hologram files actually arrive in input dir)
    
  

  // ==========================================================
  // set parameters alike parsed configuration file
  // ==========================================================

  parse_file(argv[1], map_dims, map_bools, map_paths, map_physics);
  argc -= 2; argv += 2;

  
  // ==========================================================
  // override parameters by command-line arguments
  // ==========================================================

  parse_argts(argc, argv, \
	      map_dims, map_bools, map_paths, map_physics, \
	      movie_settings, ReconMode);
  

  // ==========================================================
  // perform computations depending from parameter values
  // ==========================================================

  g_log_file.open( (map_paths["OUTPUT_DIR"] + "/" + map_paths["LOGFILE"]).c_str() );
  

  compute_vals(map_dims, map_physics, map_experiment);
  check_vals(map_dims, map_bools, map_paths, map_experiment);

  size_t cube_edge = (size_t) map_experiment["cube_edge"];
  size_t Nxmax = (size_t) map_experiment["nxmax"];
  size_t Nymax = Nxmax;
  size_t Nxmax_Rf = map_experiment["nxmax_rf"];
  float rayon = map_experiment["rayon"];
  float delta_zmax = map_experiment["delta_zmax"];


  // chemin de lecture fichiers entrée + radix de fichier (+ numfic.bmp)
  //   char Chemin[] = INPUT_DIR "/" INFILE_RADIX ;
  const string images_radix = map_paths["INPUT_DIR"] + (const string)"/" + map_paths["INFILE_RADIX"];

  // ===========================================================================


    
#ifdef DEBUG_MODE
  // debug output: valeurs saisies
  fprintf(stderr, "\n read: window_x %d | y %d | xmolim %d | y %d | edgex %d | y %d | imgx %d | y %d | %d | %d | %d | oil %f | ap %f | lam %f | fact %f | pixs %f | rf %f |", 
	  map_dims["window_dim_x"], map_dims["window_dim_y"], map_dims["xm0_limite"], map_dims["ym0_limite"], map_dims["window_edge_x"], map_dims["window_edge_y"], \
	  map_dims["image_dim_x"], map_dims["image_dim_y"], map_dims["angle_start"], map_dims["angle_end"], map_dims["angle_jump"], \
	  map_physics["oil_index"], map_physics["aperture"], map_physics["lambda"], map_physics["factor"], map_physics["pixel_size"], map_physics["focal_ratio"]);

  fprintf(stderr, "\n paths: | mask %s", map_paths["MASKPATH"].c_str());

  fprintf(stderr, "\n HDP mode %d | OFF-AXIS mode %d", map_bools["psh_mode"], map_bools["off_axis_mode"]);
#endif    
    
  

  size_t batch_size = 50;
  cout << endl << "**************************************************";





  
  // ===========================================================================
  // ordre lancement calcul
  // ===========================================================================


  // ===========================================================================
  // MODE CPU
  

#ifndef CUDA



  switch( ReconMode )
    {

      // ------------------------------------------------------------
      // on lit tous les fichiers en une fois: reconstruction rapide à partir de fichiers déjà sur disque
    case _offline:

      cout << endl << "offline non-batch mode CPU";
      map_bools["read_from_mem"] = false;
      

      if (! map_bools["read_from_mem"])
	compute_cpu_offline(map_experiment, images_radix, (unsigned char**) NULL, \
			    map_bools, map_dims, map_paths);
      else
	compute_memtest( map_experiment, images_radix, \
			 map_bools, map_dims, map_paths);

      break;


      // ------------------------------------------------------------
      // on affiche une reconstruction "temps-réel" (comme au dessus) mais, en plus, 
      // on est capable d'attendre le début de l'acquisition et que des nouveaux angles arrivent
      // attention, version CPU ici
    case _online_batch:

      cout << endl << "online batch mode CPU";
      compute_cpu_online_batch(map_experiment, images_radix, (unsigned char**) NULL, \
			       map_bools, map_dims, map_paths);
      break;


      // ------------------------------------------------------------
      // on remplace la reconstruction du volume complet (r, i) par une stack des reconstructions d'une coupe fixe 
      // chaque coupe est extraite du volume complet (r, i) calculé à divers nombres d'angles 
      // (on spécifie un intervalle entre angles, 1 par défaut)
    case _offline_movie:
      
      cout << endl << "offline movie: slice=" << movie_settings.movie_on_slice << " every " << movie_settings.movie_angle_gap << " angle";
      switch( movie_settings.movie_cut_view ) {
      case _xy: 
	cout << ", xy cut"; break;
      case _xz:
	cout << ", xz cut"; break;
      case _yz:
	cout << ", yz cut"; break;
      }
      cout << endl;

      compute_cpu_offline_movie( cube_edge, Nxmax, Nymax, Nxmax_Rf,	\
				 rayon, delta_zmax,			\
				 movie_settings.movie_cut_view, movie_settings.movie_on_slice, movie_settings.movie_angle_gap, \
				 images_radix, map_bools, map_dims, map_paths);
      
      break;

      
      // ------------------------------------------------------------
      //  DEPRECATED CPU SECTION      
      // ------------------------------------------------------------


      // ------------------------------------------------------------
      // afin de tester le mode online, on considère que tous les fichiers à traiter 
      // sont déjà sur disque, mais on les traite et reconstruit par lots
    case _offline_batch:
      cout << endl << "offline batch mode CPU";
      MSG_ASSERT(false, "retiré car inutile");
      /*
	compute_cpu_offline_batch( cube_edge, Nxmax, Nymax, Nxmax_Rf,	\
	rayon, delta_zmax,			\
	images_radix, 50, map_bools, map_dims);
      */
      break;
      
      
      // ------------------------------------------------------------
    default:
      MSG_ASSERT(false, "unsupported mode in CPU version");
      break;
    }


#else


  // ===========================================================================
  // MODE GPU
  
  switch( ReconMode )
    {
          
    case _online_batch:
      cout << endl << "online batch mode GPU";  //
      compute_gpu_online_batch( cube_edge, Nxmax, Nymax, Nxmax_Rf,	\
				window_edge_x, window_edge_y, lolno, lolno,\  // image_dim_x, image_dim_y, \
				xm0_limite, ym0_limite, rayon, delta_zmax, \
				50, images_radix);
      
      break;

    default:
      MSG_ASSERT(false, "unsupported mode in GPU version");
    }



#endif

  g_log_file.close();

  return EXIT_SUCCESS;
}






// ===========================================================================
// ===========================================================================  
// ===========================================================================


























// cube_edge: taille du volume cubique du volume final et du volume fourier 3D
// image_dim_x/y: taille des images de la caméra, à savoir les interférogrammes
// window_edge_x: position de la fenêtre de découpe dans les images précédemment définies
// window_dim_x/y: paramètres cachés (à tort) indicant la taille de la fenêtre de découpe
//                   supposée carrée ici
// angle_start: premier angle à traiter (4 interférogrammes par angle)
// images_radix: racine du nom des interférogrammes (nom fichier, pas répertoire)
// input_dir: paramètre caché, nom du répertoire 




/*
// template instanciation
void faitchier()
{
double* plan_reel;
ARRAY_ALLOC(plan_reel, 256 * 256, double);
double* plan_imag;
ARRAY_ALLOC(plan_imag, 256 * 256, double);

  
FFTW_Image<double, fftw_complex> *VF_Fourier2D;  
VF_Fourier2D = new FFTW_Image<double, fftw_complex>(256, 256, 3);
VF_Fourier2D -> allocate();
  
VF_Fourier2D -> import_from(plan_reel, plan_imag, true);
VF_Fourier2D -> set_fourier_forward();
VF_Fourier2D -> set_fourier_backward();
VF_Fourier2D -> export_to(plan_reel, plan_imag, true); 
  
VF_Fourier2D -> unallocate();

}


// ------------------------------------------------------------
// comme offline classique, avec un zoom centré en plus
case _offline_zoom:
      
cout << endl << "offline non-batch mode CPU + centered zoom (NB: c'était une fausse-bonne idée) ";
compute_cpu_zoom_offline( cube_edge, Nxmax, Nymax, Nxmax_Rf,		\
window_edge_x, window_edge_y, image_dim_x, image_dim_y,	\
xm0_limite, ym0_limite, rayon, delta_zmax,	\
angle_start, angle_count, angle_jump,	\
off_axis_mode, off_axis_circle_cx, off_axis_circle_cy, \
autofocus_mode, \
images_radix, proper_numbers_mode, modref_asked_mode, modref_file);
break;

*/


/*
  #ifdef CUDA
  // doit-on employer cuda pour la TF3D (si mémoire vidéo le permet)
  bool fft3d_cuda;
  #endif

  // debug
  bool read_saved_volume;
  // produire ou non dans /ramdisk des volumes de contrôle
  bool save_sup_redon;
*/



/* DEPRECATED

   #ifdef CUDA
   printf("opt: DEPRECATED -cuda_3d <0/1> to use CUDA for computing the 3D fft if graphic card has enough memory\n");
   #endif
   printf("opt: -zoom_mode <dummy> pour effectuer un zoom centré x2 dans Fourier (pour l'instant, seulement en mode cpu, offline, classique. Inutile ce jour.) \n");
*/
