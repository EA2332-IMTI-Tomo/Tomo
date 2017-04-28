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

#include "Camera.h"

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


// instructions sse, d'où usage de -msse dans g++
//#include "emmintrin.h"

#include "../recon_includes.h"
#include "../util.h"
#include "../Compute.h"
#include "../Holo_Process.h"
#include "../recon_parser.h"

// #include "FFTW_Image.h"
// #include "FFTW_Image.cc"
#include "FFTW_Image.hpp"




using namespace std;



/* -------------------------------------------------------------------------- */
// Usage 
/* -------------------------------------------------------------------------- */


static void
usage(int argc, char **argv) 
{

  if ((argc - 1) < 1)
    {
      printf("prend en argt un fichier de reconstruction\n");
      printf("affiche le flux caméra et l'hologramme recalé correspondant, selon paramètres de reconstruction\n");
      printf("mandatory: <file> configuration file \n");
      printf("opt: <dir> camera simulation directory \n");

      exit(EXIT_FAILURE);
    }
}



/* -------------------------------------------------------------------------- */
// 
/* -------------------------------------------------------------------------- */


int main(int argc, char *argv[])
{
  usage(argc, argv);
 
  bool simu_cam_p = false;
  char* simu_cam_dir;

  if ( argc - 1 >= 2)
    {
      simu_cam_p = true;
      simu_cam_dir = argv[2];
    }


  std::map<string, bool> map_bools;
  std::map<string, size_t> map_dims;
  std::map<string, string> map_paths;
  std::map<string, float> map_physics;
  // données expérimentales calculées
  std::map<string, float> map_experiment;


  init_vals(map_dims, map_bools, map_paths, map_physics);

    

  // ==========================================================
  // set parameters alike parsed configuration file
  // ==========================================================

  parse_file(argv[1], map_dims, map_bools, map_paths, map_physics);
  

  compute_vals(map_dims, map_physics, map_experiment);
  check_vals(map_dims, map_bools, map_paths, map_experiment);



  // ===========================================================================
  // pompé de compute1.cc en raccourci
  // ===========================================================================

  // --------------------------------------------------
  // récupération de paramètres parsés

  size_t window_edge_x = map_dims["window_edge_x"];
  size_t window_edge_y = map_dims["window_edge_y"];
  size_t image_dim_x = map_dims["image_dim_x"];
  size_t image_dim_y = map_dims["image_dim_y"];
  size_t xm0_limite = map_dims["xm0_limite"];
  size_t ym0_limite = map_dims["ym0_limite"];
  //   size_t angle_start = map_dims["angle_start"];
  //   size_t angle_last = map_dims["angle_end"];
  //   size_t angle_count = angle_last - angle_start + 1;

  size_t offaxis_cx = map_dims["off_axis_circle_cx"];
  size_t offaxis_cy = map_dims["off_axis_circle_cy"];
  size_t offaxis_r = map_dims["off_axis_radius"];

  const char* modref_file = map_paths["modref_file"].c_str();

  size_t cube_edge = (size_t) map_experiment["cube_edge"];
  size_t Nxmax = (size_t) map_experiment["nxmax"];
  size_t Nymax = Nxmax;
  size_t Nxmax_Rf = (size_t) map_experiment["nxmax_rf"];
  float rayon = map_experiment["rayon"];
  float delta_zmax = map_experiment["delta_zmax"];

  float pixel_size_recon = map_experiment["pixel_size_recon"];

  // peut être squizé sans affecter la qualité d'image, plus rapide.
  //   bool fourier_normalize_p = map_bools["fourier_normalize"];



  // --------------------------------------------------
  // Définition d'un processeur d'hologrammes 

  Holo_Process Hologram_Processor;
  ASSERT(map_dims["window_dim_x"] == map_dims["window_dim_y"]);
  Hologram_Processor.set_path_settings(map_paths);
  //Hologram_Processor.set_output_volumes(V_RealPart, V_ImagPart, V_SupRedon); // pas requis dans ce cas précis
  Hologram_Processor.set_image_settings(image_dim_x, image_dim_y, window_edge_x, window_edge_y, map_dims["window_dim_x"]);
  Hologram_Processor.set_rayon_settings(Nxmax, Nymax, Nxmax_Rf, xm0_limite, ym0_limite, rayon, delta_zmax);

    
  // --------------------------------------------------
  // renseignement du mode d'holographie
  bool psh_mode = map_bools["psh_mode"];
  bool off_axis_mode = map_bools["off_axis_mode"];


  if ( ! psh_mode && ! off_axis_mode)
    MSG_ASSERT(false, "You neither want PSH nor HA: can't perform holography!");

  if ( psh_mode && ! off_axis_mode)
    {
      cout << endl << "PSH MODE accounted (no off-axis)";
      Hologram_Processor.set_psh();
    }
  
  if ( off_axis_mode && psh_mode )
    {
      cout << endl << "OFF AXIS + PSH MODE accounted:";
      cout << off_axis_mode << " " << offaxis_cx << " " << offaxis_cy << " " << offaxis_r << endl; 
      Hologram_Processor.set_off_axis_and_psh( offaxis_cx, offaxis_cy, offaxis_r );
    }
  
  if ( off_axis_mode && ! psh_mode )
    {
      cout << endl << "OFF AXIS accounted (no PSH):";
      cout  << " " << offaxis_cx << " " << offaxis_cy << " " << offaxis_r << endl; 
      Hologram_Processor.set_off_axis_single( offaxis_cx, offaxis_cy, offaxis_r );
    }
  


  // ===========================================================================
  // on configure la caméra à lire des images de la taille indiquée dans le fichier de reconstruction

  Camera GEV_cam;
  GEV_cam.set_ImgDepth(_8b);
  
  if (simu_cam_p)
    {
      GEV_cam.set_simulation_mode( simu_cam_dir );
      cout << endl << "CAMERA SIMULATION FROM: " << simu_cam_dir;
      GEV_cam.set_Exposure(10); // l'exposition est obligatoire. à définir si simu
      //float cam_img_expo = GEV_cam.read_Exposure();
    }
  else
    cout << endl << "ETHERNET CAMERA SEEKED";

  size_t cam_x = Hologram_Processor.get_camera_dimx();
  size_t cam_y = Hologram_Processor.get_camera_dimy();
  GEV_cam.set_ImgX(cam_x);
  GEV_cam.set_ImgY(cam_y);

  // on connecte l'instance GEV_cam à la caméra physique. 
  ASSERT( GEV_cam.connect_firstCamera() );  
  GEV_cam.configure_snapshot(); 


  Mat* frame1 = GEV_cam.allocate_cv8image_Mat();
  Mat* frame2 = GEV_cam.allocate_cv8image_Mat();
  Mat* frame3 = GEV_cam.allocate_cv8image_Mat();
  Mat* frame4 = GEV_cam.allocate_cv8image_Mat();

  // IplImage* frame2 = GEV_cam . allocate_cv8image();



  // ===========================================================================
  // on alloue les données d'entrée pour un interférogramme et on le donne à la classe Hologram_Processor

  unsigned char **images_list;
  ARRAY_ALLOC( images_list, 4, unsigned char*);

  map_bools["read_from_mem"] = true;
  
  //   images_list[0] = (unsigned char*) frame1 -> imageData;
  images_list[0] = (unsigned char*) frame1 -> data;
  images_list[1] = (unsigned char*) frame2 -> data;
  images_list[2] = (unsigned char*) frame3 -> data;
  images_list[3] = (unsigned char*) frame4 -> data;

  Hologram_Processor.set_memory_settings(images_list);



  // --------------------------------------------------
  // fin de la préparation du processeur d'hologrammes

  Hologram_Processor.set_prepare();
  Hologram_Processor.display();
  //  Hologram_Processor.batch_launch();




  // ===========================================================================
  // ===========================================================================
  // on lit des images en boucle

  size_t nxmax = Hologram_Processor.get_nxmax();
  size_t nymax = Hologram_Processor.get_nymax();
  

  // on crée une image de la taille de l'hologramme recalé, et on y associe un visualiseur OpenCV
  double* holorecal_tmp;
  ARRAY_ALLOC(holorecal_tmp, 2 * nxmax * 2 * nymax, double);
  cvDisplayArray<double>* showRecal_tmp;
  showRecal_tmp = new cvDisplayArray<double>("Hologramme recalé", holorecal_tmp, 2 * nxmax, 2 * nymax);


  // on crée un visualiseur pour l'image caméra
  cvDisplayArray<unsigned char>* showCamera;
  showCamera = new cvDisplayArray<unsigned char>("Image Caméra", frame1);
  //showCamera = new cvDisplayArray<unsigned char>("Image Caméra",  (unsigned char*) frame1 -> data, cam_x, cam_y);

  char c;
  size_t iters = 0;

  showCamera -> set_tempo(10);

  while(true)
    {
      // MAJ des données sur images_list[0] (HA) ou images_list[0][1][2][3] (HDP)
      if ( off_axis_mode && ! psh_mode )
	GEV_cam.take_cv_snapshot(frame1);
      else
	{
	  GEV_cam.take_cv_snapshot(frame1);
	  GEV_cam.take_cv_snapshot(frame2);
	  GEV_cam.take_cv_snapshot(frame3);
	  GEV_cam.take_cv_snapshot(frame4);
	}

      Hologram_Processor.hologram_show_debug(1, holorecal_tmp);
      showRecal_tmp -> updateImage(); // quantifie, donc indispensable!
      showCamera -> updateImage(); // devrait être neutralisé ici			

      showRecal_tmp -> showImage();
      showCamera -> showImage();
      
      // un seul waitKey pour les deux, détection e/ESC/fermeture
      if (showCamera -> wait_continue_p('e')) break;
      // on détecte aussi la ferméture de l'autre fenêtre
      if (showRecal_tmp -> is_closed_p()) break;

      iters++;
      if (! (iters % 50) )
	cout << endl << iters;
    }
  

  // ===========================================================================

  

  Hologram_Processor.set_free();
  return EXIT_SUCCESS;
}





