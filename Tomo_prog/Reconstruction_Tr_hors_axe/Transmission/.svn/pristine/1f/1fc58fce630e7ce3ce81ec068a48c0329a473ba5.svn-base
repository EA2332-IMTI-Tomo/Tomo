#include "util_Image.h"

bool gl_proper_names_p;
size_t gl_skipend;



void 
local_compute_hologram_filename(const char* images_radix, char* filename, size_t cpt_angle, size_t cpt_dp)
{
  ASSERT((cpt_dp >= 1) && (cpt_dp <= 4));
  if (! gl_proper_names_p)
    sprintf(filename, "%s%d-%03d." INPUT_FORMAT, images_radix, cpt_angle, cpt_dp);
  else
    sprintf(filename, "%s%03d-%03d." INPUT_FORMAT, images_radix, cpt_angle, cpt_dp);
}



// Fait une reconstruction classique, sauf que les images sont pré-chargées en mémoire puis passées à la reconstruction en mémoire.
// but: tester que ça marche et que ça va assez vite.

// res: 3.80 contre 4.20 
// pas encore au point: l'image est altérée!!


void
compute_memtest(std::map<string, float> my_map_experiment,	\
		const string &images_radix, 
		std::map<string, bool> my_map_bools,	\
		std::map<string, size_t> my_map_dims,	\
		std::map<string, string> my_map_paths)
{
  size_t window_edge_x = my_map_dims["window_edge_x"];
  size_t window_edge_y = my_map_dims["window_edge_y"];
  size_t window_size = my_map_dims["window_dim_x"];
  size_t xm0_limite = my_map_dims["xm0_limite"];
  size_t ym0_limite = my_map_dims["ym0_limite"];
  size_t angle_start = my_map_dims["angle_start"];
  size_t angle_last = my_map_dims["angle_end"];
  size_t angle_count = angle_last - angle_start + 1;

  size_t offaxis_cx = my_map_dims["off_axis_circle_cx"];
  size_t offaxis_cy = my_map_dims["off_axis_circle_cy"];
  size_t offaxis_r = my_map_dims["off_axis_radius"];

  size_t cube_edge = (size_t) my_map_experiment["cube_edge"];
  size_t Nxmax = (size_t) my_map_experiment["nxmax"];
  size_t Nymax = Nxmax;
  size_t Nxmax_Rf = (size_t) my_map_experiment["nxmax_rf"];
  float rayon = my_map_experiment["rayon"];
  float delta_zmax = my_map_experiment["delta_zmax"];

  gl_proper_names_p = my_map_bools["proper_numbers_mode"];
  gl_skipend = 10;

  // ==========================================================
  // Allocation des volumes de données
  // ==========================================================


  AIR_Volume<RECON_TYPE> V_RealPart(cube_edge, cube_edge, cube_edge, 1.0);
  AIR_Volume<RECON_TYPE> V_ImagPart(cube_edge, cube_edge, cube_edge, 1.0);

  V_RealPart.allocate();
  V_ImagPart.allocate();


  // binding to the 2 real and imag volumes
  COMPLEX_Volume<RECON_TYPE> VC_Volume(V_RealPart, V_ImagPart);

  // on crée un volume calculable au format fftw. 
  FFTW_Volume<RECON_TYPE> VF_Volume(cube_edge, cube_edge, cube_edge);
  VF_Volume.use_wisdom_directory("/usr/local/phd/fftw_wisdom");
  VF_Volume.set_threads_initialized(true); // car déjà initialisé avant, echec sinon
  VF_Volume.set_nb_threads(my_map_dims["fftw_threads"]);
  VF_Volume.set_report_time(true); //déconne



  // pour normaliser par rapport au nombre de fois où l'on remplit la frequence
  unsigned short int *sup_redon;


  AIR_Volume<unsigned short int> V_SupRedon(cube_edge, cube_edge, cube_edge, 1.0);
  V_SupRedon.allocate();
  sup_redon = V_SupRedon.get_data();


  // on va devoir conserver une copie circshiftée pour normaliser sur FFTW_Volume, à moins de modifier l'import en sauvage
  AIR_Volume<unsigned short int> V_SupRedon_S(cube_edge, cube_edge, cube_edge, 1.0);
  V_SupRedon_S.allocate();



  // pour mettre la position des centres translatés, on crée une
  // variable 2D de la taille d'un plan apres tomo

  int *centre;
  //centre = new int[4*Nxmax*Nymax];


    
  // ==========================================================
  // Début du Programme effectif
  // ==========================================================

    

  vCHRONO_SET(begin, "Lancement programme"); 
  vCHRONO_START(begin);


  birdy();
  
  // initialize multithreaded fftw 
  if (my_map_dims["fftw_threads"] > 1)
    ASSERT(fftw_init_threads());
  



  // --------------------------------------------------
  // choix du mode d'holographie
  bool psh_mode = my_map_bools["psh_mode"];
  bool off_axis_mode = my_map_bools["off_axis_mode"];

  // *******************************************************
  // tranfert des hologrammes en mémoire pour test code

  char* holo_filename;
  ARRAY_ALLOC(holo_filename, 1000, char);


  // *******************************************************
  // on crée un dummy_batch_processor pour extraire le nb d'images de la liste

  Holo_Process dummy;
  ASSERT(my_map_dims["window_dim_x"] == my_map_dims["window_dim_y"]);
  dummy.set_output_volumes(&V_RealPart, &V_ImagPart, &V_SupRedon); 
  dummy.set_file_settings(images_radix.c_str(), my_map_bools["proper_numbers_mode"]);
  dummy.set_image_settings(window_edge_x, window_edge_y, my_map_dims["window_dim_x"], my_map_bools["window_autocenter"]);
  dummy.set_rayon_settings(Nxmax, Nymax, Nxmax_Rf, xm0_limite, ym0_limite, rayon, delta_zmax);
  dummy.set_path_settings(my_map_paths);
  //skip modref
  

  if ( ! psh_mode && ! off_axis_mode)
    MSG_ASSERT(false, "You neither want PSH nor HA: can't perform holography!");

  if ( psh_mode && ! off_axis_mode)
    {
      cout << endl << "PSH MODE accounted, no HA combined:";
      dummy.set_psh();
    } 
  
  if ( off_axis_mode && psh_mode )
    {
      cout << endl << "OFF AXIS + PSH MODE accounted:";
      cout << off_axis_mode << " " << offaxis_cx << " " << offaxis_cy << endl; 
      dummy.set_off_axis_and_psh( offaxis_cx, offaxis_cy, offaxis_r );
    }
  
  if ( off_axis_mode && ! psh_mode )
    {
      cout << endl << "OFF AXIS and NO PSH MODE accounted:";
      cout  << " " << offaxis_cx << " " << offaxis_cy << endl; 
      dummy.set_off_axis_single( offaxis_cx, offaxis_cy, offaxis_r );
    }



  
  dummy.set_prepare();
  MSG_ASSERT( dummy.batch_prepare(angle_start, angle_count, 100), "eedvd");
  

  // extraction terminée!!!
  size_t nb_angles_effectifs =  dummy.v_batch_angles.size();
  cout << endl << endl << "nombre d'angles considérés corrects" << nb_angles_effectifs << endl;


  // *******************************************************
  // maintenant on remplit une liste d'images avec les images inférées par le dummy_batch_processor
  

  // if( ! my_map_bools["read_from_mem"] )
  /*
  unsigned char** images_list;
  ARRAY_ALLOC(images_list, nb_angles_effectifs * 4, unsigned char*);
  for (size_t j = 0; j < nb_angles_effectifs * 4; j++)
    ARRAY_ALLOC(images_list[j], image_dim_x * image_dim_y, unsigned char);
  */


  // la pile d'hologrammes ne contient que la fenêtre centrale!
  size_t nb_holograms = nb_angles_effectifs;
  if ( psh_mode && ! off_axis_mode ) nb_holograms *= 4;
  AIR_Volume<unsigned char> V_hologramStack( my_map_dims["window_dim_x"], my_map_dims["window_dim_y"], nb_holograms);
  V_hologramStack.allocate();
  unsigned char* stack_data = V_hologramStack.get_data();


  /*
  for (size_t cursor = 0; cursor < nb_angles_effectifs; cursor++)
    {
      size_t current_angle = dummy.v_batch_angles[cursor];
      cout << " | " << current_angle;

      for (size_t k = 0; k < 4; k++)
	{
	  local_compute_hologram_filename(images_radix.c_str(), holo_filename, current_angle, k+1);
	  // cout << endl << "nb angles eff:" << nb_angles_effectifs << "holo_filename:" << holo_filename;
	  remplit_tableau_cv(images_list[cursor * 4 + k], holo_filename, image_dim_x, image_dim_y);
	}
    }
  */


  vImg<unsigned char> I_reader( dummy.get_camera_dimx(), dummy.get_camera_dimy() );
  vImg<unsigned char> I_cropper( stack_data, my_map_dims["window_dim_x"], my_map_dims["window_dim_y"] );
  size_t img_size = my_map_dims["window_dim_x"] * my_map_dims["window_dim_y"];

  // remplissage effectif de la pile
  if ( psh_mode && ! off_axis_mode ) // décalage de phase
    {
      for (size_t i = 0; i < nb_holograms; i++)
	{
	  unsigned char* write_ptr = V_hologramStack[i]; //stack_data + i * img_size;
	  for (size_t j = 0; j < 4; j++)
	    {
	      write_ptr += j = img_size;
	      I_cropper.update_pointer_single( write_ptr );
	      local_compute_hologram_filename(images_radix.c_str(), holo_filename, i, j+1);
	      I_reader.reload_pgm( holo_filename );
	      I_cropper.copy_from_bigger( I_reader );
	    }
	}
    }
  else // mode hors-axe pur
    {
      for (size_t i = 0; i < nb_holograms; i++)
	{
	  unsigned char* write_ptr = V_hologramStack[i];  
	  I_cropper.update_pointer_single( write_ptr );
	  local_compute_hologram_filename(images_radix.c_str(), holo_filename, i, 1);
	  I_reader.reload_pgm( holo_filename );
	  I_cropper.copy_from_bigger( I_reader );
	}
    }
     

  cout << endl << "============================================================";
  cout << endl << "PILE D'IMAGES REMPLIE";
  cout << endl << "============================================================";


  // *******************************************************
  // lecture des hologrammes 


  vCHRONO_SET(readloop, "Boucle de lecture fichiers"); 
  vCHRONO_START(readloop);



  Holo_Process Hologram_Processor;
  ASSERT(my_map_dims["window_dim_x"] == my_map_dims["window_dim_y"]);
  //Hologram_Processor.set_memory_settings(images_list, my_map_dims["image_dim_x"], my_map_dims["image_dim_y"]);
  Hologram_Processor.set_memory_settings( &V_hologramStack );

  Hologram_Processor.set_output_volumes(&V_RealPart, &V_ImagPart, &V_SupRedon); 
  Hologram_Processor.set_rayon_settings(Nxmax, Nymax, Nxmax_Rf, xm0_limite, ym0_limite, rayon, delta_zmax);
  Hologram_Processor.set_path_settings(my_map_paths);

  // debug
  //   Hologram_Processor.d_record_centers(nb_angles_effectifs , "./centre_dist.txt");
  //   Hologram_Processor.d_record_fourier_spots("./svg_stack_fourier_spots");


  if ( ! psh_mode && ! off_axis_mode)
    MSG_ASSERT(false, "You don't either want PSH or HA: can't perform holography!");

  if ( psh_mode && ! off_axis_mode)
    {
      cout << endl << "PSH MODE accounted, no HA combined:";
      Hologram_Processor.set_psh();
    }
  
  if ( off_axis_mode && psh_mode )
    {
      cout << endl << "OFF AXIS + PSH MODE accounted:";
      cout << off_axis_mode << " " << offaxis_cx << " " << offaxis_cy << endl; 
      Hologram_Processor.set_off_axis_and_psh( offaxis_cx, offaxis_cy, offaxis_r );
    }
  
  if ( off_axis_mode && ! psh_mode )
    {
      cout << endl << "OFF AXIS and NO PSH MODE accounted:";
      cout  << " " << offaxis_cx << " " << offaxis_cy << endl; 
      Hologram_Processor.set_off_axis_single( offaxis_cx, offaxis_cy, offaxis_r );
    }

  // --------------------------------------------------


  if ( my_map_bools["autofocus_mode"] )
    {
      cout << endl << "AUTOFOCUS MODE accounted:";
      Hologram_Processor.set_autofocus_gradient( 20 );
      /* if (my_map_bools["autofocus_mode_debug"]) {
	 cout << " + debug files on ";
	 Hologram_Processor.focus_debug_mode("/ramdisk", angle_count ); 
      
	 }*/
    }

  
  Hologram_Processor.set_prepare();
  Hologram_Processor.display();

  MSG_ASSERT( Hologram_Processor.batch_prepare(angle_start, nb_angles_effectifs, 100), \
	      "failed to prepare a batch of files, error > 100 files"); 

  cout << endl << ">>> Boucle lecture hologrammes: "; cout.flush();
  Hologram_Processor.batch_launch();
  
  vCHRONO_STOP(readloop);

  //centre = Hologram_Processor.get_centres();
  //write_array<int>(size_t(4 * Nxmax * Nymax), centre, str_concat(g_OUTPUT_DIR, OUTPUT_CENTER_FILENAME));
  //Hologram_Processor.set_free();

  // ==========================================================
  // Exécution de FFT 3D précédée et suivie de shifts circulaires
  // ==========================================================

  
  vCHRONO_ASTART(overall_fft, "whole FFT");


  VF_Volume.copy_from(VC_Volume); // copie et shift

  vCHRONO_ASTART(circ_supred, "circshift sup_redon (1 vol. int)");
  V_SupRedon.circshift_to(V_SupRedon_S);
  vCHRONO_STOP(circ_supred);

  vCHRONO_ASTART(norm_supred, "normalisation par sup_redon: (vol cplx / int)");
  VF_Volume /= V_SupRedon_S;
  vCHRONO_STOP(norm_supred);

  vCHRONO_ASTART(inner_fft, "inner fft");
  VF_Volume.set_fourier_backward();
  vCHRONO_STOP(inner_fft);
  VF_Volume.copy_to(VC_Volume); // copie et shift

  vCHRONO_STOP(overall_fft);



  // *******************************************************
  // visualisation asynchrone (n'empiète pas sur la svg)
    
  AIR_Volume<RECON_TYPE> *VC_View_I = VC_Volume.get_imag_vol();
  cvDisplayVolume<RECON_TYPE> cv_VolSlicer("Reconstruction - Partie Imaginaire", VC_View_I);
  cv_VolSlicer.set_cv_tempo(30);
  cv_VolSlicer.enable_mouse();
  cv_VolSlicer.enable_histogram(512, 200, 256);

  if(my_map_bools["end_visu_mode"])
    {
      MSG_ASSERT(my_map_dims["window_dim_x"] == my_map_dims["window_dim_y"], \
		 "la visu requiert images carrées pour le moment");
      
      boost::thread workerThread(BoostThreadFunc, &cv_VolSlicer);
    }

  // *******************************************************
  // sauvegardes sur disque

  //bool g_save_sup_redon = false;

  vCHRONO_ASTART(savedisk, "sauvegarde sur disque des données");
	
  string s_file_radix(my_map_paths["OUTPUT_DIR"]); 
  s_file_radix = s_file_radix + "/" + my_map_paths["OUTPUT_RADIX"]; 
  VC_Volume.change_files(s_file_radix);
  VC_Volume.write_files();
	
  if (my_map_bools["save_sup_redon"])
    {
      string s_filename(my_map_paths["OUTPUT_DIR"]); 
      s_filename += my_map_paths["save_sup_redon_filename"];
      V_SupRedon.change_files(s_filename.c_str());
      V_SupRedon.write_hdr(); V_SupRedon.write_data();
    }
	
  
  vCHRONO_STOP(savedisk);
      


  // *******************************************************
  // fin
  
  //delete centre;

  //libération memoire allouée pour les threads
  void fftw_cleanup_threads(void);
  
  vCHRONO_STOP(begin);


  // stats additionnelles:
  size_t count = V_SupRedon.nzero_count();
  float coverage = (float) count / (float)( V_SupRedon.get_x_dim() * V_SupRedon.get_y_dim() * V_SupRedon.get_z_dim());
  cout << endl << "% de voxels allumés dans SupRedon:" << coverage;
  cout << endl << "nb de voxels allumés dans SupRedon:" << count;
  cout << endl << endl << "recon:: waiting for C-c to exit";
  cout.flush();

  if (my_map_bools["end_visu_mode"])
    char c = cvWaitKey(0);
}
