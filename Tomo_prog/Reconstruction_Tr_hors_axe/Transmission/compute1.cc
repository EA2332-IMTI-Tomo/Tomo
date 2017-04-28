#ifdef SSE4
#include "FFTW_Volume_SSE4.h"
#include "AIR_Volume_SSE4.h"
#endif


void
compute_cpu_offline(std::map<string, float> my_map_experiment,	\
		    const string &images_radix, unsigned char **images_list, \
		    std::map<string, bool> my_map_bools, \
		    std::map<string, size_t> my_map_dims, \
		    std::map<string, string> my_map_paths)
{
  size_t window_edge_x = my_map_dims["window_edge_x"];
  size_t window_edge_y = my_map_dims["window_edge_y"];
  //   size_t image_dim_x = my_map_dims["image_dim_x"];
  //   size_t image_dim_y = my_map_dims["image_dim_y"];
  size_t xm0_limite = my_map_dims["xm0_limite"];
  size_t ym0_limite = my_map_dims["ym0_limite"];
  size_t angle_start = my_map_dims["angle_start"];
  size_t angle_last = my_map_dims["angle_end"];
  size_t angle_count = angle_last - angle_start + 1;

  size_t offaxis_cx = my_map_dims["off_axis_circle_cx"];
  size_t offaxis_cy = my_map_dims["off_axis_circle_cy"];
  size_t offaxis_r = my_map_dims["off_axis_radius"];

  const char* modref_file = my_map_paths["modref_file"].c_str();

  size_t cube_edge = (size_t) my_map_experiment["cube_edge"];
  size_t Nxmax = (size_t) my_map_experiment["nxmax"];
  size_t Nymax = Nxmax;
  size_t Nxmax_Rf = (size_t) my_map_experiment["nxmax_rf"];
  float rayon = my_map_experiment["rayon"];
  float delta_zmax = my_map_experiment["delta_zmax"];

  float pixel_size_recon = my_map_experiment["pixel_size_recon"];

  // peut être squizé sans affecter la qualité d'image, plus rapide.
  bool fourier_normalize_p = my_map_bools["fourier_normalize"];


  // ==========================================================
  // Allocation des volumes de données
  // ==========================================================


  AIR_Volume<RECON_TYPE> *V_RealPart, *V_ImagPart;

#ifdef SSE4
  V_RealPart = new AIR_Volume_SSE4<RECON_TYPE>(cube_edge, cube_edge, cube_edge, 1.0);
  V_ImagPart = new AIR_Volume_SSE4<RECON_TYPE>(cube_edge, cube_edge, cube_edge, 1.0);
#else
  V_RealPart = new AIR_Volume<RECON_TYPE>(cube_edge, cube_edge, cube_edge, 1.0);
  V_ImagPart = new AIR_Volume<RECON_TYPE>(cube_edge, cube_edge, cube_edge, 1.0);
#endif

  V_RealPart -> allocate();
  V_ImagPart -> allocate();
  


  // binding to the 2 real and imag volumes
  COMPLEX_Volume<RECON_TYPE> VC_Volume(*V_RealPart, *V_ImagPart);

  // on crée un volume calculable au format fftw. 
  FFTW_Volume<RECON_TYPE>* VF_Volume;

#ifdef SSE4
  VF_Volume = new FFTW_Volume_SSE4<RECON_TYPE>(cube_edge, cube_edge, cube_edge);
#else
  VF_Volume = new FFTW_Volume<RECON_TYPE>(cube_edge, cube_edge, cube_edge);
#endif

  VF_Volume -> use_wisdom_directory("/usr/local/phd/fftw_wisdom");
  VF_Volume -> set_threads_initialized(true); // car déjà initialisé avant, echec sinon
  VF_Volume -> set_nb_threads(my_map_dims["fftw_threads"]);
  VF_Volume -> set_report_time(true); //déconne



  // pour normaliser par rapport au nombre de fois où l'on remplit la frequence
  unsigned short int *sup_redon;
  AIR_Volume<unsigned short int> *V_SupRedon;
  
#ifdef SSE4
  V_SupRedon = new AIR_Volume_SSE4<unsigned short int>(cube_edge, cube_edge, cube_edge, 1.0);
#else
  V_SupRedon = new AIR_Volume<unsigned short int>(cube_edge, cube_edge, cube_edge, 1.0);
#endif

  V_SupRedon -> allocate();
  sup_redon = V_SupRedon -> get_data();


  // on va devoir conserver une copie circshiftée pour normaliser sur FFTW_Volume, à moins de modifier l'import en sauvage
  AIR_Volume<unsigned short int> *V_SupRedon_S;
  
#ifdef SSE4
  V_SupRedon_S = new AIR_Volume_SSE4<unsigned short int>(cube_edge, cube_edge, cube_edge, 1.0);
#else
  V_SupRedon_S = new AIR_Volume<unsigned short int>(cube_edge, cube_edge, cube_edge, 1.0);
#endif

  V_SupRedon_S -> allocate();
  


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
  
  
  // *******************************************************
  // lecture des hologrammes 



  vCHRONO_SET(readloop, "Boucle de lecture fichiers"); 
  vCHRONO_START(readloop);

  
  /*
    readLoop_angleImages(&V_RealPart, &V_ImagPart, sup_redon, centre,	\
    Nxmax, Nymax, Nxmax_Rf, window_edge_x, window_edge_y,		\
    image_dim_x, image_dim_y,				\
    xm0_limite, ym0_limite, rayon, delta_zmax,	\
    angle_start, angle_count, angle_jump,		\
    images_radix.c_str());
  */


  Holo_Process Hologram_Processor;
  ASSERT(my_map_dims["window_dim_x"] == my_map_dims["window_dim_y"]);

  if( ! my_map_bools["read_from_mem"] )
    Hologram_Processor.set_file_settings(images_radix.c_str(), my_map_bools["proper_numbers_mode"]);
  else
    MSG_ASSERT(false, "not implemented yet, check compute6.cc");

    //Hologram_Processor.set_memory_settings(images_list, my_map_dims["image_dim_x"], my_map_dims["image_dim_y"]);
  
  Hologram_Processor.set_path_settings(my_map_paths);
  Hologram_Processor.set_output_volumes(V_RealPart, V_ImagPart, V_SupRedon); 
  Hologram_Processor.set_image_settings(window_edge_x, window_edge_y, my_map_dims["window_dim_x"], my_map_bools["window_autocenter"]);
  Hologram_Processor.set_rayon_settings(Nxmax, Nymax, Nxmax_Rf, xm0_limite, ym0_limite, rayon, delta_zmax);
    

  if (my_map_bools["modref_correction_mode"])
    Hologram_Processor.set_modref_correction(modref_file);
  

  // debug
  string outdir_name(my_map_paths["OUTPUT_DIR"]);
  if (my_map_bools["debug_record_centers"])
    { string s1 = outdir_name + "/centre_dist.txt"; Hologram_Processor.d_record_centers(angle_count, s1); }
  if (my_map_bools["debug_record_fourier_spots"])
    { string s1 = outdir_name + "/svg_stack_fourier_spots"; Hologram_Processor.d_record_fourier_spots(s1); }
  if (my_map_bools["debug_record_fourier_windows"])
    { string s1 = outdir_name + "/svg_stack_fourier_windows"; Hologram_Processor.d_record_fourier_windows(s1); }
  if (my_map_bools["debug_record_final_holograms"])
    { string s1 = outdir_name + "/svg_stack_holos_recales"; Hologram_Processor.d_record_final_holograms( s1); }
  //Hologram_Processor.d_show_peaks(true);


  // --------------------------------------------------
  // choix du mode d'holographie
  bool psh_mode = my_map_bools["psh_mode"];
  bool off_axis_mode = my_map_bools["off_axis_mode"];


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

  // --------------------------------------------------


  if ( my_map_bools["autofocus_mode"] )
    {
      cout << endl << "AUTOFOCUS MODE accounted:";
      Hologram_Processor.set_autofocus_gradient( 20 );

      size_t count = 0;
      if (my_map_bools["autofocus_mode_debug"]) {
	cout << " + debug files on ";
	count = angle_count;
      }
      
      Hologram_Processor.focus_debug_mode(outdir_name, count, my_map_bools["autofocus_mode_show_vals"], my_map_bools["autofocus_mode_show_pics"]);
      
    }


  // nécessaire pour que les show_hdp y fonctionnent
  //Hologram_Processor.set_show_hdp();
  
  Hologram_Processor.set_prepare();
  Hologram_Processor.display();

  size_t angle_err_max = my_map_dims["angle_error_limit"];
  cout << endl << endl << "------------" << angle_start << " " << angle_count;
  MSG_ASSERT( Hologram_Processor.batch_prepare(angle_start, angle_count, angle_err_max), \
	      "failed to prepare a batch of files, error > 100 files"); 

  cout << endl << ">>> Boucle lecture hologrammes: "; cout.flush();
  Hologram_Processor.batch_launch();
  
  vCHRONO_STOP(readloop);

  //centre = Hologram_Processor.get_centres();
  //write_array<int>(size_t(4 * Nxmax * Nymax), centre, str_concat(g_OUTPUT_DIR, OUTPUT_CENTER_FILENAME));
  Hologram_Processor.set_free();


  // ==========================================================
  // Exécution de FFT 3D précédée et suivie de shifts circulaires
  // ==========================================================

  
  vCHRONO_ASTART(overall_fft, "whole FFT sequence (+shift & normalization) ");


  VF_Volume -> copy_from(VC_Volume); // copie et shift

  vCHRONO_ASTART(circ_supred, "circshift sup_redon (1 vol. int)");
  V_SupRedon -> circshift_to(*V_SupRedon_S);
  vCHRONO_STOP(circ_supred);

  vCHRONO_ASTART(norm_supred, "normalisation par sup_redon: (vol cplx / int)");
  *VF_Volume /= *V_SupRedon_S;
  vCHRONO_STOP(norm_supred);


  if(my_map_bools["save_fourier"])
    {
      //VF_Volume -> copy_to(VC_Volume); // copie et shift pour centrage sur 255,255,255
      string s_filename(my_map_paths["OUTPUT_DIR"] + "/"); 
      s_filename += my_map_paths["save_fourier_filename"];
      VF_Volume -> dump_to(s_filename); 
      //VC_Volume.change_files(s_filename);
      //VC_Volume.write_files();
      cout << endl << "sauvegarde fréquences de fourier (zéro en [0,0,0]) après normalisation sup_redon";
    }


  vCHRONO_ASTART(inner_fft, "inner fft");
  VF_Volume -> set_fourier_backward();
  vCHRONO_STOP(inner_fft);

  if (fourier_normalize_p)
    {
      VF_Volume -> set_fourier_normalize();
      cout << endl << "fourier post-normalization";
    }

  VF_Volume -> copy_to(VC_Volume); // copie et shift 

  vCHRONO_STOP(overall_fft);



  // *******************************************************
  // visualisation asynchrone (n'empiète pas sur la svg)
  
  // on peut sacrifier VC_Volume
  AIR_Volume<RECON_TYPE> *VC_View_R = VC_Volume.get_real_vol();
  AIR_Volume<RECON_TYPE> *VC_View_I = VC_Volume.get_imag_vol();
  cvDisplayVolume<RECON_TYPE> cv_VolSlicer("Reconstruction - Partie Imaginaire", VC_View_I);
  cv_VolSlicer.set_cv_tempo(30);
  cv_VolSlicer.enable_mouse();
  cv_VolSlicer.enable_histogram(512, 200, 256);
  cv_VolSlicer.enable_scale(1000000000 * pixel_size_recon, 10); // nm, um
  
  if(my_map_bools["end_visu_mode"])
    {
      MSG_ASSERT(my_map_dims["window_dim_x"] == my_map_dims["window_dim_y"], \
		 "la visu requiert images carrées pour le moment");
      
      boost::thread workerThread(BoostThreadFunc, &cv_VolSlicer);
    }
  

  // *******************************************************
  // sauvegardes sur disque


  vCHRONO_ASTART(savedisk, "sauvegarde sur disque des données");
	
  string s_file_radix(my_map_paths["OUTPUT_DIR"]); 
  s_file_radix = s_file_radix + "/" + my_map_paths["OUTPUT_RADIX"]; 
  VC_Volume.change_files(s_file_radix);
  VC_Volume.write_files();
  
	
  if (my_map_bools["save_sup_redon"])
    {
      string s_filename(my_map_paths["OUTPUT_DIR"] + "/"); 
      s_filename += my_map_paths["save_sup_redon_filename"];
      V_SupRedon -> change_files(s_filename);
      V_SupRedon -> write_hdr(); 
      V_SupRedon -> write_data();
    }
       
  
  vCHRONO_STOP(savedisk);    


  // *******************************************************
  // fin
  
  //delete centre;

  //libération memoire allouée pour les threads
  void fftw_cleanup_threads(void);
  
  vCHRONO_STOP(begin);


  // stats additionnelles:
  size_t count = V_SupRedon -> nzero_count();
  float coverage = (float) count / (float)( V_SupRedon -> get_x_dim() * V_SupRedon -> get_y_dim() * V_SupRedon -> get_z_dim());
  cout << endl << "% de voxels allumés dans SupRedon:" << coverage;
  cout << endl << "nb de voxels allumés dans SupRedon:" << count;
  cout << endl << endl << "recon:: waiting for C-c to exit";
  cout.flush();

  if (my_map_bools["end_visu_mode"])
    { char c = cvWaitKey(0); c = cvWaitKey(0); }
}
