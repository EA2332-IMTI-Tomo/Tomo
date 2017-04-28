void
compute_cpu_offline_movie(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
			  float rayon, float delta_zmax,		\
			  SlicePlanType cut, size_t slice_saved, size_t slice_gap, \
			  const string &images_radix, \
			  std::map<string, bool> my_map_bools, \
			  std::map<string, size_t> my_map_dims, \
			  std::map<string, string> my_map_paths)
{
  size_t window_edge_x = my_map_dims["window_edge_x"];
  size_t window_edge_y = my_map_dims["window_edge_y"];
  size_t image_dim_x = my_map_dims["image_dim_x"];
  size_t image_dim_y = my_map_dims["image_dim_y"];
  size_t xm0_limite = my_map_dims["xm0_limite"];
  size_t ym0_limite = my_map_dims["ym0_limite"];
  size_t angle_first = my_map_dims["angle_start"];
  size_t angle_last = my_map_dims["angle_end"];

  
  size_t offaxis_cx = my_map_dims["off_axis_circle_cx"];
  size_t offaxis_cy = my_map_dims["off_axis_circle_cy"];
  size_t offaxis_r = my_map_dims["off_axis_radius"];
  
  // ca c'est bon, virer ces arguments de merde de l'entête de la fonction plutôt
  //
  //   size_t cube_edge = (size_t) my_map_experiment["cube_edge"];
  //   size_t Nxmax = (size_t) my_map_experiment["nxmax"];
  //   size_t Nymax = Nxmax;
  //   size_t Nxmax_Rf = (size_t) my_map_experiment["nxmax_rf"];
  //   float rayon = my_map_experiment["rayon"];
  //   float delta_zmax = my_map_experiment["delta_zmax"];

  size_t batch_size = slice_gap;


  // ==========================================================
  // Allocation des volumes de données
  // ==========================================================


  // MSG_ASSERT(false, "au boulot hollywoud");
  cout << endl << "MOVIE RECORD MODE" << endl;

  // Volumes d'entrée ?
  AIR_Volume<RECON_TYPE> V_RealPart(cube_edge, cube_edge, cube_edge, 1.0);
  AIR_Volume<RECON_TYPE> V_ImagPart(cube_edge, cube_edge, cube_edge, 1.0);

  V_RealPart.allocate();
  V_ImagPart.allocate();

  // binding to the 2 real and imag volumes
  COMPLEX_Volume<RECON_TYPE> VC_Volume(V_RealPart, V_ImagPart);
  

  // Volume recontruit affichable ?
  // linéaire et 1.0 par défaut
  COMPLEX_Volume<RECON_TYPE> VC_View(cube_edge, cube_edge, cube_edge);
  VC_View.allocate();


  // on crée un volume calculable au format fftw. 
  FFTW_Volume<RECON_TYPE> VF_Volume(cube_edge, cube_edge, cube_edge);
  VF_Volume.use_wisdom_directory("/usr/local/phd/fftw_wisdom");
  VF_Volume.set_threads_initialized(true); // car déjà initialisé avant, echec sinon
  VF_Volume.set_nb_threads(my_map_dims["fftw_threads"]);
  VF_Volume.set_report_time(false); 



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
  centre = new int[4*Nxmax*Nymax];


    
  // ==========================================================
  // Début du Programme effectif
  // ==========================================================

    
 
  vCHRONO_SET(begin, "Lancement programme"); 
  vCHRONO_START(begin);


  birdy();
  
  // initialize multithreaded fftw 
  if (my_map_dims["fftw_threads"] > 1)
    ASSERT(fftw_init_threads());
  
  
  // ---------------------------------------------------
  // PRÉPARATION visu
  ASSERT(my_map_dims["window_dim_x"] == my_map_dims["window_dim_y"]);

  AIR_Volume<RECON_TYPE> *VC_View_R = VC_View.get_real_vol();
  AIR_Volume<RECON_TYPE> *VC_View_I = VC_View.get_imag_vol();

  

  // *******************************************************
  // lecture des hologrammes 



  // hophop, faut repasser en mode offline non batch , la
  
  Holo_Process Hologram_Processor;
  Hologram_Processor.set_path_settings(my_map_paths);
  Hologram_Processor.set_output_volumes(&V_RealPart, &V_ImagPart, &V_SupRedon); 
  Hologram_Processor.set_image_settings(window_edge_x, window_edge_y, my_map_dims["window_dim_x"], my_map_bools["window_autocenter"]);
  Hologram_Processor.set_rayon_settings(Nxmax, Nymax, Nxmax_Rf, xm0_limite, ym0_limite, rayon, delta_zmax);


 // --------------------------------------------------
  // choix du mode d'holographie
  bool psh_mode = my_map_bools["psh_mode"];
  bool off_axis_mode = my_map_bools["off_axis_mode"];


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


  Hologram_Processor.set_prepare();

  //Hologram_Processor.display();


  // Je pense que ça foire s'il y a des sauts d'angle dans les enregistrements
  size_t nb_batchs = angle_last / batch_size;
  size_t angle_lastread = angle_first;


  
  // Volume final: stack de slices
  // linéaire et 1.0 par défaut
  COMPLEX_Volume<RECON_TYPE> VC_Movie(cube_edge, cube_edge, nb_batchs);
  VC_Movie.allocate();
  AIR_Volume<RECON_TYPE> *V_Movie_R = VC_Movie.get_real_vol();
  AIR_Volume<RECON_TYPE> *V_Movie_I = VC_Movie.get_imag_vol();
  //VC_View_R,I

    
     
  for (size_t batch_num = 0; batch_num < nb_batchs; batch_num++)
    {
      cout << endl << endl << "movie batch n°" << batch_num + 1 << " || first angle: " << angle_lastread << "/ last angle: " << angle_lastread + batch_size;
      cout.flush();

      if(!  Hologram_Processor.batch_prepare(angle_lastread, batch_size, 4) )
	{ 
	  cout << endl << "batch break condition at n°" << angle_lastread ;
	  break; 
	}
      else
	Hologram_Processor.batch_launch();
      
      // ==========================================================
      // Exécution de FFT 3D précédée et suivie de shifts circulaires (2.5s cpu)
      // ==========================================================

      
 
      // VC_Volume /= V_SupRedon; non, destructif! on est pas en GPU

      VF_Volume.copy_from(VC_Volume); // copie et shift
      
      V_SupRedon.circshift_to(V_SupRedon_S);
      VF_Volume /= V_SupRedon_S;

      VF_Volume.set_fourier_backward();
 
      VF_Volume.copy_to(VC_View); // copie et shift
 
      //cv_VolSlicer.notify_data_changed();

      // ici; copie des slices
      switch( cut )
	{
	case _xy:
	  VC_View_R -> copy_slice( slice_saved, *V_Movie_R, batch_num );
	  VC_View_I -> copy_slice( slice_saved, *V_Movie_I, batch_num );
	  break;
	case _xz:
	  MSG_ASSERT(false, "pas implémenté encore");
	  break;
	case _yz:
	  VC_View_R -> copy_yz_slice( slice_saved, *V_Movie_R, batch_num );
	  VC_View_I -> copy_yz_slice( slice_saved, *V_Movie_I, batch_num );
	  break;
	}


      angle_lastread += batch_size;
      
    }
  // FIN BOUCLE BATCHES
  // ===========================================================================


  // *******************************************************
  // sauvegardes sur disque   // ici, remplacer VC_View par la stack des coupes sauvées


  vCHRONO_ASTART(savedisk, "sauvegarde sur disque des données");

  
  string s_file_radix(my_map_paths["OUTPUT_DIR"]); 
  s_file_radix = s_file_radix + "/" + my_map_paths["OUTPUT_RADIX"]; 
  VC_Movie.change_files(s_file_radix);
  VC_Movie.write_files();

     
  if (my_map_bools["save_sup_redon"])
    {
      string s_filename(my_map_paths["OUTPUT_DIR"]); 
      s_filename += my_map_paths["save_sup_redon_filename"].c_str();
      V_SupRedon.change_files(s_filename.c_str());
      V_SupRedon.write_hdr(); V_SupRedon.write_data();
    }

  vCHRONO_STOP(savedisk);
      
  // *******************************************************
  // fin
  
  delete centre;

  //libération memoire allouée pour les threads
  void fftw_cleanup_threads(void);
  
  vCHRONO_STOP(begin);

}

