void
compute_cpu_online_batch(std::map<string, float> my_map_experiment,	\
			 const string &images_radix, unsigned char **images_list, \
			 std::map<string, bool> my_map_bools,  \
			 std::map<string, size_t> my_map_dims,	\
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
  
  size_t cube_edge = (size_t) my_map_experiment["cube_edge"];
  size_t Nxmax = (size_t) my_map_experiment["nxmax"];
  size_t Nymax = Nxmax;
  size_t Nxmax_Rf = (size_t) my_map_experiment["nxmax_rf"];
  float rayon = my_map_experiment["rayon"];
  float delta_zmax = my_map_experiment["delta_zmax"];

  float pixel_size_recon = my_map_experiment["pixel_size_recon"];

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
  

  // linéaire et 1.0 par défaut
  //*********** aie aie aie
  COMPLEX_Volume<RECON_TYPE> VC_View(cube_edge, cube_edge, cube_edge);
  VC_View.allocate();
  

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
  VF_Volume -> set_report_time(false);
  


  // pour normaliser par rapport au nombre de fois où l'on remplit la frequence
  unsigned short int *sup_redon;


  AIR_Volume<unsigned short int> *V_SupRedon;
#ifdef SSE4
  V_SupRedon = new AIR_Volume_SSE4<unsigned short int> (cube_edge, cube_edge, cube_edge, 1.0);
#else
  V_SupRedon = new AIR_Volume<unsigned short int> (cube_edge, cube_edge, cube_edge, 1.0);
#endif


  V_SupRedon -> allocate();
  sup_redon = V_SupRedon -> get_data();


  // on va devoir conserver une copie circshiftée pour normaliser sur FFTW_Volume, à moins de modifier l'import en sauvage
  AIR_Volume<unsigned short int> *V_SupRedon_S;
#ifdef SSE4
  V_SupRedon_S = new AIR_Volume_SSE4<unsigned short int> (cube_edge, cube_edge, cube_edge, 1.0);
#else    
  V_SupRedon_S = new AIR_Volume<unsigned short int> (cube_edge, cube_edge, cube_edge, 1.0);
#endif


  V_SupRedon_S -> allocate();


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
  cvDisplayVolume<RECON_TYPE> cv_VolSlicer("Reconstruction - Partie Imaginaire", VC_View_I);
  cv_VolSlicer.set_cv_tempo(30);
  cv_VolSlicer.enable_mouse();
  cv_VolSlicer.enable_histogram(256, 100, 256);
  cv_VolSlicer.enable_scale(1000000000 * pixel_size_recon, 10); // nm, um
  

  // *******************************************************
  // lecture des hologrammes 

  
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
  Hologram_Processor.display();



  size_t max_batch_error = 2;
  size_t errors = 0;
  size_t batch_num = 0;
  bool break_batches = false;
  char* str_angles_done; 
  ARRAY_ALLOC(str_angles_done, 512, char);
  unsigned char c;


  // ===========================================================================
  // LANCEMENT DU THREAD GUI en // 
  boost::thread workerThread(BoostThreadFunc, &cv_VolSlicer);
  
  

  vCHRONO_SET(processing_time, "effective processing time besides disk save and wait");
  bool chrono_started = false;

  vChronos timer_loop("batch processing loop");
  vChronos timer_recon("batch to volume processing");
  vChronos timer_cv("prepare & show current layer");
  

  // ===========================================================================
  // BOUCLE DE TRAITEMENT DES BATCHES
  while(! break_batches)
    {


      // --------------------------------------------------
      // création d'un batch temps-réel
      while (! Hologram_Processor.batch_prepare_realtime(200, 15))
	{
	  errors++;
	  if (errors > max_batch_error)
	    {
	      cerr << endl << "recon:: probably reached end of acquisition, halting"; 
	      cerr.flush();
	      break_batches = true;
	      break;
	    }
	}
      
      // only exit point possible
      if (break_batches) 
	{
	  vCHRONO_STOP(processing_time);
	  cerr << endl << "recon:: out with it";
	  
	  cout << endl << endl << "end of image reconstruction: ";
	  time_t now; time(&now); printf("%s\n", ctime(&now));

	  cout << endl << "recon:: waiting for keypress to exit";
	  c = cvWaitKey(0); //cin.get();
	  break;
	}

      
      // --------------------------------------------------
      // si la création a réussi

      if (! chrono_started)
	{
	  vCHRONO_START(processing_time);
	  chrono_started = true;
	}

      batch_num++;
      size_t _batch_size = Hologram_Processor.batch_query_size();
      size_t _batch_start = Hologram_Processor.batch_query_angle_last_done() + 1;
      size_t _batch_end = _batch_start - 1 + _batch_size;
      cout << endl << endl << "batch n°" << batch_num << " [" <<  _batch_start\
	   << "; " << _batch_end << "]  ||" << _batch_size << "||";
      sprintf(str_angles_done, "angles: %d", _batch_end);
      
      cout.flush();

      errors = 0;

      // batch launch!
      timer_loop.clear(); timer_loop.start();
      Hologram_Processor.batch_launch();
      timer_loop.stop();

      
      // ==========================================================
      // Exécution de FFT 3D précédée et suivie de shifts circulaires (2.5s cpu)
      // ==========================================================

      //MSG_ASSERT(false, "stop, la");

      timer_recon.clear(); timer_recon.start();
 
      //VC_Volume /= V_SupRedon; non, destructif! on est pas en GPU

      VF_Volume -> copy_from(VC_Volume); // copie et shift
      
      V_SupRedon -> circshift_to(*V_SupRedon_S);
      *VF_Volume /= *V_SupRedon_S;

      VF_Volume -> set_fourier_backward();
 
      VF_Volume -> copy_to(VC_View); // copie et shift
 
      cv_VolSlicer.notify_data_changed();
      cv_VolSlicer.set_Label(str_angles_done);

      timer_recon.stop();
      
    }
  // FIN BOUCLE BATCHES
  // ===========================================================================


  // *******************************************************
  // sauvegardes sur disque


  vCHRONO_ASTART(savedisk, "sauvegarde sur disque des données");

  string s_file_radix(my_map_paths["OUTPUT_DIR"]); 
  s_file_radix = s_file_radix + "/" + my_map_paths["OUTPUT_RADIX"]; 
  VC_View.change_files(s_file_radix);
  VC_View.write_files();
   
  
  if (my_map_bools["save_sup_redon"])
    {
      string s_filename(my_map_paths["OUTPUT_DIR"]); 
      s_filename += my_map_paths["save_sup_redon_filename"];
      V_SupRedon -> change_files(s_filename.c_str());
      V_SupRedon -> write_hdr(); 
      V_SupRedon -> write_data();
    }

  vCHRONO_STOP(savedisk);
      
  // *******************************************************
  // fin
  
  delete centre;

  //libération memoire allouée pour les threads
  void fftw_cleanup_threads(void);
  
  vCHRONO_STOP(begin);

}

