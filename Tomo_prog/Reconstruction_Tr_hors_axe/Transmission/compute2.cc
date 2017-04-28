void
compute_cpu_offline_batch(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
			  float rayon, float delta_zmax,		\
			  const string &images_radix, \
			  size_t batch_size, \
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
  

  // ==========================================================
  // Allocation des volumes de données
  // ==========================================================


  AIR_Volume<RECON_TYPE> *V_RealPart, *V_ImagPart;

#ifdef SSE4
  V_RealPart = new AIR_Volume_SSE4<RECON_TYPE> (cube_edge, cube_edge, cube_edge, 1.0);
  V_ImagPart = new AIR_Volume_SSE4<RECON_TYPE> (cube_edge, cube_edge, cube_edge, 1.0);
#else
  V_RealPart = new AIR_Volume<RECON_TYPE> (cube_edge, cube_edge, cube_edge, 1.0);
  V_ImagPart = new AIR_Volume<RECON_TYPE>(cube_edge, cube_edge, cube_edge, 1.0);
#endif


  V_RealPart -> allocate();
  V_ImagPart -> allocate();

  // binding to the 2 real and imag volumes
  COMPLEX_Volume<RECON_TYPE> VC_Volume(*V_RealPart, *V_ImagPart);


  // linéaire et 1.0 par défaut
  COMPLEX_Volume<RECON_TYPE> VC_View(cube_edge, cube_edge, cube_edge);


  // on crée un volume calculable au format fftw. 
  FFTW_Volume<RECON_TYPE> *VF_Volume;
#ifdef SSE4
  VF_Volume = new FFTW_Volume_SSE4<RECON_TYPE> (cube_edge, cube_edge, cube_edge);
#else
  VF_Volume = new FFTW_Volume<RECON_TYPE> (cube_edge, cube_edge, cube_edge);
#endif

  VF_Volume -> use_wisdom_directory("/usr/local/phd/fftw_wisdom");
  VF_Volume -> set_threads_initialized(true); // car déjà initialisé avant, echec sinon
  VF_Volume -> set_nb_threads(my_map_dims["fftw_threads"]);
  VF_Volume -> set_report_time(true); //déconne



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
  
  
  // *******************************************************
  // lecture des hologrammes 



  
  Holo_Process Hologram_Processor;
  //ASSERT(my_map_dims["window_dim_x"] == my_map_dims["window_dim_y"]);
  Hologram_Processor.set_path_settings(my_map_paths);
  Hologram_Processor.set_output_volumes(V_RealPart, V_ImagPart, V_SupRedon); 
  Hologram_Processor.set_image_settings(window_edge_x, window_edge_y, my_map_dims["window_dim_x"], my_map_bools["window_autocenter"]);
  Hologram_Processor.set_rayon_settings(Nxmax, Nymax, Nxmax_Rf, xm0_limite, ym0_limite, rayon, delta_zmax);
  Hologram_Processor.set_file_settings(images_radix.c_str(), my_map_bools["proper_numbers_mode"]);
  Hologram_Processor.set_prepare();

  Hologram_Processor.display();



  size_t nb_batchs = angle_last / batch_size;
  size_t angle_lastread = angle_first;


  for (size_t batch_num = 0; batch_num < nb_batchs; batch_num++)
    {
      cout << endl << endl << "batch n°" << batch_num + 1 << " || first angle: " << angle_lastread << "/ last angle: " << angle_lastread + batch_size;
      cout.flush();

      // on fixe une limite d'erreur basse pour passer vite en sommeil si fichiers pas là
      ASSERT( Hologram_Processor.batch_prepare(angle_lastread, batch_size, 4) );
      Hologram_Processor.batch_launch();
    
      
      // ==========================================================
      // Exécution de FFT 3D précédée et suivie de shifts circulaires
      // ==========================================================

      VF_Volume -> copy_from(VC_Volume); // copie et shift
      
      V_SupRedon -> circshift_to(*V_SupRedon_S);
      *VF_Volume /= *V_SupRedon_S;

      VF_Volume -> set_fourier_backward();
 
      VF_Volume -> copy_to(VC_View); // copie et shift
      
 
      angle_lastread += batch_size;

    }




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
