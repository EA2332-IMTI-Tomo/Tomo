
void
compute_cpu_zoom_offline(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
			 size_t window_edge_x, size_t window_edge_y, size_t image_dim_x, size_t image_dim_y, \
			 size_t xm0_limite, size_t ym0_limite, float rayon, float delta_zmax, \
			 size_t  angle_start, size_t angle_count, size_t angle_jump, \
			 bool off_axis_and_psh, size_t offaxis_cx, size_t offaxis_cy, \
			 bool autofocus_gradient,			\
			 const string &images_radix, bool proper_numbers, bool modref_asked, const char* modref_file)
{

  // ==========================================================
  // Allocation des volumes de données
  // ==========================================================


  AIR_Volume<RECON_TYPE> V_RealPart(cube_edge, cube_edge, cube_edge, 1.0);
  AIR_Volume<RECON_TYPE> V_ImagPart(cube_edge, cube_edge, cube_edge, 1.0);

  V_RealPart.set_data_linear_mode(true); 
  V_RealPart.allocate();
  V_ImagPart.set_data_linear_mode(true); 
  V_ImagPart.allocate();


  // binding to the 2 real and imag volumes
  COMPLEX_Volume<RECON_TYPE> VC_Volume(V_RealPart, V_ImagPart);
  
  
  // version de taille double pour padding
  COMPLEX_Volume<RECON_TYPE> VC_Volume_padded(2 * cube_edge, 2 * cube_edge, 2 * cube_edge);
  VC_Volume_padded.allocate();
  VC_Volume_padded.fill(0);

  
  // on crée un volume calculable au format fftw. 
  // il est de taille double
  FFTW_Volume<RECON_TYPE> VF_Volume(2 * cube_edge, 2 * cube_edge, 2 * cube_edge);
  // VF_Volume.use_wisdom_directory("/usr/local/phd/fftw_wisdom"); // pas encore de wisdom pour ca
  VF_Volume.set_threads_initialized(true); // car déjà initialisé avant, echec sinon
  VF_Volume.set_nb_threads(g_fftw_threads);
  VF_Volume.set_report_time(true); //déconne



  // pour normaliser par rapport au nombre de fois où l'on remplit la frequence
  unsigned short int *sup_redon;


  AIR_Volume<unsigned short int> V_SupRedon(cube_edge, cube_edge, cube_edge, 1.0);
  V_SupRedon.set_data_linear_mode(true);
  V_SupRedon.allocate();
  sup_redon = V_SupRedon.get_data_linear();


  // on va devoir conserver une copie circshiftée pour normaliser sur FFTW_Volume, à moins de modifier l'import en sauvage
  AIR_Volume<unsigned short int> V_SupRedon_S(cube_edge, cube_edge, cube_edge, 1.0);
  V_SupRedon_S.set_data_linear_mode(true);
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
  if (g_fftw_threads > 1)
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
  ASSERT(g_window_dim_x == g_window_dim_y);
  Hologram_Processor.set_output_volumes(&V_RealPart, &V_ImagPart, &V_SupRedon); 
  Hologram_Processor.set_image_settings(window_edge_x, window_edge_y, g_window_dim_x, my_map_bools["window_autocenter"]);
  Hologram_Processor.set_rayon_settings(Nxmax, Nymax, Nxmax_Rf, xm0_limite, ym0_limite, rayon, delta_zmax);
  Hologram_Processor.set_file_settings(images_radix.c_str(), proper_numbers, modref_asked, modref_file);

  // debug
  Hologram_Processor.d_record_centers(angle_count, "./centre_dist.txt");
  Hologram_Processor.d_record_fourier_spots("./svg_stack_fourier_spots");


  if ( off_axis_and_psh )
    {
      cout << endl << "OFF AXIS + PSH MODE accounted:";
      cout << off_axis_and_psh << " " << offaxis_cx << " " << offaxis_cy << endl; 
      Hologram_Processor.set_off_axis_and_psh( offaxis_cx, offaxis_cy );
    }
  if ( autofocus_gradient )
    {
      cout << endl << "AUTOFOCUS MODE accounted:";
      Hologram_Processor.set_autofocus_gradient( 20 );
      if (g_map_bools["autofocus_mode_debug"]) {
	cout << " + debug files on ";
	Hologram_Processor.focus_debug_mode("/ramdisk", angle_count);
      }
    }

  
  Hologram_Processor.set_prepare();
  Hologram_Processor.display();

  ASSERT( Hologram_Processor.batch_prepare(angle_start, angle_count, 100) );
  cout << endl << ">>> Boucle lecture hologrammes: "; cout.flush();
  Hologram_Processor.batch_launch();
  
  vCHRONO_STOP(readloop);

  //centre = Hologram_Processor.get_centres();
  //write_array<int>(size_t(4 * Nxmax * Nymax), centre, str_concat(g_OUTPUT_DIR, OUTPUT_CENTER_FILENAME));
  Hologram_Processor.set_free();

  // ==========================================================
  // Exécution de FFT 3D précédée et suivie de shifts circulaires
  // ==========================================================


  vCHRONO_ASTART(overall_fft, "whole FFT");


  VC_Volume /= V_SupRedon; // normalisation
  V_SupRedon.unallocate();

  VC_Volume_padded.copy_from_smaller(VC_Volume); // zero padding dans fourier
 

  VF_Volume.copy_from(VC_Volume_padded); // copie centrée et shift

  vCHRONO_ASTART(inner_fft, "inner fft");
  VF_Volume.set_fourier_backward();
  vCHRONO_STOP(inner_fft);
  VF_Volume.copy_to(VC_Volume_padded); // copie et shift

  VC_Volume.copy_from_bigger(VC_Volume_padded);


  vCHRONO_STOP(overall_fft);



  // *******************************************************
  // sauvegardes sur disque

  //bool g_save_sup_redon = false;

  vCHRONO_ASTART(savedisk, "sauvegarde sur disque des données");
	
  string s_file_radix(g_OUTPUT_DIR); s_file_radix = s_file_radix +  "/" + g_OUTPUT_RADIX; 
  VC_Volume.change_files(s_file_radix);
  VC_Volume.write_files();

	
  
  vCHRONO_STOP(savedisk);
      
  // *******************************************************
  // fin
  
  //delete centre;

  //libération memoire allouée pour les threads
  void fftw_cleanup_threads(void);
  
  vCHRONO_STOP(begin);


}
