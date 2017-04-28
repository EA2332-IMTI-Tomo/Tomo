
// -------------------------------
#include "cu_DisplayVolume.h"
#include "cpu_COMPLEXVolume.h"
#include "cu_ComplexVolume.h"
#include "cu_Volume.h"
#include "cuFFT_Volume.h"

// -------------------------------



void
compute_cpu_online_batch(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
			 size_t window_edge_x, size_t window_edge_y, size_t image_dim_x, size_t image_dim_y, \
			 size_t xm0_limite, size_t ym0_limite, float rayon, float delta_zmax, \
			 size_t batch_size, \
			 const string &images_radix)
{
  cout << endl << "ATTENTION, version cuda bricolee";


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
  

  // linéaire et 1.0 par défaut
  COMPLEX_Volume<RECON_TYPE> VC_View(cube_edge, cube_edge, cube_edge);
  VC_View.allocate();


  // on crée un volume calculable au format fftw. 
  FFTW_Volume<RECON_TYPE> VF_Volume(cube_edge, cube_edge, cube_edge);
  VF_Volume.use_wisdom_directory("/usr/local/phd/fftw_wisdom");
  VF_Volume.set_threads_initialized(true); // car déjà initialisé avant, echec sinon
  VF_Volume.set_nb_threads(g_fftw_threads);
  VF_Volume.set_report_time(false); 



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
  centre = new int[4*Nxmax*Nymax];


    
  // ==========================================================
  // Début du Programme effectif
  // ==========================================================

    

  vCHRONO_SET(begin, "Lancement programme"); 
  vCHRONO_START(begin);


  birdy();
  
  // initialize multithreaded fftw 
  if (g_fftw_threads > 1)
    ASSERT(fftw_init_threads());
  
  
  // ---------------------------------------------------
  // PRÉPARATION visu
  ASSERT(g_window_dim_x == g_window_dim_y);

  AIR_Volume<RECON_TYPE> *VC_View_R = VC_View.get_real_vol();
  AIR_Volume<RECON_TYPE> *VC_View_I = VC_View.get_imag_vol();
  cvDisplayVolume<RECON_TYPE> cv_VolSlicer("kilébo", VC_View_I);
  

  

  // *******************************************************
  // lecture des hologrammes 

  
  Holo_Process Hologram_Processor;
  Hologram_Processor.set_output_volumes(&V_RealPart, &V_ImagPart, &V_SupRedon); 
  Hologram_Processor.set_image_settings(image_dim_x, image_dim_y, window_edge_x, window_edge_y, g_window_dim_x);
  Hologram_Processor.set_rayon_settings(Nxmax, Nymax, Nxmax_Rf, xm0_limite, ym0_limite, rayon, delta_zmax);
  Hologram_Processor.set_file_settings(images_radix.c_str());
  Hologram_Processor.set_prepare();

  Hologram_Processor.display();



  size_t max_batch_error = 2;
  size_t errors = 0;
  size_t batch_num = 0;
  bool break_batches = false;
  char c;



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
      cout << endl << endl << "batch n°" << batch_num << " [" <<  _batch_start\
	   << "; " << _batch_start - 1 + _batch_size << "]  ||" << _batch_size << "||";
      cout.flush();

      errors = 0;

      // batch launch!
      timer_loop.clear(); timer_loop.start();
      Hologram_Processor.batch_launch();
      timer_loop.stop();

      
      // ==========================================================
      // Exécution de FFT 3D précédée et suivie de shifts circulaires (2.5s cpu)
      // ==========================================================

      timer_recon.clear(); timer_recon.start();


      VF_Volume.copy_from(VC_Volume); // copie et shift
      
      V_SupRedon.circshift_to(V_SupRedon_S);
 
      VF_Volume /= V_SupRedon_S;

      VF_Volume.set_fourier_backward();
 
      VF_Volume.copy_to(VC_View); // copie et shift
 
      cv_VolSlicer.notify_data_changed();

      timer_recon.stop();
      
    }
  // FIN BOUCLE BATCHES
  // ===========================================================================


  // *******************************************************
  // sauvegardes sur disque

  bool g_save_sup_redon = false;


  vCHRONO_ASTART(savedisk, "sauvegarde sur disque des données");

  string s_file_radix(g_OUTPUT_DIR); s_file_radix = s_file_radix +  "/" + g_OUTPUT_RADIX; 
  VC_View.change_files(s_file_radix);
  VC_View.write_files();
   
  
  if (g_save_sup_redon)
    {
      string s_filename(g_OUTPUT_DIR); s_filename += OUTPUT_REDON_FILENAME;
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

