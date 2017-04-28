#include <iostream>
#include <string>
#include <stdlib.h>

#include "macros.h"
#include "vChronos.h" // requiert boost

#include <boost/thread.hpp>  
#include <boost/date_time.hpp>

// -------------------------------
#include "main.h"
#include "util.h"
//#include "Holo_Process.h"
#include "cu_Holo_Process.h"

// -------------------------------
#include "cu_DisplayVolume.h"
#include "cpu_COMPLEXVolume.h"
#include "cu_ComplexVolume.h"
#include "cu_Volume.h"
#include "cuFFT_Volume.h"

// -------------------------------
using namespace std;

#include "cv.h"
#include "highgui.h"

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
using namespace cv;



//******************************************************************************
// 
//******************************************************************************


extern char *g_OUTPUT_DIR;
extern char *g_OUTPUT_RADIX;
extern size_t g_fftw_threads;
//extern bool g_save_sup_redon; 

extern size_t g_window_dim_x;
extern size_t g_window_dim_y;


//==============================================================================
// Thread de visualisation de volume
//==============================================================================


void 
BoostThreadFuncGPU(cu_DisplayVolume<RECON_TYPE> *slicer)
{
  char c;

  while(true)
    {
      slicer -> updateImage();
      slicer -> showImage();
      c = cvWaitKey(15);
    }
  
         
  cout << "Worker: finished" << endl;
}


// =============================================================================
//  traite les batches d'hologrammes sur CPU, envoie sup_redon et le volume des fréquences sur le GPU pour calculer la reconstruction courante
//
//  La visualisation s'effectue à partir du GPU
//
// =============================================================================


void
compute_gpu_online_batch(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
			 size_t window_edge_x, size_t window_edge_y, size_t image_dim_x, size_t image_dim_y, \
			 size_t xm0_limite, size_t ym0_limite, float rayon, float delta_zmax, \
			 size_t batch_size, \
			 const string &images_radix)
{

  // ==========================================================
  // Allocation des volumes de données
  // ==========================================================


  // CPU: volume des fréquences de fourier, rempli par la boucle de lecture
  cpu_COMPLEXVolume<RECON_TYPE> VC_Fourier_h(cube_edge, cube_edge, cube_edge);
  VC_Fourier_h.allocate();
  // pointeurs sur les sous-volumes
  AIR_Volume<RECON_TYPE>* V_FourierR_h = VC_Fourier_h.get_real_vol();
  AIR_Volume<RECON_TYPE>* V_FourierI_h = VC_Fourier_h.get_imag_vol();


  // GPU: le même
  cu_ComplexVolume<RECON_TYPE> VC_Fourier_d(cube_edge, cube_edge, cube_edge);
  VC_Fourier_d.allocate();
  VC_Fourier_d.fill(0);
  cu_Volume<RECON_TYPE>* V_FourierR_d = VC_Fourier_d.get_real_vol();
  cu_Volume<RECON_TYPE>* V_FourierI_d = VC_Fourier_d.get_imag_vol();


  // GPU: volume reconstruit à l'itération i, côté GPU
  cu_ComplexVolume<RECON_TYPE> VC_Visu_d(cube_edge, cube_edge, cube_edge);
  VC_Visu_d.allocate();
  VC_Visu_d.fill(0); // pour éviter le ghosting de la dernière acquisition
    // pointeurs sur les sous-volumes
  cu_Volume<RECON_TYPE>* V_VisuR_d = VC_Visu_d.get_real_vol();
  cu_Volume<RECON_TYPE>* V_VisuI_d = VC_Visu_d.get_imag_vol();



  // CPU: à la fin du calcul, volume qui va réceptionner la dernière reconstruction et la sauver sur disque
  cpu_COMPLEXVolume<RECON_TYPE> VC_Svg_h(cube_edge, cube_edge, cube_edge);
  VC_Svg_h.allocate();
  

  // GPU: on crée un volume calculable au format cuFFT
  cuFFT_Volume VF_Fourier_d(cube_edge, cube_edge, cube_edge);
  VF_Fourier_d.allocate();
  VF_Fourier_d.fill(0, 0);
  

  // CPU: on crée sup_redon pour normalisation dans fourier
  cpu_AIRVolume<unsigned short int> V_SupRedon_h(cube_edge, cube_edge, cube_edge);
  V_SupRedon_h.set_data_linear_mode(true);
  V_SupRedon_h.allocate();
  
  

  // GPU: le même
  cu_Volume<unsigned short int> V_SupRedon_d(cube_edge, cube_edge, cube_edge);
  V_SupRedon_d.allocate();
  V_SupRedon_d.fill(0);
  cu_Volume<unsigned short int> V_SupRedon2_d(cube_edge, cube_edge, cube_edge);
  V_SupRedon2_d.allocate();
  V_SupRedon2_d.fill(0);

  // ==========================================================
  // Début du Programme effectif
  // ==========================================================
    

  birdy();
  
  
  // ---------------------------------------------------
  // PRÉPARATION visu
  ASSERT(g_window_dim_x == g_window_dim_y);
  cu_DisplayVolume<RECON_TYPE> cu_VolSlicer("kilébo", V_VisuI_d);
  

  

  // *******************************************************
  // lecture des hologrammes sur CPU

  
  cu_Holo_Process Hologram_Processor;
  // Hologram_Processor.set_output_volumes(V_FourierR_h, V_FourierI_h, &V_SupRedon_h); 

  Hologram_Processor.set_output_volumes_gpu(V_FourierR_d -> get_data(), V_FourierI_d -> get_data(), V_SupRedon_d.get_data()); 
  
  Hologram_Processor.set_image_settings(image_dim_x, image_dim_y, window_edge_x, window_edge_y, g_window_dim_x);
  Hologram_Processor.set_rayon_settings(Nxmax, Nymax, Nxmax_Rf, xm0_limite, ym0_limite, rayon, delta_zmax);
  Hologram_Processor.set_file_settings(images_radix.c_str());
  Hologram_Processor.set_prepare();

  Hologram_Processor.display();



  size_t max_batch_error = 2; // was 2
  size_t errors = 0;
  size_t batch_num = 0;
  bool break_batches = false;
  char c;



  // ===========================================================================
  // LANCEMENT DU THREAD GUI en // 
  boost::thread workerThread(BoostThreadFuncGPU, &cu_VolSlicer);
  


  // ===========================================================================
  // BOUCLE DE TRAITEMENT DES BATCHES
  while(! break_batches)
    {

      // --------------------------------------------------
      // création d'un batch temps-réel (seuil max au-dela duquel on traite de force
      //, taille mini sous laquelle on attend, intervalle de dodo
      while (! Hologram_Processor.batch_prepare_realtime(400, 3, 50))
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
	  //vCHRONO_STOP(processing_time);
	  cerr << endl << "recon:: out with it";
	  
	  cout << endl << "recon:: waiting for keypress to exit";
	  cin.get(); //c = cvWaitKey(0); //
	  break;
	}

      
      // --------------------------------------------------
      // si la création a réussi

      batch_num++;
      size_t _batch_size = Hologram_Processor.batch_query_size();
      size_t _batch_start = Hologram_Processor.batch_query_angle_last_done() + 1;
      cout << endl << "batch n°" << batch_num << " [" <<  _batch_start\
	   << "; " << _batch_start - 1 + _batch_size << "]  ||" << _batch_size << "||";
      cout.flush();

      errors = 0;

      // batch launch!
      Hologram_Processor.batch_launch();

      
      // ==========================================================
      // Exécution de FFT 3D précédée et suivie de shifts circulaires (2.5s cpu)
      // ==========================================================


      // on a donné les volumes de VC_Fourier_d à Holo_Compute (gpu)
      //VF_Fourier_d.import_from( VC_Fourier_d );
      VF_Fourier_d.import_from_andfriend( VC_Fourier_d, V_SupRedon_d, V_SupRedon2_d );

      // et aussi V_SupRedon_d
      //V_SupRedon_d.circshift_to( V_SupRedon2_d );
      V_SupRedon2_d.nozero(1);
      
      VF_Fourier_d /= V_SupRedon2_d;

      VF_Fourier_d.set_fourier_backward();

      VC_Visu_d.fill(0);
      VF_Fourier_d.export_to( VC_Visu_d );

      cu_VolSlicer.notify_data_changed();      
      //cudaThreadSynchronize();  


    }
  // FIN BOUCLE BATCHES
  // ===========================================================================


  // *******************************************************
  // rapatriement du dernier volume reconstruit en GPU vers CPU
  VC_Visu_d.push_to_cpu( VC_Svg_h );
  

  // *******************************************************
  // sauvegardes sur disque

  bool g_save_sup_redon = false;
  
  

  string s_file_radix(g_OUTPUT_DIR); s_file_radix = s_file_radix +  "/" + g_OUTPUT_RADIX; 
  VC_Svg_h.change_files(s_file_radix);
  VC_Svg_h.write_files();
   
  
  if (g_save_sup_redon)
    {
      cerr << "unsupported";
      /*
      string s_filename(g_OUTPUT_DIR); s_filename += OUTPUT_REDON_FILENAME;
      V_SupRedon_h.change_files(s_filename.c_str());
      V_SupRedon_h.write_files();
      */
    }
      
  // *******************************************************
  // fin
  
  


}





      /*
      // ceci fonctionne très bien: Sup_Redon et Fréquences sur le CPU

      // send supredon to GPU
      V_SupRedon_d.pull_from_cpu( V_SupRedon_h );
      V_SupRedon_d.copy_to( V_SupRedon2_d );
      V_SupRedon2_d.nozero(1);

      // also send frequencies volume
      VC_Fourier_d.pull_from_cpu( VC_Fourier_h );
      

      // normalisation
      VC_Fourier_d /= V_SupRedon2_d;


      // and circshift and convert to cuFFT
      VF_Fourier_d.import_from( VC_Fourier_d );

      VF_Fourier_d.set_fourier_backward();
      
      // circshift et conversion inverse
      // mais copie sur le volume GPU de visu
      VF_Fourier_d.export_to( VC_Visu_d );
      cu_VolSlicer.notify_data_changed();

      */
