  
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

// on crée un volume calculable au format fftw. 
FFTW_Volume<RECON_TYPE> VF_Volume(cube_edge, cube_edge, cube_edge);
VF_Volume.use_wisdom_directory("/usr/local/phd/fftw_wisdom");
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
  

// *******************************************************
// lecture des hologrammes 



vCHRONO_SET(readloop, "Boucle de lecture fichiers"); 
vCHRONO_START(readloop);


   
readLoop_angleImages(&V_RealPart, &V_ImagPart, sup_redon, centre,	\
		     Nxmax, Nymax, Nxmax_Rf, coin_x, coin_y,		\
		     image_dimx, image_dimy,				\
		     xm0_limite, ym0_limite, rayon, delta_zmax,	\
		     premier_plan, NbAngle, SautAngle,		\
		     images_radix.c_str());
  

vCHRONO_STOP(readloop);
cout << endl;

     
    
// ==========================================================
// Préparation TF3D
// ==========================================================
  

// --------------------------------------------------
// hacker ponctuellement le volume prêt à reconstruire à fins de test
#ifdef SAVE_FFT_VOLUME
#include "hack_volumes.cc"
#endif
// --------------------------------------------------


      
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


  
  

  
// ===========================================================================
// ===========================================================================
// ===========================================================================
// conclusion partie commune




  
// *******************************************************
// sauvegardes sur disque


TIMES("Ecriture des résultats",
      {
	string s_file_real(g_OUTPUT_DIR); s_file_real = s_file_real +  "/" + g_OUTPUT_RADIX + "_R";
	cout << endl << s_file_real; cout.flush();
	V_RealPart.change_files(s_file_real.c_str());
	V_RealPart.write_hdr(); V_RealPart.write_data(); 

	string s_file_imag(g_OUTPUT_DIR); s_file_imag = s_file_imag + "/" + g_OUTPUT_RADIX + "_I";
	V_ImagPart.change_files(s_file_imag.c_str());
	V_ImagPart.write_hdr(); V_ImagPart.write_data(); 


	write_array<int>(size_t(4 * Nxmax * Nymax), centre, str_concat(g_OUTPUT_DIR, OUTPUT_CENTER_FILENAME));

	
	if (g_save_sup_redon)
	  {
	    string s_filename(g_OUTPUT_DIR); s_filename += OUTPUT_REDON_FILENAME;
	    V_SupRedon.change_files(s_filename.c_str());
	    V_SupRedon.write_hdr(); V_SupRedon.write_data();
	  }
      });
  
   
      
// *******************************************************
// fin
  
delete centre;

  
  
//libération memoire allouée pour les threads
void fftw_cleanup_threads(void);
  



vCHRONO_STOP(begin);

cout << endl;
