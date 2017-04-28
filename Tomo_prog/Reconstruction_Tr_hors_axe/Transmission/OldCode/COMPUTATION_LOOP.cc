  
// ==========================================================
// Allocation des volumes de données
// ==========================================================

// volumes R et I où on accumule les projections non normalisées

AIR_Volume<RECON_TYPE> V_Projected_R(cube_edge, cube_edge, cube_edge, 1.0);
AIR_Volume<RECON_TYPE> V_Projected_I(cube_edge, cube_edge, cube_edge, 1.0);

V_Projected_R.set_data_linear_mode(true); 
V_Projected_R.allocate();
V_Projected_I.set_data_linear_mode(true); 
V_Projected_I.allocate();


// volumes R et I reconstruits

AIR_Volume<RECON_TYPE> V_Recon_R(cube_edge, cube_edge, cube_edge, 1.0);
AIR_Volume<RECON_TYPE> V_Recon_I(cube_edge, cube_edge, cube_edge, 1.0);

V_Recon_R.set_data_linear_mode(true); 
V_Recon_R.allocate();
V_Recon_I.set_data_linear_mode(true); 
V_Recon_I.allocate();


// pour compatibilité, pointeur sur les données en tableau 1D
RECON_TYPE* reel_arc = V_Projected_R.get_data_linear();
RECON_TYPE* imag_arc = V_Projected_I.get_data_linear();

// binding to the 2 real and imag volumes
COMPLEX_Volume<RECON_TYPE> VC_ProjectedVolume(V_Projected_R, V_Projected_I);
COMPLEX_Volume<RECON_TYPE> VC_ReconVolume(V_Recon_R, V_Recon_I);

// on crée un volume calculable au format fftw. 
FFTW_Volume<RECON_TYPE> VF_Volume(cube_edge, cube_edge, cube_edge);
VF_Volume.use_wisdom_directory("/usr/local/phd/fftw_wisdom");
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

    
TIME_START(begin, "Lancement programme\n");

birdy();
  
// initialize multithreaded fftw 
if (g_fftw_threads > 1)
  assert(fftw_init_threads());
  

// *******************************************************
// lecture des hologrammes 



vChrono<boost::chrono::system_clock> TIMER_read_loop;
cout << endl << "::Boucle de lecture START\n";
TIMER_read_loop.reset();

assert (!read_single_volume);
     
readLoop_angleImages(&V_Projected_R, &V_Projected_I, sup_redon, centre,		\
		     Nxmax, Nymax, Nxmax_Rf, coin_x, coin_y,	\
		     xm0_limite, ym0_limite, rayon, delta_zmax,	\
		     premier_plan, NbAngle, SautAngle,		\
		     Chemin);
  
cout << endl << "::Boucle de lecture STOP in " << TIMER_read_loop.seconds() << " (s)\n";

  

  

// --------------------------------------------------
// hacker ponctuellement le volume prêt à reconstruire à fins de test
#ifdef SAVE_FFT_VOLUME
#include "hack_volumes.cc"
#endif
// --------------------------------------------------


// *******************************************************
// écriture d'un "papillon" binaire à partir de reel_arc normalisé

/* a rétablir, on en a besoin pour convolution

   TIMES("papillon", {
   write_butterfly(reel_arc, N_tab, str_concat(g_OUTPUT_DIR, OUTPUT_BUTTERFLY_FILENAME));
   });

*/

      
// ==========================================================
// Exécution de FFT 3D précédée et suivie de shifts circulaires
// ==========================================================

  

// Si on choisit CUDA, le volume sera converti de double en float pour
// le calcul afin d'économiser de la place en GPU

#ifdef CUDA
fourier3DInverse_cshift(reel_arc, imag_arc, N_tab, cube_edge, g_fft3d_cuda);
#else

// version moderne, basée sur classe FFTW<T>. gère les circshifts

//fourier3DInverse_cshift(reel_arc, imag_arc, N_tab, cube_edge, false);
vChrono<boost::chrono::system_clock> TIMER_inner_fft;
vChrono<boost::chrono::system_clock> TIMER_overall_fft;
vChrono<boost::chrono::system_clock> TIMER_normalize;


TIMER_overall_fft.reset();
VF_Volume.copy_from(VC_ProjectedVolume);
TIMER_normalize.reset();
V_SupRedon.circshift_to(V_SupRedon_S);
VF_Volume /= V_SupRedon_S;
cout << endl << ">normalisation sup_redon: " << TIMER_normalize.seconds() << " (s)\n";
TIMER_inner_fft.reset();
VF_Volume.set_fourier_backward();
cout << endl << "::inner FFT STOP in " << TIMER_inner_fft.seconds() << " (s)\n";
VF_Volume.copy_to(VC_ReconVolume);
cout << endl << "::whole FFT STOP in " << TIMER_overall_fft.seconds() << " (s)\n";

#endif
  
// ca pointe toujours sur les volumes originaux, bindés par VC
RECON_TYPE *final_reel_shift = VC_ReconVolume.get_real_data();
RECON_TYPE *final_imag_shift = VC_ReconVolume.get_imag_data(); 
  

  
// ===========================================================================
// ===========================================================================
// ===========================================================================
// conclusion partie commune

  
// *******************************************************
// calcul du module carre (en double, ca ne dérange rien même si on était en float)

double *valeur_module_shift;

if (g_compute_module)
  {
    ARRAY_NEW(valeur_module_shift, N_tab, double);

    for(int cpt = 0, l_reel, l_imag; cpt < N_tab; cpt++)
      {
	l_reel = final_reel_shift[cpt];
	l_imag = final_imag_shift[cpt];
	valeur_module_shift[cpt] = l_reel * l_reel + l_imag * l_imag;
      }
  }
    
// *******************************************************
// sauvegardes sur disque


TIMES("Ecriture des résultats",
      {
	  
	// sizeof est compile-time, et precision est de type float en ce moment
	int data_size = sizeof(OUTPUT_TYPE);
	  
	// conversion et sauvegarde sur disque
	// (passage de 64 à 32 bits pour gain de place)
	convert_write_results(N_tab, str_concat(g_OUTPUT_DIR, OUTPUT_REAL_FILENAME), str_concat(g_OUTPUT_DIR, OUTPUT_IMAG_FILENAME), str_concat(g_OUTPUT_DIR, OUTPUT_MOD_FILENAME), final_reel_shift, final_imag_shift, valeur_module_shift, !(g_save_module && g_compute_module) );
	write_array<int>(size_t(4 * Nxmax * Nymax), centre, str_concat(g_OUTPUT_DIR, OUTPUT_CENTER_FILENAME));
	if (g_save_sup_redon)
	  write_array<SUPREDON_OUTPUT_TYPE>(N_tab, sup_redon, str_concat(g_OUTPUT_DIR, OUTPUT_REDON_FILENAME));

      });

  
  
// ===========================================================================
// ===========================================================================
// ===========================================================================
// essais de bertrand
  
      
// *******************************************************
// fin

  
// free(final_reel_shift);
// free(final_imag_shift);
if (g_save_module)
  free(valeur_module_shift);
delete centre;
free(sup_redon);
  
  
//libération memoire allouée pour les threads
void fftw_cleanup_threads(void);
  
TIME_END(begin, "Lancement programme\n");
