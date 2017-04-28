#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>



#include "cv.h"
#include "highgui.h"


using namespace std;
// using namespace Magick;

#include "vChrono.h"

#include "main.h"
#include "util.h"
// #include "util_Magick.h"
#include "util_Image.h"
#include "readloop.h"
#include "memory.h"
#include "vectra.h"

#include "fourier2D.h"
#ifdef CUDA 
#include "cuFourier.h"
#endif


#ifndef OMP_THREADS
#define OMP_THREADS 8
#endif




extern char *g_OUTPUT_DIR;
extern size_t g_fftw_threads;


// =============================================================================


inline double 
carre(double n)
{
  return n * n;
}


// =============================================================================


//double* reel_arc, double* imag_arc
void 
readLoop_angleImages(AIR_Volume<RECON_TYPE>* V_RealPart, AIR_Volume<RECON_TYPE>* V_ImagPart, \
		     unsigned short int* sup_redon, int* centre,	\
		     const int Nxmax, const int Nymax, int Nxmax_Rf,	\
		     const int window_edge_x, const int window_edge_y, const size_t image_dim_x, const size_t image_dim_y, \
		     const int xm0_limite, const int ym0_limite,	\
		     float rayon, float delta_zmax,			\
		     int angle_start, int angle_count, int angle_jump,	\
		     const char* images_radix)
{    
  // pour plus de facilité
  RECON_TYPE* reel_arc = V_RealPart -> get_data_linear();
  RECON_TYPE* imag_arc = V_ImagPart -> get_data_linear();

  //nombre de pixel dans l'image de destination
  const int N = g_window_dim_x * g_window_dim_y;
  

  //dimx_arc
  //const int dv0rf=4*Nxmax/Rf; //dimension totale apres correction par Rf
  const int dv0rf=4*Nxmax_Rf;
  //dimx_arc/2
  const int dv0s2rf=2*Nxmax_Rf;
  //const int dv0s2rf=2*(Nxmax)/Rf; ancienne valeur posant des problemes d'arrondi : N_tab#64*Nxmax^3 causant des débordment de tableau
  const int dv1s2rf=2*Nxmax_Rf;
  //const int dv1s2rf=2*(Nymax)/Rf;
  //dimx_arc*dimy_arc
  const int dv1xdv0rf=16*Nxmax_Rf*Nxmax_Rf;
  //const int dv1xdv0rf=16*Nxmax*Nymax/(Rf*Rf);
  //dimz_arc
  const int dv2s2rf=2*Nxmax_Rf;
  //const int dv2s2rf=2*Nxmax/Rf;



  // pour mesure cumulative du temps d'exécution
  clock_t temps_initial, temps_final;
  float temps_cpu = 0, temps_total = 0;

  // valeur maximale contenue dans sup_redon, tableau comptant le
  // nombre d'overlaps par voxel afin de permettre de normaliser avant TF
  size_t sup_redon_max;
  

  // ----------------------------------------------------------------------
  // allocations  
  // ----------------------------------------------------------------------

  IplImage *tampon_image_cv = cvCreateImage( cvSize(g_window_dim_x, g_window_dim_y), IPL_DEPTH_8U, 1); 
  

  unsigned char *masque, *cache_jumeau, *reference, *phase1, *phase2, *phase3, *phase4;
  double *phased1, *phased2, *phased3, *phased4;
  unsigned char *tampon_image;

  ARRAY_ALLOC(masque, N, unsigned char);     //masque hamming
  //reservation reference et "noir" camera pour une "tentative" de correction de la ref
  ARRAY_ALLOC(reference, N, unsigned char);
  
  /// 4 tableaux pour transtypage des images 8 bits vers 64 bits (double).
  ARRAY_ALLOC(phased1, N, double); ARRAY_ALLOC(phased2, N, double);
  ARRAY_ALLOC(phased3, N, double); ARRAY_ALLOC(phased4, N, double);
  
  //reservation de 4 tableaux (image 2D)  pour le phase shifting
  ARRAY_ALLOC(phase1, N, unsigned char); ARRAY_ALLOC(phase2, N, unsigned char);
  ARRAY_ALLOC(phase3, N, unsigned char); ARRAY_ALLOC(phase4, N, unsigned char);

  //reservation d'un tampon pour la découpe d'image, à la taille du fichier image
  ARRAY_ALLOC(tampon_image, (image_dim_x * image_dim_y), unsigned char);

    
  double *plan_reel, *plan_imag, *plan_reel_shift, *plan_imag_shift,	\
    *fft_reel_shift_norm, *fft_imag_shift_norm, *fft_reel_shift,	\
    *fft_imag_shift, *fft_module_shift, *fft_reel_tmp, *fft_imag_tmp, *fft_reel, *fft_imag;

  //int N=g_window_dim_x*g_window_dim_y;
  ARRAY_ALLOC(plan_reel, N, double);  ARRAY_ALLOC(plan_imag, N, double);
  ARRAY_ALLOC(plan_reel_shift, N, double);  ARRAY_ALLOC(plan_imag_shift, N, double);

  size_t NxNy4 = 4 * Nxmax * Nymax;
  ARRAY_ALLOC(fft_reel_shift, NxNy4, double);
  ARRAY_ALLOC(fft_imag_shift, NxNy4, double);
  ARRAY_ALLOC(fft_module_shift, NxNy4, double);
  ARRAY_ALLOC(fft_reel_shift_norm, NxNy4, double);
  ARRAY_ALLOC(fft_imag_shift_norm, NxNy4, double);

  
  //variables utilisées pour la TF2D
  
  //Avant Crop à NxMax;
  ARRAY_ALLOC(fft_reel_tmp, N, double); ARRAY_ALLOC(fft_imag_tmp, N, double);
  //Après Crop à NxMax;
  ARRAY_ALLOC(fft_reel, NxNy4, double); ARRAY_ALLOC(fft_imag, NxNy4, double);


  // ==================================================
  // init FFT2D

#ifndef CUDA 
  fftw_complex* temp_fftw1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N); //
  fftw_complex* temp_fftw2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N); //
#endif  
  

  // gestion des wisdoms fftw: cherche si disponible, sinon crée et sauve
  ARRAY_DEC_ALLOC(wisdom_filename, 256, char);
  compute_wisdom_filename_d(wisdom_filename, g_window_dim_x, g_window_dim_y, true, g_fftw_threads); // true = forw
  bool have_wisdom = importWisdom_f(wisdom_filename);
  if (! have_wisdom) 
    TF2D_compute_wisdom_d(g_window_dim_x, g_window_dim_y, true, g_fftw_threads);
  
  //
  // ==================================================


  //double *qualite_phase_shift=new double[N];
  //espace3D reel, imaginaire et support de redondance
  
  //tester l'existence des fichiers au cas ou
  char *holo_filename;
  assert(holo_filename = (char *) calloc(strlen(images_radix) + GENERATED_FILENAME_LEN, sizeof(char)));
  //char test_FinNom[15];


  // ----------------------------------------------------------------------
  // initialisations
  // ----------------------------------------------------------------------

  int nb_proj = 0;
  int points_faux = 0;
  int jumeau_elimine = 0;  
  const int rayon_inf = xm0_limite * xm0_limite + ym0_limite * ym0_limite;
  int centres_exclus = 0;


  // taille totale du masque en pixel pour elimination objet jumeau;
  const int t_mask = 30;
  char* chaine = str_concat(g_MASKPATH, "/k_30x30.pgm");
  ARRAY_ALLOC(cache_jumeau, 30 * 30, unsigned char);
  remplit_tableau_cv(cache_jumeau, chaine, t_mask, t_mask);

  //str_concat(g_MASKPATH, "/k_30x30.bmp"), 


  switch(g_window_dim_x)
    {
    case 512:   
      remplit_tableau_cv(masque, str_concat(g_MASKPATH, "/kmasq_512x512_30.pgm"), \
			 g_window_dim_x, g_window_dim_y);
      printf("masque 512x512\n");
      break;
    case 256:
      remplit_tableau_cv(masque, str_concat(g_MASKPATH, "/kmasq_256x256.pgm"), \
			 g_window_dim_x, g_window_dim_y);
      printf("masque 256x256\n");
      break;
    default:
      fprintf(stderr, "unsupported value for variable: g_window_dim_x");
      exit(EXIT_FAILURE);
    }


  /****************************************************************************/
  /****************************************************************************/
  // boucle principale de lecture
  /****************************************************************************/
  /****************************************************************************/


  // INIT --------------------------------------------------

  
  vChrono<boost::chrono::system_clock> TIMER_tf2d; double TIME_tf2d = 0;
  vChrono<boost::chrono::system_clock> TIMER_circ2d; double TIME_circ2d = 0;
  vChrono<boost::chrono::system_clock> TIMER_load2d; double TIME_load2d = 0;
  vChrono<boost::chrono::system_clock> TIMER_sonde; double TIME_sonde = 0;

  // pre-calcul du tableau des images à traiter
  int* ANGLES_TO_PROCESS;
  int ANGLES_size_alloc = 1 + (angle_count / angle_jump);
  int ANGLES_cardinal = 0;
  ARRAY_ALLOC(ANGLES_TO_PROCESS, ANGLES_size_alloc, int);
  for (int cpt_angle = angle_start; cpt_angle < angle_start + angle_count; cpt_angle = cpt_angle + angle_jump)
    {

      compute_hologram_filename(images_radix, holo_filename, cpt_angle, 1);
      if (vectra::file_exists_p(holo_filename))
	{
	  ANGLES_TO_PROCESS[ANGLES_cardinal] = cpt_angle;
	  ANGLES_cardinal++;
	}
      else
	cerr << endl << "not found: " << holo_filename;
      
    }
  cout << endl << "Nombre de fichiers à traiter: " << ANGLES_cardinal;
  if (! ANGLES_cardinal) {exit(0);}

  
  ///SAVE_INPUT_VOLUME
  FILE *output_vol;
  if (SAVE_INPUT_VOLUME) 
    assert(output_vol = fopen(str_concat(g_OUTPUT_DIR, "/volume_realImg.raw"), "wb"));
  OUTPUT_TYPE* image_out;
  if (SAVE_INPUT_VOLUME) 
    ARRAY_ALLOC(image_out, N, OUTPUT_TYPE);
  ///SAVE_INPUT_VOLUME
  

  // BODY --------------------------------------------------
  for (int I = 0; I < ANGLES_cardinal; I++)	 
    {
      int cpt_angle = ANGLES_TO_PROCESS[I];

      if((cpt_angle - 100 * (cpt_angle / 100)) == 0)
	printf("\ncpt_angle=%i",cpt_angle);

      //----------------------------------------------------------------------
      //////////////////////////phase shifting////////////////////////////////
      //----------------------------------------------------------------------

      //**  On teste jusqu'à  angle_end fichiers
      //    On peut faire la reconstruction avec relativement peu de données,
      //    donc l'absence de certaines séries de coupes n'a aucune incidence
      //    grave et ne doit pas interrompre l'analyse
          

      
      TIMER_load2d.reset();  // 0.66s cumulé

      //       sprintf(holo_filename, "%s%d-%03d." INPUT_FORMAT, images_radix, cpt_angle, 1);

      compute_hologram_filename(images_radix, holo_filename, cpt_angle, 1);
      charge_decoupe_image_cv(phase1, holo_filename, tampon_image_cv, window_edge_x, window_edge_y, g_window_dim_x, g_window_dim_y);
      compute_hologram_filename(images_radix, holo_filename, cpt_angle, 2);
      charge_decoupe_image_cv(phase2, holo_filename, tampon_image_cv, window_edge_x, window_edge_y, g_window_dim_x, g_window_dim_y);
      compute_hologram_filename(images_radix, holo_filename, cpt_angle, 3);
      charge_decoupe_image_cv(phase3, holo_filename, tampon_image_cv, window_edge_x, window_edge_y, g_window_dim_x, g_window_dim_y);
      compute_hologram_filename(images_radix, holo_filename, cpt_angle, 4);
      charge_decoupe_image_cv(phase4, holo_filename, tampon_image_cv, window_edge_x, window_edge_y, g_window_dim_x, g_window_dim_y);

      TIME_load2d += TIMER_load2d.seconds();
      


      //** remise à zéro (non nécessaire)
      /*
	zerofill_double(plan_reel, N);
	zerofill_double(plan_imag, N);
	zerofill_double(plan_reel_shift, N);
	zerofill_double(plan_imag_shift, N);
      */

      
#pragma omp parallel for num_threads( OMP_THREADS ) 
      for(int pixel = 0; pixel < N; pixel++)
	{
	  plan_reel[pixel] = (double)((phase1[pixel] - phase3[pixel]) * (masque[pixel]));
	  plan_imag[pixel] = (double)((phase4[pixel] - phase2[pixel]) * (masque[pixel]));
	}

      if (SAVE_INPUT_VOLUME)
	{
	  convert_float(plan_reel, image_out, N);
	  fwrite(image_out, sizeof(OUTPUT_TYPE), N, output_vol);
	  convert_float(plan_imag, image_out, N);
	  fwrite(image_out, sizeof(OUTPUT_TYPE), N, output_vol);
	}

      
      //**  */
      TIMER_circ2d.reset(); //0.4s cumulé
      circshift_2D_mc(plan_reel, plan_reel_shift, g_window_dim_x, g_window_dim_y); 
      circshift_2D_mc(plan_imag, plan_imag_shift, g_window_dim_x, g_window_dim_y); 	
      //       circshift_ccd2(plan_reel, plan_reel_shift, g_window_dim_x, g_window_dim_y, g_window_dim_x/2, g_window_dim_y/2);
      //       circshift_ccd2(plan_imag, plan_imag_shift, g_window_dim_x, g_window_dim_y, g_window_dim_x/2, g_window_dim_y/2);
      TIME_circ2d += TIMER_circ2d.seconds();


      // ----------------------------------------------------------------------
      // TF 2D du front d'onde: passage du plan image au plan réciproque///////
      // ----------------------------------------------------------------------

      
      TIMER_tf2d.reset(); //2.5s cumulé en CPU

#ifdef CUDA
      TF2D_cuda(plan_reel_shift, plan_imag_shift, fft_reel_tmp, fft_imag_tmp, \
		g_window_dim_x, g_window_dim_y);
#else
      TF2D_d(plan_reel_shift, plan_imag_shift, fft_reel_tmp, fft_imag_tmp, \
	     g_window_dim_x, g_window_dim_y, temp_fftw1, temp_fftw2, wisdom_filename, g_fftw_threads);
#endif
      
      TIME_tf2d += TIMER_tf2d.seconds();


      // **      delete[] plan_reel_shift;
      // **      delete[] plan_imag_shift;
      // ----------------------------------------------------------------------
      // Découpage dans fourier à Nxmax (crop dans fourier = 
      // mise à l'échelle dans plan image)
      // la frequence max correspond maintenant a 
      // l'ouverture numerique experimentale
      // ----------------------------------------------------------------------

      ///////////////fft_reel_tmp=entree, fft_reel=sortie//////////////////////
      //-----------------------------------------------------------------------

      if (!I)
	cout << endl << "valeurs: " << " | Nxmax:" << Nxmax << " | Nymax:" << Nymax << " | g_window_dim_x:" << g_window_dim_x << " | g_window_dim_y:" << g_window_dim_y << endl;
      



      // deux espaces: temps cumulé de 0.35s -> 0.16 avec OMPx3 / 0.13 x4 / augm. si x5
      size_t xi, yi;
      const size_t cg_window_dim_x = g_window_dim_x; const size_t cg_window_dim_y = g_window_dim_y;
      // premier demi-espace /////////////////////////////////////////////
#pragma omp parallel for num_threads( 4 ) 
      for (xi = 0; xi < Nxmax; xi++)
	{

	  for ( yi = 0; yi < Nymax; yi++)
	    {
	      int cpt1 = yi * cg_window_dim_x + xi;
	      int cpt2 = yi * Nxmax * 2 + xi;

	      fft_reel[cpt2] = fft_reel_tmp[cpt1];
	      fft_imag[cpt2] = fft_imag_tmp[cpt1];
	    }


	  for (int yi = g_window_dim_y - Nymax; yi < g_window_dim_y; yi++)
	    {
	      int cpt1 = yi * cg_window_dim_x + xi;
	      int cpt2 = xi + ( yi - cg_window_dim_y + 2 * Nymax) * 2 * Nxmax;

	      fft_reel[cpt2] = fft_reel_tmp[cpt1];
	      fft_imag[cpt2] = fft_imag_tmp[cpt1];
	    }
	}
#pragma omp barrier

      // deuxieme demi-espace /////////////////////////////////////////////
#pragma omp parallel for num_threads( 4 ) 
      for(int xi = cg_window_dim_x - Nxmax ;xi < cg_window_dim_x; xi++)
	{
	  for (int yi = 0; yi < Nymax; yi++)
	    {
	      int cpt1 = yi * cg_window_dim_x + xi;
	      int cpt2 = yi * Nxmax * 2 + (xi - cg_window_dim_x + 2 * Nxmax);

	      fft_reel[cpt2] = fft_reel_tmp[cpt1];
	      fft_imag[cpt2] = fft_imag_tmp[cpt1];
	    }

	  for (int yi = cg_window_dim_y - Nymax; yi < cg_window_dim_y; yi++)
	    {
	      int cpt1 = yi * cg_window_dim_x + xi;

	      int cpt2 = (xi - cg_window_dim_x + 2 * Nxmax) + (yi - cg_window_dim_y + 2 * Nymax)* 2 *(Nxmax);
	      fft_reel[cpt2] = fft_reel_tmp[cpt1];
	      fft_imag[cpt2] = fft_imag_tmp[cpt1];
	    }
	}
     

      TIMER_circ2d.reset();
      circshift_2D_mc(fft_reel, fft_reel_shift, 2 * Nxmax, 2 * Nymax); 
      circshift_2D_mc(fft_imag, fft_imag_shift, 2 * Nxmax, 2 * Nymax); 
      TIME_circ2d += TIMER_circ2d.seconds();

      //----------------Recherche du maximum non centré------------------------


      //gloubiboulga function /!\ achtung, dont touch /!!\ Caramba, Né tochas
      
      // seul cpt_max est réutilisé;
      size_t cpt_max=0; double max_val = 0;
      
      fft_module_shift[0] = carre( fft_reel_shift[0] ) + carre( fft_imag_shift[0] );
  
       
      //Recherche du MAX dansle module //0.14s cumulé -> 0.08 avec dépaysement max /0.02 avec inclusion res
      const size_t NxNyMax4 = 4 * Nxmax * Nymax;
      
#pragma omp parallel for num_threads( 3 ) 
      for(size_t cpt = 1; cpt < NxNyMax4; cpt++)
	{
	  double res = fft_module_shift[cpt] = carre(fft_reel_shift[cpt]) + carre(fft_imag_shift[cpt]); 
	  if (res > max_val)
	    {
	      cpt_max = cpt;
	      max_val = fft_module_shift[cpt_max];
	    }
	}
      //pow(fft_reel_shift[cpt],2)+pow(fft_imag_shift[cpt],2);
      // 	  if(fft_module_shift[cpt]>fft_module_shift[cpt_max])

      /////////////////////////////
      //on met 1 dans l'image centre en cpt_max
      /////////////////////////////

      //centre[cpt_max]++;


      double max_part_reel = fft_reel_shift[cpt_max];
      double max_part_imag = fft_imag_shift[cpt_max];
      double max_module = fft_module_shift[cpt_max];
      double max_mod_divisor = 1 / max_module;

      // ----------------------------------------------------------------------
      ////////////////////////////normalisation par le pic central
      // ----------------------------------------------------------------------

      


      // 0.17s -> 0.07 avec *divisor // 0.04 omp
#pragma omp parallel for num_threads( 4 )
      for(size_t cpt = 0; cpt < NxNyMax4; cpt++)
	{
	  fft_reel_shift_norm[cpt] = (fft_reel_shift[cpt] * max_part_reel + fft_imag_shift[cpt] * max_part_imag) * max_mod_divisor;
	  fft_imag_shift_norm[cpt] = (fft_imag_shift[cpt] * max_part_reel - fft_reel_shift[cpt] * max_part_imag) * max_mod_divisor;
	}

      //fft_reel[cpt]=(fft_reel[cpt]*fft_reel[cpt_max]+fft_imag[cpt]*fft_imag[cpt_max])/fft_module[cpt_max];danger!!!!
      //fft_imag[cpt]=(fft_imag[cpt]*fft_reel[cpt_max]-fft_reel[cpt]*fft_imag[cpt_max])/fft_module[cpt_max];danger!!!!


      //--------------------------------------


      //Coordonnées dans l'espace 2D à partir de l'indice 1D: (xc,yc)=(0,0)=en haut à gauche de l'image
      int xmi = cpt_max % (2 * Nxmax);
      int ymi = cpt_max / (2 * Nymax);

      int xc = xmi - Nxmax;
      int yc = ymi - Nymax;
      //coordonnée objet jumeaux
      int xmij = 2 * Nxmax-xmi;
      int ymij = 2 * Nymax-ymi;
      
      
      /* //DEBUG 
      cout << endl << ": xmi ymi xc yc xmij ymij";
      cout << xmi << " " << ymi << " " << xc << " " << yc << " " << xmij << " " << ymij;
      */ //DEBUG 
      

      ////////////////////////////////////Virer l'objet jumeaux (0.005 s cumulé)
      if((xc*xc+yc*yc) > 900)//35*35=1225;objet jumeau pas trop près de l'objet
	{ 
	  jumeau_elimine++;

	  for(int x =xmij-t_mask/2;x<xmij+t_mask/2;x++)
	    {
	      for(int y=ymij-t_mask/2;y<ymij+t_mask/2;y++)
		{

		  //attention le masque repasse du coté gauche lorsque le jumeau touche le bord droit! A corriger
		  int x_2=x;
		  int y_2=y;
		  if(x>2*Nxmax)
		    x_2=x-2*Nxmax;
		  if(y>2*Nymax)
		    y_2=y-2*Nymax;
		  if(x<0)
		    x_2=2*Nxmax+x;

		  if(y<0)
		    y_2=2*Nymax+y;

		  int cptj=2*Nymax*y_2+x_2;
		  fft_reel_shift_norm[cptj]=fft_reel_shift_norm[cptj]*cache_jumeau[t_mask*(y-ymij+t_mask/2)+(x-xmij+t_mask/2)]/255;
		  fft_imag_shift_norm[cptj]=fft_imag_shift_norm[cptj]*cache_jumeau[t_mask*(y-ymij+t_mask/2)+(x-xmij+t_mask/2)]/255;
		}//fin for y
	    }//fin for x
	}//fin if objet jumeau pas trop pres
	    

      /* //DEBUG
      Array<double> a_plan_reel(fft_reel_shift_norm, N);
      Array<double> a_plan_imag(fft_imag_shift_norm, N);
      cout << endl<< "2d fft reel imag::" << a_plan_reel.sum() << " " << a_plan_imag.sum();
      */ //DEBUG


      //----------------------------------------------------------------------
      ////////////////////////////Mapping 3D
      //----------------------------------------------------------------------


      //coordonnée dans l'image2D centrée (xm0,ym0)=(0,0)=au centre de l'image
      int xm0=(xmi-Nxmax);//
      //printf("xm0 %i \n", xm0);
      int ym0=(ymi-Nymax);

      if(xm0==0 && ym0==0)
	printf("(xm0,ym0)=(0,0) dans le plan %i\n", cpt_angle);

      
      TIMER_sonde.reset(); //1.87s cumulé!

      //cout << endl << xm0 << ": " << ym0 << " " << xm0*xm0 + ym0*ym0 << " " <<  rayon_inf;
      if((xm0*xm0+ym0*ym0) > rayon_inf)
	{
	  //printf("xm0 : %i, ym0:%i\n",xm0,ym0);
	  //printf("\nFichier %i sans signal\n", cpt_angle);
	  centres_exclus++; MSG_ASSERT(false, "tiens, jamais arrivé encore");
	}
      else
	{
	  //sauvegarde des centres en coordonnées non centrée; on met le numero d'angle
	  centre[xmi*2*Nxmax+ymi]=cpt_angle;
	  //pour vérifier les pb. A ouvrir en sizeof(int)=32 bits mais non signé
	  //printf("cpt_angle:%i, xmi: %i, ymi : %i\n",cpt_angle,xmi,ymi);

	  //création de variable pour éviter N calculs dans la boucle sur le volume 3D

	  float r2=rayon*rayon;

	  //indice 3D
	  int k=0;
	  double arg_z_arc=0;
	  double z_arc=0;
	  double zm0;

	  //printf("round(rayon*rayon-(xm0)^2-(ym0)^2: %i\n",round(rayon*rayon-(xm0)^2-(ym0)^2));

	  double zm0_carre = rayon * rayon - xm0 * xm0 - ym0 * ym0;
	  if(round(zm0_carre) > -1)
	    {

	      zm0=sqrt(zm0_carre);
	      //printf("zm0:%i\n",zm0);
	      nb_proj++;

	     
	      temps_initial = clock ();
	      //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
	      for (int x = -Nxmax; x < Nxmax; x++)
		{
		  //on balaye l'image 2D en y, centre au milieu
		  for (int y = -Nymax; y < Nymax; y++)
		    {
		      //calcul du cpt du tableau 1D de l'image 2D
		      int cpt = (y + Nymax)* 2 * Nxmax + x + Nxmax;
		      //argument sous la racine calcul dans ARC_centre
		      arg_z_arc = r2 - x * x - y * y;

		      //ne pas depasser l'ouverture numérique pour 1 hologramme
		      if(arg_z_arc > delta_zmax * delta_zmax)
			{
			  //altitude au centre du volume
			  z_arc = round(sqrt(arg_z_arc) - zm0);
			  //indice du tableau 1D du volume 3D
			  k = round(( - xm0 +x + dv0s2rf) + (- ym0 + y + dv1s2rf) * dv0rf + \
				    (z_arc + dv2s2rf) * dv1xdv0rf);

			  // sup_redon est un tableau de même taille que le volume
			  sup_redon[k] += 1;//pour calculer le support
			  //
			  if (sup_redon[k] > sup_redon_max) sup_redon_max = sup_redon[k];
			    

			  reel_arc[k] += fft_reel_shift_norm[cpt];//pour calculer l'image
			  imag_arc[k] += fft_imag_shift_norm[cpt];//pour calculer l'image

			}
		      else
			points_faux++;
		    } //fin for y
		}

	    }//fin if zm0>-1
	}//fin else xm0_limite

      TIME_sonde += TIMER_sonde.seconds();

      temps_final = clock ();
      temps_cpu = (temps_final - temps_initial) / CLOCKS_PER_SEC;
      temps_total = temps_total + temps_cpu;
      temps_initial = clock();
      
    } 

  // fin boucle principale de lecture sur tous les angles
  // END OF LOOP
  
  cout << endl << "Temps passé à invoquer la TF2D: " << TIME_tf2d;
  cout << endl << "Temps passé à circshifter en 2D: " << TIME_circ2d;  
  cout << endl << "Temps passé à lire les images 2D: " <<  TIME_load2d;
  cout << endl << "Temps passé à sonde: : " <<  TIME_sonde << endl;


  /* //DEBUG
  cout << endl << "contrôle reel_arc/imag_arc";
  size_t cub = 512 * 512 * 512;
  cout.precision(10);
  Array<RECON_TYPE> a_plan_reel(reel_arc, cub);
  Array<RECON_TYPE> a_plan_imag(imag_arc, cub);
  cout << endl<< "(class)reelarc fin::" << a_plan_reel.sum() << " " << a_plan_imag.sum();
  */ //DEBUG


  if (SAVE_INPUT_VOLUME)
    fclose(output_vol);

  // ----------------------------------------------------------------------
  // reporting, mesure du temps de calcul
  // ----------------------------------------------------------------------

  
  printf("jumeaux re-elimines : %i\n",jumeau_elimine);
  printf("nb_proj: %i\n", nb_proj);
  printf("Nombre de centre exclus : %i\n", centres_exclus);  
  printf("points_faux : %i,points_faux\n",points_faux);
  printf("temps_total proj : %f \n",temps_total);
  printf("coucou \n");

  
  temps_final = clock ();
  temps_cpu = (temps_final - temps_initial) / CLOCKS_PER_SEC / nb_proj;
  printf("temps moyen pour 1 angle: %f\n", temps_cpu);
  temps_cpu = (temps_final - temps_initial) / CLOCKS_PER_SEC;
  printf("temps total pour %i angle(s): %f\n", nb_proj, temps_cpu);


  cout << endl << "sup_redon_max: " << sup_redon_max;

  
  // ----------------------------------------------------------------------
  // desallocations
  // ----------------------------------------------------------------------


  free(masque); //free(cache_jumeau; //
  free(cache_jumeau);

  free(phased1); free(phased2); free(phased3); free(phased4); 
  free(phase1); free(phase2); free(phase3); free(phase4); 

  free(tampon_image);
  free(reference); 

  free(fft_reel); free(fft_imag); 
  free(fft_reel_tmp); free(fft_imag_tmp); 
  free(fft_reel_shift); free(fft_imag_shift); free(fft_module_shift); 
  free(fft_reel_shift_norm); free(fft_imag_shift_norm); 
  free(plan_reel); free(plan_imag); 
  free(plan_reel_shift); free(plan_imag_shift); 


#ifndef CUDA
  fftw_free( temp_fftw1 );//
  fftw_free( temp_fftw2 ); //
#endif

}



/*
  try
  {
  //masque hamming
  masque = new unsigned char[N]; //
  //cache objet jumeau
  }
  catch(bad_alloc ex)
  {cerr << endl << "memory allocation failure" << ex.what(); exit(EXIT_FAILURE); } 

*/


