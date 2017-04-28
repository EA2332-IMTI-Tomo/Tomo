#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>
#include <Magick++.h>


using namespace std;
using namespace Magick;


#include "main.h"
#include "util.h"
#include "readloop.h"
#include "memory.h"
#include "fourier.h"
#ifdef CUDA 
#include "cuFourier.h"
#endif


void 
readLoop_angleImages(double* reel_arc, double* imag_arc, unsigned short int* sup_redon, int* centre, \
	 const int Nxmax, const int Nymax, int Nxmax_Rf, \
	 const int window_edge_x, const int window_edge_y,				\
	 const int xm0_limite, const int ym0_limite,			\
	 float rayon, float delta_zmax,					\
	 int angle_start, int angle_count, int angle_jump,		\
	 char* Chemin)
{    

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
  

  // ----------------------------------------------------------------------
  // allocations  
  // ----------------------------------------------------------------------


  unsigned char *masque, *cache_jumeau, *reference, *phase1, *phase2, *phase3, *phase4;
  double *phased1, *phased2, *phased3, *phased4;
  try
    {
      //masque hamming
      masque = new unsigned char[N]; //
      //cache objet jumeau
      cache_jumeau = new unsigned char[CACHE_JUMEAU_LEN]; //
  
      //reservation reference et "noir" camera pour une "tentative" de correction de la ref
      reference = new unsigned char[N]; //
      // 	unsigned char* noir_camera=new unsigned char[N];

      /// 4 tableaux pour transtypage des images 8 bits vers 64 bits (double).
      phased1 = new double[N];  //
      phased2 = new double[N];  //
      phased3 = new double[N];  //
      phased4 = new double[N];  //
  
      //reservation de 4 tableaux (image 2D)  pour le phase shifting
      phase1 = new unsigned char[N]; //
      phase2 = new unsigned char[N]; //
      phase3 = new unsigned char[N]; //
      phase4 = new unsigned char[N]; //
    }
  catch(bad_alloc ex)
    {cerr << endl << "memory allocation failure" << ex.what(); exit(EXIT_FAILURE); } 
  

  //g_entree_shift = new double[N];
#ifndef CUDA 
  fftw_complex* temp_fftw1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N); //
  fftw_complex* temp_fftw2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N); //
#endif  
  
  double *plan_reel, *plan_imag, *plan_reel_shift, *plan_imag_shift,	\
    *fft_reel_shift_norm, *fft_imag_shift_norm, *fft_reel_shift,	\
    *fft_imag_shift, *fft_module_shift, *fft_reel_tmp, *fft_imag_tmp, *fft_reel, *fft_imag;
  try{
    /** */
    plan_reel=new double[N]; //
    plan_imag=new double[N]; //

    //int N=g_window_dim_x*g_window_dim_y;
    plan_reel_shift=new double[N]; //
    plan_imag_shift=new double[N]; //

    fft_reel_shift_norm=new double[4*Nxmax*Nymax]; //
    fft_imag_shift_norm=new double[4*Nxmax*Nymax]; //

    /** */
    fft_reel_shift=new double[4*Nxmax*Nymax]; //
    fft_imag_shift=new double[4*Nxmax*Nymax]; // 
    fft_module_shift=new double[4*Nxmax*Nymax]; //

  
    //variables utilisées pour la TF2D

    //Avant Crop à NxMax;
    fft_reel_tmp = new double[N]; //
    fft_imag_tmp = new double[N]; //
    //Après Crop à NxMax;
    fft_reel = new double[4*Nxmax*Nymax]; //
    fft_imag = new double[4*Nxmax*Nymax]; //
  }
  catch(bad_alloc ex){
    cerr << endl << "memory allocation failure" << ex.what(); exit(EXIT_FAILURE); 
  }


  //double *qualite_phase_shift=new double[N];
  //espace3D reel, imaginaire et support de redondance
  
  //tester l'existence des fichiers au cas ou
  char *test_Copie;
  assert(test_Copie = (char *) calloc(strlen(Chemin) + GENERATED_FILENAME_LEN, sizeof(char)));
  //char test_FinNom[15];


  // ----------------------------------------------------------------------
  // initialisations
  // ----------------------------------------------------------------------

  int i_nb_proj = 0;
  int points_faux = 0;
  int jumeau_elimine = 0;  
  const int rayon_inf = xm0_limite * xm0_limite + ym0_limite * ym0_limite;
  int centres_exclus = 0;


  // taille totale du masque en pixel pour elimination objet jumeau;
  const int t_mask = 30;
  char* chaine = str_concat(g_MASKPATH, "/k_30x30.bmp");
  cache_jumeau = rempli_tableau( chaine, 0, 0, t_mask, t_mask);

  //str_concat(g_MASKPATH, "/k_30x30.bmp"), 


  switch(g_window_dim_x)
    {
    case 512:   
      masque = rempli_tableau( str_concat(g_MASKPATH, "/kmasq_512x512_30.bmp"),\
			       0, 0, g_window_dim_x, g_window_dim_y);
      printf("masque 512x512\n");
      break;
    case 256:
      masque = rempli_tableau( str_concat(g_MASKPATH, "/kmasq_256x256.bmp"), \
			       0, 0, g_window_dim_x, g_window_dim_y);
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

  CHRONO_SET(tf2d_time);
  CHRONO_SET(circ2d_time);
  CHRONO_SET(load2d_time);


  // pre-calcul du tableau des images à traiter
  int* ANGLES_TO_PROCESS;
  int ANGLES_size_alloc = 1 + (angle_count / angle_jump);
  int ANGLES_cardinal = 0;
  ARRAY_ALLOC(ANGLES_TO_PROCESS, ANGLES_size_alloc, int);
  for (int cpt_angle = angle_start; cpt_angle < angle_start + angle_count; cpt_angle = cpt_angle + angle_jump)
    {

      compute_radix_filename(cpt_angle, Chemin, test_Copie);
      if (ftest(test_Copie))
	{
	  ANGLES_TO_PROCESS[ANGLES_cardinal] = cpt_angle;
	  ANGLES_cardinal++;
	}
      
    }
  cout << endl << "Nombre de fichiers à traiter: " << ANGLES_cardinal;

  
  ///SAVE_INPUT_VOLUME
  FILE *output_vol;
  assert(output_vol = fopen(str_concat(g_OUTPUT_DIR, "/volume_realImg.raw"), "wb"));
  OUTPUT_TYPE* image_out;
  ARRAY_ALLOC(image_out, image_dim_x * image_dim_y, OUTPUT_TYPE);
  ///SAVE_INPUT_VOLUME
  

  // BOUCLE DE LECTURE
  for (int I = 0; I < ANGLES_cardinal; I++)	 
    {
      int cpt_angle = ANGLES_TO_PROCESS[I];

      if((cpt_angle - 100 * (cpt_angle / 100)) == 0)
	printf("cpt_angle=%i\n",cpt_angle);

      //----------------------------------------------------------------------
      //////////////////////////phase shifting////////////////////////////////
      //----------------------------------------------------------------------

      //**  On teste jusqu'à  angle_end fichiers
      //    On peut faire la reconstruction avec relativement peu de données,
      //    donc l'absence de certaines séries de coupes n'a aucune incidence
      //    grave et ne doit pas interrompre l'analyse
          

      compute_radix_filename(cpt_angle, Chemin, test_Copie);
      CHRONO_START(load2d_time);
      load_2D_image(phase1, 1, Chemin, cpt_angle, window_edge_x, window_edge_y, g_window_dim_x, g_window_dim_y);
      load_2D_image(phase2, 2, Chemin, cpt_angle, window_edge_x, window_edge_y, g_window_dim_x, g_window_dim_y);
      load_2D_image(phase3, 3, Chemin, cpt_angle, window_edge_x, window_edge_y, g_window_dim_x, g_window_dim_y);
      load_2D_image(phase4, 4, Chemin, cpt_angle, window_edge_x, window_edge_y, g_window_dim_x, g_window_dim_y);
      CHRONO_STOP(load2d_time);


      //** remise à zéro
      zerofill_double(plan_reel, N);
      zerofill_double(plan_imag, N);
      zerofill_double(plan_reel_shift, N);
      zerofill_double(plan_imag_shift, N);


      for(int pixel = 0; pixel < N; pixel++)
	{
	  plan_reel[pixel] = ((double)phase1[pixel] - (double)phase3[pixel]) * (masque[pixel]);
	  plan_imag[pixel] = ((double)phase4[pixel] - (double)phase2[pixel]) * (masque[pixel]);
	}

      if (SAVE_INPUT_VOLUME)
	{
	  convert_float(plan_reel, image_out, image_dim_x * image_dim_y);
	  fwrite(image_out, sizeof(OUTPUT_TYPE), image_dim_x * image_dim_y, output_vol);
	  convert_float(plan_imag, image_out, image_dim_x * image_dim_y);
	  fwrite(image_out, sizeof(OUTPUT_TYPE), image_dim_x * image_dim_y, output_vol);
	}

      
      //**  */
      CHRONO_START(circ2d_time);
      circshift_ccd2(plan_reel, plan_reel_shift, g_window_dim_x, g_window_dim_y, g_window_dim_x/2, g_window_dim_y/2);
      circshift_ccd2(plan_imag, plan_imag_shift, g_window_dim_x, g_window_dim_y, g_window_dim_x/2, g_window_dim_y/2);
      CHRONO_STOP(circ2d_time);


      // ----------------------------------------------------------------------
      // TF 2D du front d'onde: passage du plan image au plan réciproque///////
      // ----------------------------------------------------------------------

      //       TF2D(plan_reel_shift,plan_imag_shift,fft_reel_tmp,fft_imag_tmp,g_window_dim_x,g_window_dim_y);
      CHRONO_START(tf2d_time);

#ifdef CUDA
      TF2D_cuda(plan_reel_shift, plan_imag_shift, fft_reel_tmp, fft_imag_tmp, \
		g_window_dim_x, g_window_dim_y);
#else
      TF2D_eco(plan_reel_shift, plan_imag_shift, fft_reel_tmp, fft_imag_tmp, \
	       g_window_dim_x, g_window_dim_y, temp_fftw1, temp_fftw2);
#endif
      
      CHRONO_STOP(tf2d_time);

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


      int xi;
      int yi;
      // premier demi-espace /////////////////////////////////////////////
      for (xi=0;xi<Nxmax;xi++)
	{
	  for ( yi=0;yi<Nymax;yi++)
	    {
	      int cpt1=yi*g_window_dim_x+xi;
	      int cpt2=yi*(Nxmax)*2+xi;

	      fft_reel[cpt2]=fft_reel_tmp[cpt1];
	      fft_imag[cpt2]=fft_imag_tmp[cpt1];
	    }

	  for (int yi=g_window_dim_y-Nymax;yi<g_window_dim_y;yi++)
	    {
	      int cpt1=yi*g_window_dim_x+xi;

	      int cpt2 = xi+(yi-g_window_dim_y+2*Nymax)*2*(Nxmax);
	      fft_reel[cpt2]=fft_reel_tmp[cpt1];
	      fft_imag[cpt2]=fft_imag_tmp[cpt1];
	    }
	}

      // deuxieme demi-espace /////////////////////////////////////////////
      for(int xi=g_window_dim_x-Nxmax;xi<g_window_dim_x;xi++)
	{
	  for (int yi=0;yi<Nymax;yi++)
	    {
	      int cpt1=yi*g_window_dim_x+xi;
	      int cpt2=yi*(Nxmax)*2+(xi-g_window_dim_x+2*Nxmax);

	      fft_reel[cpt2]=fft_reel_tmp[cpt1];
	      fft_imag[cpt2]=fft_imag_tmp[cpt1];
	    }

	  for (int yi=g_window_dim_y-Nymax;yi<g_window_dim_y;yi++)
	    {
	      int cpt1=yi*g_window_dim_x+xi;

	      int cpt2 = (xi-g_window_dim_x+2*Nxmax)+(yi-g_window_dim_y+2*Nymax)*2*(Nxmax);
	      fft_reel[cpt2]=fft_reel_tmp[cpt1];
	      fft_imag[cpt2]=fft_imag_tmp[cpt1];
	    }
	}

      circshift_ccd2(fft_reel, fft_reel_shift, 2 * Nxmax, 2 * Nymax, Nxmax, Nymax);
      circshift_ccd2(fft_imag, fft_imag_shift, 2 * Nxmax, 2 * Nymax, Nxmax, Nymax);


      //----------------Recherche du maximum non centré------------------------


      //gloubiboulga function /!\ achtung, dont touch /!!\ Caramba, Né tochas
      
      int cpt_max=0;
      fft_module_shift[0]=pow(fft_reel_shift[0],2)+pow(fft_imag_shift[0],2);

      //Recherche du MAX dansle module
      for(int cpt=1;cpt<(4*Nxmax*Nymax);cpt++)
	{
	  fft_module_shift[cpt]=pow(fft_reel_shift[cpt],2)+pow(fft_imag_shift[cpt],2);
	  if(fft_module_shift[cpt]>fft_module_shift[cpt_max])
	    {
	      cpt_max=cpt;
	    }
	}

      /////////////////////////////
      //on met 1 dans l'image centre en cpt_max
      /////////////////////////////

      //centre[cpt_max]++;


      double max_part_reel = fft_reel_shift[cpt_max];
      double max_part_imag = fft_imag_shift[cpt_max];
      double max_module = fft_module_shift[cpt_max];

      // ----------------------------------------------------------------------
      ////////////////////////////normalisation par le pic central
      // ----------------------------------------------------------------------

      for(int cpt=0;cpt<(4*Nxmax*Nymax);cpt++)
	{
	  fft_reel_shift_norm[cpt]=(fft_reel_shift[cpt]*max_part_reel+fft_imag_shift[cpt]*max_part_imag)/max_module;
	  //fft_reel[cpt]=(fft_reel[cpt]*fft_reel[cpt_max]+fft_imag[cpt]*fft_imag[cpt_max])/fft_module[cpt_max];danger!!!!
	  fft_imag_shift_norm[cpt]=(fft_imag_shift[cpt]*max_part_reel-fft_reel_shift[cpt]*max_part_imag)/max_module;
	  //fft_imag[cpt]=(fft_imag[cpt]*fft_reel[cpt_max]-fft_reel[cpt]*fft_imag[cpt_max])/fft_module[cpt_max];danger!!!!
	}

      //--------------------------------------


      //Coordonnées dans l'espace 2D à partir de l'indice 1D: (xc,yc)=(0,0)=en haut à gauche de l'image
      int xmi = cpt_max % (2 * Nxmax);
      int ymi = cpt_max / (2 * Nymax);

      int xc = xmi - Nxmax;
      int yc = ymi - Nymax;
      //coordonnée objet jumeaux
      int xmij = 2 * Nxmax-xmi;
      int ymij = 2 * Nymax-ymi;
      ////////////////////////////////////Virer l'objet jumeaux
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


      //----------------------------------------------------------------------
      ////////////////////////////Mapping 3D
      //----------------------------------------------------------------------


      //coordonnée dans l'image2D centrée (xm0,ym0)=(0,0)=au centre de l'image
      int xm0=(xmi-Nxmax);//
      //printf("xm0 %i \n", xm0);
      int ym0=(ymi-Nymax);

      if(xm0==0 && ym0==0)
	printf("(xm0,ym0)=(0,0) dans le plan %i\n", cpt_angle);

      if((xm0*xm0+ym0*ym0) > rayon_inf)
	{
	  //printf("xm0 : %i, ym0:%i\n",xm0,ym0);
	  //printf("\nFichier %i sans signal\n", cpt_angle);
	  centres_exclus++;
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
	      i_nb_proj++;

	     
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
			  sup_redon[k] += 1;//pour calculer le support

			  reel_arc[k] += fft_reel_shift_norm[cpt];//pour calculer l'image
			  imag_arc[k] += fft_imag_shift_norm[cpt];//pour calculer l'image

			}
		      else
			points_faux++;
		    } //fin for y
		}

	    }//fin if zm0>-1
	}//fin else xm0_limite

      temps_final = clock ();
      temps_cpu = (temps_final - temps_initial) / CLOCKS_PER_SEC;
      temps_total = temps_total + temps_cpu;
      temps_initial = clock();
      
    }


  // fin boucle principale de lecture sur tous les angles
  // END OF LOOP
  
  CHRONO_END(tf2d_time, "Temps passé à invoquer la TF2D");
  CHRONO_END(circ2d_time, "Temps passé à circshifter en 2D");
  CHRONO_END(load2d_time, "Temps passé à lire les images 2D");

  if (SAVE_INPUT_VOLUME)
    fclose(output_vol);

  // ----------------------------------------------------------------------
  // reporting, mesure du temps de calcul
  // ----------------------------------------------------------------------

  
  printf("jumeaux re-elimines : %i\n",jumeau_elimine);
  printf("i_nb_proj: %i\n", i_nb_proj);
  printf("Nombre de centre exclus : %i\n", centres_exclus);  
  printf("points_faux : %i,points_faux\n",points_faux);
  printf("temps_total proj : %f \n",temps_total);
  printf("coucou \n");

  
  temps_final = clock ();
  temps_cpu = (temps_final - temps_initial) / CLOCKS_PER_SEC / i_nb_proj;
  printf("temps moyen pour 1 angle: %f\n", temps_cpu);
  temps_cpu = (temps_final - temps_initial) / CLOCKS_PER_SEC;
  printf("temps total pour %i angle(s): %f\n", i_nb_proj, temps_cpu);


  
  // ----------------------------------------------------------------------
  // desallocations
  // ----------------------------------------------------------------------


  delete masque; delete cache_jumeau; //

  delete phased1; delete phased2; delete phased3; delete phased4; //

  delete phase1; delete phase2; delete phase3; delete phase4; //
  delete[] reference; //
 
  delete fft_reel; delete fft_imag; //
  delete fft_reel_tmp; delete fft_imag_tmp; //

  delete fft_reel_shift; delete fft_imag_shift; delete fft_module_shift; //
  
  delete fft_reel_shift_norm; delete fft_imag_shift_norm; //

  delete plan_reel; delete plan_imag; //

  delete plan_reel_shift; delete plan_imag_shift; //


#ifndef CUDA
  fftw_free( temp_fftw1 );//
  fftw_free( temp_fftw2 ); //
#endif

}
