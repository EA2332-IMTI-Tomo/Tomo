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
#include "util_Magick.h"
#include "util.h"
#include "readloop.h"
#include "memory.h"
#include "fourier.h"
#ifdef CUDA 
#include "cuFourier.h"
#endif



extern bool read_saved_volume;

void 
readLoop_binaryVolume(RECON_TYPE* reel_arc, RECON_TYPE* imag_arc, unsigned short int* sup_redon, int* centre, \
		      const int Nxmax, const int Nymax, int Nxmax_Rf, \
		      const int coin_x, const int coin_y,				\
		      const int xm0_limite, const int ym0_limite,			\
		      float rayon, float delta_zmax, int NbAngle,                    \
		      int image_dimx, int image_dimy,                                \
		      char* Chemin)
{    

  //nombre de pixel (pour le tableau)
  const size_t s_IMG_CCD = g_dimx_ccd * g_dimy_ccd;
  const size_t s_IMG_INPUT = image_dimx * image_dimy;
  const size_t s_FOURIER_SUPPORT = 4 * Nxmax * Nymax; //Nxmax = 1/2 ouverture fourier

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
  RECON_TYPE *phased1, *phased2, *phased3, *phased4;
  
  //cache objet jumeau
  ARRAY_ALLOC(cache_jumeau, CACHE_JUMEAU_LEN, unsigned char);
  //cache_jumeau = new unsigned char[CACHE_JUMEAU_LEN]; //
  
  //reservation reference et "noir" camera pour une "tentative" de correction de la ref
  ARRAY_ALLOC(reference, s_IMG_CCD, unsigned char);
  //reference = new unsigned char[s_IMG_CCD]; //
  // 	unsigned char* noir_camera=new unsigned char[N];

  //g_entree_shift = new RECON_TYPE[s_IMG_CCD];
#ifndef CUDA 
  fftw_complex* temp_fftw1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * s_IMG_CCD); //
  fftw_complex* temp_fftw2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * s_IMG_CCD); //
#endif  
  
  RECON_TYPE *plan_reel, *plan_imag, *plan_reel_shift, *plan_imag_shift,	\
    *fft_reel_shift_norm, *fft_imag_shift_norm, *fft_reel_shift,	\
    *fft_imag_shift, *fft_module_shift, *fft_reel_tmp, *fft_imag_tmp, *fft_reel, *fft_imag;

  ARRAY_ALLOC(plan_reel, s_IMG_CCD, RECON_TYPE);
  ARRAY_ALLOC(plan_imag, s_IMG_CCD, RECON_TYPE);
  ARRAY_ALLOC(plan_reel_shift, s_IMG_CCD, RECON_TYPE);
  ARRAY_ALLOC(plan_imag_shift, s_IMG_CCD, RECON_TYPE);

  ARRAY_ALLOC(fft_reel_shift_norm, s_FOURIER_SUPPORT, RECON_TYPE);
  ARRAY_ALLOC(fft_imag_shift_norm, s_FOURIER_SUPPORT, RECON_TYPE);
  ARRAY_ALLOC(fft_reel_shift, s_FOURIER_SUPPORT, RECON_TYPE);
  ARRAY_ALLOC(fft_imag_shift, s_FOURIER_SUPPORT, RECON_TYPE);
  ARRAY_ALLOC(fft_module_shift, s_FOURIER_SUPPORT, RECON_TYPE);

  
  //variables utilisées pour la TF2D
  
  //Avant Crop à NxMax;
  ARRAY_ALLOC(fft_reel_tmp, s_IMG_CCD, RECON_TYPE);
  ARRAY_ALLOC(fft_imag_tmp, s_IMG_CCD, RECON_TYPE);

  //Après Crop à NxMax;
  ARRAY_ALLOC(fft_reel, s_FOURIER_SUPPORT, RECON_TYPE);
  ARRAY_ALLOC(fft_imag, s_FOURIER_SUPPORT, RECON_TYPE);
  //fft_reel = new RECON_TYPE[s_FOURIER_SUPPORT]; //
  //fft_imag = new RECON_TYPE[s_FOURIER_SUPPORT]; //
  

  //RECON_TYPE *qualite_phase_shift=new RECON_TYPE[N];
  //espace3D reel, imaginaire et support de redondance
  
  //tester l'existence des fichiers au cas ou
  char *test_Copie;
  assert(test_Copie = (char *) calloc(strlen(Chemin) + GENERATED_FILENAME_LEN, sizeof(char)));
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
  char* chaine = str_concat(g_MASKPATH, "/k_30x30.bmp");
  cache_jumeau = rempli_tableau( chaine, 0, 0, t_mask, t_mask);

  //str_concat(g_MASKPATH, "/k_30x30.bmp"), 


  switch(g_dimx_ccd)
    {
    case 512:   
      masque = rempli_tableau( str_concat(g_MASKPATH, "/kmasq_512x512_30.bmp"),\
			       0, 0, g_dimx_ccd, g_dimy_ccd);
      printf("masque 512x512\n");
      break;
    case 256:
      masque = rempli_tableau( str_concat(g_MASKPATH, "/kmasq_256x256.bmp"), \
			       0, 0, g_dimx_ccd, g_dimy_ccd);
      printf("masque 256x256\n");
      break;
    default:
      fprintf(stderr, "unsupported value for variable: g_dimx_ccd");
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

  /**

   // pre-calcul du tableau des images à traiter
   int* ANGLES_TO_PROCESS;
   int ANGLES_size_alloc = 1 + (NbAngle / SautAngle);
   int ANGLES_cardinal = 0;
   ARRAY_ALLOC(ANGLES_TO_PROCESS, ANGLES_size_alloc, int);
   for (int cpt_angle = premier_plan; cpt_angle < premier_plan + NbAngle; cpt_angle = cpt_angle + SautAngle)
   {

   compute_radix_filename(cpt_angle, Chemin, test_Copie);
   if (ftest(test_Copie))
   {
   ANGLES_TO_PROCESS[ANGLES_cardinal] = cpt_angle;
   ANGLES_cardinal++;
   }
      
   }
   cout << endl << "Nombre de fichiers à traiter: " << ANGLES_cardinal;
  
  */


  /**

   // BOUCLE DE LECTURE
   for (int I = 0; I < ANGLES_cardinal; I++)	 
   {
   int cpt_angle = ANGLES_TO_PROCESS[I];

   if((cpt_angle - 100 * (cpt_angle / 100)) == 0)
   printf("cpt_angle=%i\n",cpt_angle);


   //----------------------------------------------------------------------
   //////////////////////////phase shifting////////////////////////////////
   //----------------------------------------------------------------------

   //**  On teste jusqu'à  Num_Angle_final fichiers
   //    On peut faire la reconstruction avec relativement peu de données,
   //    donc l'absence de certaines séries de coupes n'a aucune incidence
   //    grave et ne doit pas interrompre l'analyse
          

   compute_radix_filename(cpt_angle, Chemin, test_Copie);
   CHRONO_START(load2d_time);
   load_2D_image(phase1, 1, Chemin, cpt_angle, coin_x, coin_y, g_dimx_ccd, g_dimy_ccd);
   load_2D_image(phase2, 2, Chemin, cpt_angle, coin_x, coin_y, g_dimx_ccd, g_dimy_ccd);
   load_2D_image(phase3, 3, Chemin, cpt_angle, coin_x, coin_y, g_dimx_ccd, g_dimy_ccd);
   load_2D_image(phase4, 4, Chemin, cpt_angle, coin_x, coin_y, g_dimx_ccd, g_dimy_ccd);
   CHRONO_STOP(load2d_time);


   //** remise à zéro
   zerofill_RECON_TYPE(plan_reel, N);
   zerofill_RECON_TYPE(plan_imag, N);

   for(int pixel = 0; pixel < N; pixel++)
   {
   plan_reel[pixel] = ((RECON_TYPE)phase1[pixel] - (RECON_TYPE)phase3[pixel]) * (masque[pixel]);
   plan_imag[pixel] = ((RECON_TYPE)phase4[pixel] - (RECON_TYPE)phase2[pixel]) * (masque[pixel]);
   }

      
   zerofill_RECON_TYPE(plan_reel_shift, N);
   zerofill_RECON_TYPE(plan_imag_shift, N);

  **/





  //************************************************************************************************************************

  FILE* fp;
  assert(fp = fopen(Chemin, "rb"));
  INPUT_TYPE* readImage;
  RECON_TYPE* readImageD; 
  ARRAY_ALLOC(readImage, image_dimx * image_dimy, INPUT_TYPE); 
  ARRAY_ALLOC(readImageD, image_dimx * image_dimy, RECON_TYPE); 
  int cpt_angle;
  

  // BOUCLE DE LECTURE
  for (int I = 0; I < NbAngle; I++)
    {
      cpt_angle = I;

      zerofill_double(plan_reel, s_IMG_CCD);
      zerofill_double(plan_imag, s_IMG_CCD);
      zerofill_double(plan_reel_shift, s_IMG_CCD);
      zerofill_double(plan_imag_shift, s_IMG_CCD);

      if (read_saved_volume)
	assert(s_IMG_INPUT == s_IMG_CCD);

      CHRONO_START(load2d_time);
      // lecture plan reel
      assert(fread(readImage, sizeof(INPUT_TYPE), s_IMG_INPUT, fp) != EOF);
      convert_double(readImage, readImageD, s_IMG_INPUT);
      extract_subImage(readImageD, plan_reel, image_dimx, image_dimy, coin_x, coin_y, g_dimx_ccd, g_dimy_ccd);
      // lecture plan image
      assert(fread(readImage, sizeof(INPUT_TYPE), s_IMG_INPUT, fp) != EOF);
      convert_double(readImage, readImageD, s_IMG_INPUT);
      extract_subImage(readImageD, plan_imag, image_dimx, image_dimy, coin_x, coin_y, g_dimx_ccd, g_dimy_ccd);
      CHRONO_STOP(load2d_time);
	    

      // à sauter si on recharge un volume déjà sauvegardé: la sauvegarde 
      // s'est déjà déroulée avec application du masque
      if (! read_saved_volume) {
	//multiplication par masqueuniquement
	for (int pixel = 0; pixel < s_IMG_CCD; pixel++)
	  {
	    plan_reel[pixel] *= masque[pixel];
	    plan_imag[pixel] *= masque[pixel];
	  }
      }

      //************************************************************************************************************************

      //**  */
      CHRONO_START(circ2d_time);
      circshift2D_memcpy(plan_reel, plan_reel_shift, g_dimx_ccd, g_dimy_ccd);//, g_dimx_ccd/2, g_dimy_ccd/2);
      circshift2D_memcpy(plan_imag, plan_imag_shift, g_dimx_ccd, g_dimy_ccd);//, g_dimx_ccd/2, g_dimy_ccd/2);
      CHRONO_STOP(circ2d_time);



      // ----------------------------------------------------------------------
      // TF 2D du front d'onde: passage du plan image au plan réciproque///////
      // ----------------------------------------------------------------------

      //       TF2D(plan_reel_shift,plan_imag_shift,fft_reel_tmp,fft_imag_tmp,g_dimx_ccd,g_dimy_ccd);
      CHRONO_START(tf2d_time);      
#ifdef CUDA
      TF2D_cuda(plan_reel_shift, plan_imag_shift, fft_reel_tmp, fft_imag_tmp, \
		g_dimx_ccd, g_dimy_ccd);
#else
      
      TF2D_eco(plan_reel_shift, plan_imag_shift, fft_reel_tmp, fft_imag_tmp, \
	       g_dimx_ccd, g_dimy_ccd, temp_fftw1, temp_fftw2);
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
	      int cpt1=yi*g_dimx_ccd+xi;
	      int cpt2=yi*(Nxmax)*2+xi;

	      fft_reel[cpt2]=fft_reel_tmp[cpt1];
	      fft_imag[cpt2]=fft_imag_tmp[cpt1];
	    }

	  for (int yi=g_dimy_ccd-Nymax;yi<g_dimy_ccd;yi++)
	    {
	      int cpt1=yi*g_dimx_ccd+xi;

	      int cpt2 = xi+(yi-g_dimy_ccd+2*Nymax)*2*(Nxmax);
	      fft_reel[cpt2]=fft_reel_tmp[cpt1];
	      fft_imag[cpt2]=fft_imag_tmp[cpt1];
	    }
	}

      // deuxieme demi-espace /////////////////////////////////////////////
      for(int xi=g_dimx_ccd-Nxmax;xi<g_dimx_ccd;xi++)
	{
	  for (int yi=0;yi<Nymax;yi++)
	    {
	      int cpt1=yi*g_dimx_ccd+xi;
	      int cpt2=yi*(Nxmax)*2+(xi-g_dimx_ccd+2*Nxmax);

	      fft_reel[cpt2]=fft_reel_tmp[cpt1];
	      fft_imag[cpt2]=fft_imag_tmp[cpt1];
	    }

	  for (int yi=g_dimy_ccd-Nymax;yi<g_dimy_ccd;yi++)
	    {
	      int cpt1=yi*g_dimx_ccd+xi;

	      int cpt2 = (xi-g_dimx_ccd+2*Nxmax)+(yi-g_dimy_ccd+2*Nymax)*2*(Nxmax);
	      fft_reel[cpt2]=fft_reel_tmp[cpt1];
	      fft_imag[cpt2]=fft_imag_tmp[cpt1];
	    }
	}

      circshift2D_memcpy(fft_reel, fft_reel_shift, 2 * Nxmax, 2 * Nymax);//, Nxmax, Nymax);
      circshift2D_memcpy(fft_imag, fft_imag_shift, 2 * Nxmax, 2 * Nymax);//, Nxmax, Nymax);


      //----------------Recherche du maximum non centré------------------------


      //gloubiboulga function /!\ achtung, dont touch /!!\ Caramba, Né tochas
      
      int position_max_module=0;
      fft_module_shift[0]=pow(fft_reel_shift[0],2)+pow(fft_imag_shift[0],2);

      //Recherche du MAX dansle module
      for(int cpt = 1; cpt < s_FOURIER_SUPPORT; cpt++)
	{
	  fft_module_shift[cpt]=pow(fft_reel_shift[cpt],2)+pow(fft_imag_shift[cpt],2);
	  if(fft_module_shift[cpt]>fft_module_shift[position_max_module])
	    position_max_module=cpt;
	  
	}

      /////////////////////////////
      //on met 1 dans l'image centre en position_max_module
      /////////////////////////////

      //centre[position_max_module]++;


      RECON_TYPE max_part_reel = fft_reel_shift[position_max_module];
      RECON_TYPE max_part_imag = fft_imag_shift[position_max_module];
      RECON_TYPE max_module = fft_module_shift[position_max_module];


      // ----------------------------------------------------------------------
      ////////////////////////////normalisation par le pic central
      // ----------------------------------------------------------------------
      
      assert(max_module);
      for(int cpt = 0; cpt < s_FOURIER_SUPPORT; cpt++)
	{
	  fft_reel_shift_norm[cpt]=(fft_reel_shift[cpt]*max_part_reel+fft_imag_shift[cpt]*max_part_imag) / max_module;
	  fft_imag_shift_norm[cpt]=(fft_imag_shift[cpt]*max_part_reel-fft_reel_shift[cpt]*max_part_imag) / max_module;
	}

      //fft_reel[cpt]=(fft_reel[cpt]*fft_reel[position_max_module]+fft_imag[cpt]*fft_imag[position_max_module])/fft_module[position_max_module];danger!!!!
      //fft_imag[cpt]=(fft_imag[cpt]*fft_reel[position_max_module]-fft_reel[cpt]*fft_imag[position_max_module])/fft_module[position_max_module];danger!!!!

      //--------------------------------------


      //Coordonnées dans l'espace 2D à partir de l'indice 1D: (xc,yc)=(0,0)=en haut à gauche de l'image
      assert(Nxmax * Nymax);
      int xmi = position_max_module % (2 * Nxmax);
      int ymi = position_max_module / (2 * Nymax);

      int xc = xmi - Nxmax;
      int yc = ymi - Nymax;
      //coordonnée objet jumeaux
      int xmij = 2 * Nxmax-xmi;
      int ymij = 2 * Nymax-ymi;
      ////////////////////////////////////Virer l'objet jumeaux
      if((xc*xc+yc*yc) > 900)//35*35=1225;objet jumeau pas trop près de l'objet
	{ 
	  jumeau_elimine++;

	  for(int x = xmij - t_mask/2; x < xmij + t_mask/2; x++)
	    {
	      for(int y = ymij - t_mask/2; y < ymij + t_mask/2; y++)
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
  
  //*************************************************************************************************
  fclose(fp);
  free(readImage);


  CHRONO_END(tf2d_time, "Temps passé à invoquer la TF2D");
  CHRONO_END(circ2d_time, "Temps passé à circshifter en 2D");
  CHRONO_END(load2d_time, "Temps passé à lire les images 2D");

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
  /** */ 
  if ( CLOCKS_PER_SEC * nb_proj)
    {
      temps_cpu = (temps_final - temps_initial) / CLOCKS_PER_SEC / nb_proj;
      printf("temps moyen pour 1 angle: %f\n", temps_cpu);
      temps_cpu = (temps_final - temps_initial) / CLOCKS_PER_SEC;
      printf("temps total pour %i angle(s): %f\n", nb_proj, temps_cpu);
    }

  
  // ----------------------------------------------------------------------
  // desallocations
  // ----------------------------------------------------------------------

  
#ifndef CUDA
  fftw_free( temp_fftw1 );
  fftw_free( temp_fftw2 );
#endif


  free(cache_jumeau); //
  free(reference); //
  free(plan_reel); free(plan_imag); //
  free(plan_reel_shift); free(plan_imag_shift); //

  free(fft_reel_shift); free(fft_imag_shift); free(fft_module_shift); //
  free(fft_reel_shift_norm); 
  free(fft_imag_shift_norm); //


 
  free(fft_reel); free(fft_imag); //
  free(fft_reel_tmp); free(fft_imag_tmp); //


  cerr << "desallocations";
    
}
