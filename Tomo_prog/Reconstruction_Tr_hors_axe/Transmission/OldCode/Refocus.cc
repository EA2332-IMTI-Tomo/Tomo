#include "vChrono.h" // requiert boost
#include "vChronos.h" 

#include <iostream>
#include <boost/thread.hpp>   
#include <boost/date_time.hpp>

#include "Refocus.h"


#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>
#include <exception>


#include "main.h"
#include "util.h"
//#include "readloop.h"
#include "IO.h"
#include "Holo_Process.h"
#include "Compute.h"
#include "cvDisplayVolume.h"


// fait tout déconner à l'édition de lien
/*
  #include <AIR_Volume.h>
  #include <COMPLEX_Volume.h>
  #include <FFTW_Volume.h>
  #include <FFTW_Image.h>
*/


using namespace std;

#include "cv.h"
#include "highgui.h"

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
using namespace cv;




//******************************************************************************
// fonctions locales à implémenter 
//******************************************************************************


// circshift 2D non centré (permutation des rectangles au lieu de carrés)
void circshift2(double* entree, double* result, Var2D dim,Var2D decal);

// possible de passer par FFTW_Image et d'importer-exporter sans shift
void TF2D_INV(double entree_reelle[],double entree_imag[], double sortie_reelle[],double sortie_imag[],int taille_x,int taille_y);
void TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y);


// vol: la stack des plans pour diverses valeurs de retropropag(delta, spectre)
void Plan_ds_Vol(double *Vol3DRe, double *Vol3DIm, double *plan2DRe, double *plan2DIm, Var3D dimVol, Var2D dimPlan, int z3Di);

void read_dataM(char *inputfile,Mat* Data, int image_size, int IMAGE_NUM);



//******************************************************************************
// 
//******************************************************************************


void AutofocusGrad(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe, double *sortieIm, size_t dimX, size_t dimY, double rayon, size_t SpecPosX, size_t SpecPosY, int cpt_angle)
{
  Var2D size = {(int)dimX, (int)dimY};
  Var2D specPos = {(int)SpecPosX, (int)SpecPosY};
  AutofocusGrad(fft_shift_normRe, fft_shift_normIm, sortieRe, sortieIm, size, rayon, specPos, cpt_angle);
}


// **************************************************


//  Var2D dimInit={2*NXMAX,2*NYMAX},dimFin={4*NXMAX,4*NYMAX};
void AutofocusGrad(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe, double *sortieIm, Var2D dimInit, double rayon, Var2D CoordSpecI, int cpt_angle)
{

  int zInf=-5, zSup=5, NbPlanZ=zSup-zInf;
  
  int xm0=(CoordSpecI.x-dimInit.x/2), ym0=(CoordSpecI.y-dimInit.x/2);
  Var2D decal={-xm0+dimInit.x/2,-ym0+dimInit.x/2};  ///Eliminer les franges d'inclinaision
  Var2D plandecal={dimInit.x/2, dimInit.x/2};

  Var3D dimVol={dimInit.x,dimInit.y,NbPlanZ};


  double *spectrePropag_normRe=new double[dimInit.x*dimInit.x];
  double *spectrePropag_normIm=new double[dimInit.x*dimInit.x];
  double *spectrePropag_normRe_shift=new double[dimInit.x*dimInit.x];
  double *spectrePropag_normIm_shift=new double[dimInit.x*dimInit.x];

  double *planObjetRe = new double[dimInit.x*dimInit.x];
  double *planObjetIm = new double[dimInit.x*dimInit.x];
  double *planObjetRe_shift = new double[dimInit.x*dimInit.x];
  double *planObjetIm_shift = new double[dimInit.x*dimInit.x];

  // stack d'images pour toutes les valeurs de focus
  double *stack_planObjets_Re = new double[dimInit.x*dimInit.x*NbPlanZ];
  double *stack_planObjets_Im = new double[dimInit.x*dimInit.x*NbPlanZ];
  //
  double *stack_planObjets_refocal_Re = new double[dimInit.x*dimInit.x];
  double *stack_planObjets_refocal_Im = new double[dimInit.x*dimInit.x];

  double *spectreRe_refocal=new double[dimInit.x*dimInit.x];
  double *spectreIm_refocal=new double[dimInit.x*dimInit.x];
  double *Gradient= new double[NbPlanZ];

  

  for(int deltaZ = zInf; deltaZ < zSup ; deltaZ++)
    {
      // --------------------------------------------------
      // calcul du spectre rétropropagé courant
      // paramètre: deltaZ
      retroPropagSA(deltaZ, fft_shift_normRe, fft_shift_normIm, spectrePropag_normRe, spectrePropag_normIm, \
		    dimInit, rayon, CoordSpecI );
      
      ///Eliminer les franges d'inclinaision
      
      //  --------------------------------------------------
      // on repasse en espace image: tf2d inv, précédée de circshift
      circshift2(spectrePropag_normRe, spectrePropag_normRe_shift, dimInit, decal);
      circshift2(spectrePropag_normIm, spectrePropag_normIm_shift, dimInit, decal);

      TF2D_INV(spectrePropag_normRe_shift, spectrePropag_normIm_shift, planObjetRe, planObjetIm, dimInit.x, dimInit.y);

      // circshift final différent 
      circshift2(planObjetRe, planObjetRe_shift ,dimInit,plandecal);
      circshift2(planObjetIm, planObjetIm_shift ,dimInit,plandecal);
      // on est alors en espace image
      //  --------------------------------------------------

      // on ajoute le plan courant (f. de deltaZ) dans la stack des plans 
      Plan_ds_Vol(stack_planObjets_Re, stack_planObjets_Im,  \
		  planObjetRe_shift, planObjetIm_shift,	     \
		  dimVol, dimInit, deltaZ);

      // ecrit dans Gradient[] la valeur du plan courant
      Grad(planObjetRe_shift, planObjetIm_shift, Gradient, dimInit.x, NbPlanZ, deltaZ+zSup);
    }
  delete[] planObjetRe, planObjetIm, planObjetRe_shift, planObjetIm_shift;


  ///############################################
  // on cherche la position du minimum de Gradient[]
  int cpt_min=0;
  for(int cpt=0;cpt<NbPlanZ;cpt++)
    {
      if(Gradient[cpt] < Gradient[cpt_min])
	{
	  Gradient[cpt_min]=Gradient[cpt];
	  cpt_min=cpt;
	}
    }
      

  cout<<"cpt_min="<<cpt_min+1<<","<<"Gradient="<<Gradient[cpt_min]<<endl;
  //cout<<"cpt_max="<<cpt_max+1<<","<<"Gradient="<<Gradient[cpt_max]<<endl;
  delete[] Gradient;


  double *objet3DReTmp=new double[dimInit.x*dimInit.x*NbPlanZ];
  double *objet3DImTmp=new double[dimInit.x*dimInit.x*NbPlanZ];
  //memset(objet3DReTmp, 0, dimInit.x*dimInit.x*8);
  //memset(objet3DImTmp, 0, dimInit.x*dimInit.x*8);
  for(int cpt=0; cpt<dimInit.x*dimInit.x*NbPlanZ; cpt++)
    {
      objet3DReTmp[cpt]=stack_planObjets_Re[cpt];
      objet3DImTmp[cpt]=stack_planObjets_Im[cpt];
    }
  //SAV(objet3DReTmp, dimInit.x*dimInit.x*NbPlanZ, "/home/hui/maniptomo/IDP/champ_2d/Objet3DPlanRepro.bin", FLOAT,"a+b");


  FILE* fichier_objet3DRe = NULL;
  FILE* fichier_objet3DIm = NULL;
  fichier_objet3DRe = fopen("/ramdisk/objet3DRe_refocal.bin", "wb");
  fichier_objet3DIm = fopen("/ramdisk/objet3DIm_refocal.bin", "wb");

  ASSERT( fichier_objet3DRe );
  ASSERT( fichier_objet3DIm );

  int cpt=0;
  float precision;
  for(int z=cpt_min; z < cpt_min+1; z++)
    {
      for(int y=0; y < dimInit.x; y++)
	{
	  for(int x=0; x < dimInit.x; x++)
	    {
	      cpt= x + y * dimInit.x + z * dimInit.x * dimInit.x;
	      precision=objet3DReTmp[cpt];
	      fwrite(&precision,sizeof(precision),1,fichier_objet3DRe);
	      precision=objet3DImTmp[cpt];
	      fwrite(&precision,sizeof(precision),1,fichier_objet3DIm);
	    }
	}
    }
  fclose(fichier_objet3DRe);
  fclose(fichier_objet3DIm);
  //printf("ecriture de objet3DRe_refocal.bin et objet3DIm_refocal.bin : OK \n");
  delete[] objet3DReTmp, objet3DImTmp, stack_planObjets_Re, stack_planObjets_Im;
  delete[] spectrePropag_normRe, spectrePropag_normIm;


  Mat *Re= new Mat[NbPlanZ];
  Mat *Im= new Mat[NbPlanZ];
  read_dataM("/ramdisk/objet3DRe_refocal.bin", Re, dimInit.x, NbPlanZ);
  read_dataM("/ramdisk/objet3DIm_refocal.bin", Im, dimInit.x, NbPlanZ);


  ///Copier matrice de opencv vers double///
  for(int y=0; y<dimInit.x; y++)
    {
      for(int x=0; x<dimInit.x; x++)
	{
	  int cpt=y*dimInit.x+x;
	  stack_planObjets_refocal_Re[cpt]=Re[0].at<double>(y,x);
	  stack_planObjets_refocal_Im[cpt]=Im[0].at<double>(y,x);
	}
    }
  delete[] Re, Im;


  //SAV(objet3DRe_refocal, dimInit.x*dimInit.y, "/home/hui/maniptomo/IDP/champ_2d/objet3DPlanRecal.bin", FLOAT,"a+b");
  TF2D(stack_planObjets_refocal_Re, stack_planObjets_refocal_Im, spectreRe_refocal, spectreIm_refocal, dimInit.x, dimInit.y);
  delete[] stack_planObjets_refocal_Re, stack_planObjets_refocal_Im;


  Var2D decalI={xm0-dimInit.x/2,ym0-dimInit.x/2};
  circshift2(spectreRe_refocal, sortieRe, dimInit, decalI); ///Ajouter les franges d'inclinaision
  circshift2(spectreIm_refocal, sortieIm, dimInit, decalI);
  delete[] spectreRe_refocal, spectreIm_refocal;
}



// *****************************************************************************
//  pas utilisé!
// *****************************************************************************


int 
retroPropag(double* spectre3D_Re, double* spectre3D_Im, double *sup_redon, \
	    int dim_final, double* Spectre2D_Re, double* Spectre2D_Im, \
	    Var2D posSpec, Var3D decal,Var2D NMAX, double rayon)
{

  int k_T=0, k_R=0;
  int ptsExcluHolo=0;
  int dimPlanFinal = round(dim_final*dim_final),
    xm0 = posSpec.x, 
    ym0 = posSpec.y,
    dimVolX = round(dim_final);

  ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
  double r2=rayon*rayon;
  double arg_z_arc=0,z_arc=0,zm0;
  
  //printf("round(rayon*rayon-(xm0)^2-(ym0)^2: %i\n",round(rayon*rayon-(xm0)^2-(ym0)^2));
  double zm0_carre = rayon*rayon-xm0*xm0-ym0*ym0;
  
  if( round(zm0_carre)> -1)
    {
      zm0=sqrt(zm0_carre);
      int NXMAX_CARRE=NMAX.x*NMAX.x;
      for (int y = -NMAX.y; y < NMAX.y; y++)
	{ //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
	  int y_carre=y*y;
	  for (int x = -NMAX.x; x < NMAX.x; x++)
	    { //on balaye l'image 2D en y, centre au milieu
	      int cpt=(y+NMAX.y)*2*NMAX.x+x+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D
	      if(x*x+y_carre<NXMAX_CARRE)
		{ //ne pas depasser l'ouverture numérique pour 1 hologramme
		  double z_carre=r2-x*x-y_carre; //altitude au carré des données
		  double z_T=round(sqrt(z_carre)-zm0);
		  double altitude_T=(z_T+decal.z)*dimPlanFinal; //donne n'importequoi sans l'arrondi sur z!!
		  k_T=(-xm0+x+decal.x)+(-ym0+y+decal.y)*dimVolX+round(altitude_T);//indice du tableau 1D du volume 3D
		  spectre3D_Re[k_T]+=Spectre2D_Re[cpt];//pour calculer l'image;//
		  spectre3D_Im[k_T]+=Spectre2D_Im[cpt];//pour calculer l'image
		  sup_redon[k_T]+=1;//pour calculer le support

		}
	      else
		ptsExcluHolo++;
	    }
	}

    }
  return ptsExcluHolo;
}



// *****************************************************************************
//
// *****************************************************************************

//  Var2D dimInit={2*NXMAX,2*NYMAX},dimFin={4*NXMAX,4*NYMAX};
// NxMax = NyMax = 121 (242)
//retroPropag(deltaZ, <images>, dimInit, rayon, CoordSpecI );

void  
retroPropagSA(float deltaZ, double *fft_shift_normRe, double *fft_shift_normIm, \
	      double* spectrePropag_normRe, double* spectrePropag_normIm, \
	      Var2D dimSpctHolo, double rayon, Var2D CoordSpecI)
{
  deltaZ=deltaZ/20;
  Var2D KMAX={round(dimSpctHolo.x/2),round(dimSpctHolo.y/2)};
  const int NbPixHolo=dimSpctHolo.x*dimSpctHolo.y;
  double kz[dimSpctHolo.x*dimSpctHolo.y];
  nbCplx rephase[dimSpctHolo.x*dimSpctHolo.y];
  nbCplx spectrePropag[dimSpctHolo.x*dimSpctHolo.y];
  Var2D CoordSpecC={CoordSpecI.x-KMAX.x,CoordSpecI.y-KMAX.y};

  // nbCplx planObjet_I[dimSpctHolo.x*dimSpctHolo.y];


  ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
  double r2=rayon*rayon;
  double arg_kz_arc=0,kz_arc=0;
  //printf("round(rayon*rayon-(xm0)^2-(ym0)^2: %i\n",round(rayon*rayon-(xm0)^2-(ym0)^2));

  double kz0 = sqrt(rayon*rayon-CoordSpecC.x*CoordSpecC.x-CoordSpecC.y*CoordSpecC.y);
  //cout<<"kz0="<<kz0<<endl;
  int NXMAX_CARRE=KMAX.x*KMAX.x;

  ///CALCULER LES KZ-----------------------------------
  //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
  for (int ky = -KMAX.y; ky < KMAX.y; ky++) { 
    int ky_carre=ky*ky;
    //on balaye l'image 2D en y, centre au milieu
    for (int kx = -KMAX.x; kx < KMAX.x; kx++) { 
      //kyi*dim+kxi ...calcul du cpt du tableau 1D de l'image 2D
      int cpt=(ky+KMAX.y)*dimSpctHolo.x+kx+KMAX.x;
      //ne pas depasser l'ouverture numérique pour 1 hologramme
      if(kx*kx+ky_carre<NXMAX_CARRE)
	{
	  double kz_carre=r2-kx*kx-ky_carre; //altitude au carré des données
	  kz[cpt]=sqrt(kz_carre)-kz0;
	  //kz[cpt]=sqrt(kz_carre);
	}
      else  {

	kz[cpt]=0;
      }
    } //fin for y
  }
  //SAV(kz, NbPixHolo, "/home/hui/maniptomo/IDP/SA/kz.bin", FLOAT,"a+b");
  ///---------------------CALCULER LE SPECTRE REPROPAGE-----------------------------------
  double energie=0;
  for(int pix=0;pix<NbPixHolo;pix++)
    {
      rephase[pix].Re=cos(kz[pix]*deltaZ);
      rephase[pix].Im=sin(kz[pix]*deltaZ);

      spectrePropag[pix].Re=fft_shift_normRe[pix]*rephase[pix].Re-fft_shift_normIm[pix]*rephase[pix].Im;
      spectrePropag[pix].Im=fft_shift_normRe[pix]*rephase[pix].Im+fft_shift_normIm[pix]*rephase[pix].Re;
      // energie=sqrt(pow(spectrePropag[pix].Re,2)+pow(spectrePropag[pix].Im,2))+energie;
    }
  // cout<<"energie moyenne par pixel="<<energie/taillePlan<<endl;
  //SAV_Re(rephase, 4*KMAX.x*KMAX.x, "/home/mat/tomo_test/rephase.bin", FLOAT,"a+b");
  //SAV_Re(spectrePropag, NbPixHolo, "/home/mat/tomo_test/spectrePropag.bin", FLOAT,"a+b");
  //decal2DCplxGen(spectrePropag,spectrePropag_I,dimSpctHolo,KMAX);
  //SAV_Re(spectrePropag, 4*KMAX.x*KMAX.x, "/home/mat/tomo_test/spectrePropag.bin", FLOAT,"a+b");

  ///--------Renormalisation de la phase du fond après repropagation
  ///le spéculaire étant centré, on indique les coordonénes du milieu de l'image

  //cout<<"kmax="<<KMAX.x<<endl;
  int cptSpec=(dimSpctHolo.x)*CoordSpecI.y+CoordSpecI.x;
  //cout<<"cptSpec="<<cptSpec<<endl;
  //cout<<"valeur="<<spectrePropag[cptSpec].Re<<endl;

  double max_module=pow(spectrePropag[cptSpec].Re,2)+pow(spectrePropag[cptSpec].Im,2),
    max_part_reel=spectrePropag[cptSpec].Re,max_part_imag=spectrePropag[cptSpec].Im;
  for(int cpt=0; cpt<NbPixHolo; cpt++) {
    spectrePropag_normRe[cpt]=(spectrePropag[cpt].Re*max_part_reel+spectrePropag[cpt].Im*max_part_imag)/max_module;
    spectrePropag_normIm[cpt]=(spectrePropag[cpt].Im*max_part_reel-spectrePropag[cpt].Re*max_part_imag)/max_module;
    //spectrePropag_norm[cpt].Re=spectrePropag[cpt].Re;
    //spectrePropag_norm[cpt].Im=spectrePropag[cpt].Im;
  }

  // SAV_Re(spectrePropag_norm, NbPixHolo, "/home/mat/tomo_test/spectrePropag_Norm.bin", FLOAT,"a+b");
}



// *****************************************************************************
//
// *****************************************************************************


void 
Grad(double* entreeRe, double* entreeIm, double* Gradient, int taille, int NbPlan, int cptPlan)
{
  ///calcul amplitude
  double *ampli=new double[taille*taille];
  double *ampliGrad=new double[taille*taille];
  for(int pixel=0; pixel<taille*taille; pixel++)
    {
      ampli[pixel]=sqrt(pow(entreeRe[pixel],2)+pow(entreeIm[pixel],2));
    }
  //SAV(ampli, taille*taille, "/home/hui/maniptomo/IDP/SA/ampli3D.bin", FLOAT,"a+b");

  double *phase=new double[taille*taille];
  double *phaseGrad=new double[taille*taille];

  ///calcul phase de -PI à PI
  for(int pixel=0;pixel<taille*taille;pixel++)
    {
      double denominPhase=entreeRe[pixel];
      double numeraPhase=entreeIm[pixel];

      if(denominPhase>0)
        {
	  double racine=sqrt(pow(numeraPhase,2)+pow(denominPhase,2));
	  phase[pixel]=asin(numeraPhase/racine);
        }
      else
        {
	  double racine=sqrt(pow(numeraPhase,2)+pow(denominPhase,2));
	  phase[pixel]=M_PI-asin(numeraPhase/racine);
        }
    }
  //SAV(phase, taille*taille, "/home/hui/maniptomo/IDP/champ2d/phase3D.bin", FLOAT,"a+b");

  ///Copier double vers matrice de opencv///
  Mat img(taille,taille,CV_64F);
  for(int y=0; y<taille; y++)
    {
      for(int x=0; x<taille; x++)
	{
	  int cpt=y*taille+x;
	  img.at<double>(y,x)=ampli[cpt];
	}
    }

  /// sobel
  int scale = 1;
  int delta = 0;
  int ddepth = CV_64F;

  Mat grad_x, grad_y, imgResult;

  /// Gradient X
  //Scharr( img1[0], grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
  Sobel( img, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
  /// Gradient Y
  //Scharr( img1[0], grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
  Sobel( img, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
  /// Total Gradient (approximate)
  addWeighted( grad_x, 0.5, grad_y, 0.5, 0, imgResult );

  //    /// laplacian
  //    Mat imgResult;
  //    int kernel_size = 3;
  //    int scale = 1;
  //    int delta = 0;
  //    int ddepth = CV_64F;
  //    Laplacian(img, imgResult, ddepth, kernel_size, scale, delta, BORDER_DEFAULT );

  ///Copier matrice de opencv vers double///
  for(int y=0; y<imgResult.rows; y++)
    {
      for(int x=0; x<imgResult.cols; x++)
	{
	  int cpt=y*imgResult.cols+x;
	  ampliGrad[cpt]=imgResult.at<double>(y,x);
	}
    }
  //SAV(ampliGrad, taille*taille, "/home/hui/maniptomo/IDP/SA/ampli3D_Grad.bin", FLOAT,"a+b");

  double *ampli1= new double[taille*taille];
  double *ampli2= new double[taille*taille];
  Var2D decal1={0,-1};
  Var2D decal2={-1,0};
  Var2D dimFin={taille, taille};

  circshift2(ampliGrad,ampli1,dimFin,decal1);
  circshift2(ampliGrad,ampli2,dimFin,decal2);

  double *Grad=new double[taille*taille];
  double sommegrad=0;
  for(int cpt=0;cpt<taille*taille;cpt++)
    {
      Grad[cpt]=sqrt(pow(ampliGrad[cpt]-ampli1[cpt],2)+pow(ampliGrad[cpt]-ampli2[cpt],2));
      sommegrad=sommegrad+Grad[cpt];
    }
  //SAV(Grad, taille*taille, "/home/hui/maniptomo/IDP/champ_2d/Grad.bin", FLOAT,"a+b");

  double *sommeGrad=new double[NbPlan];
  for(int cpt=0;cpt<NbPlan;cpt++)
    {
      sommeGrad[cpt]=sommegrad;
    }
  Gradient[cptPlan]=sommeGrad[cptPlan];
  //cout<<"cptPlan= "<<cptPlan+1<<" sommegradient="<<sommegrad<<endl;
}




// =============================================================================
// fonctions outil 
// =============================================================================


// *****************************************************************************
//
// *****************************************************************************


void circshift2(double* entree, double* result, Var2D dim,Var2D decal)
{

  //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
  decal.y=decal.y%dim.y;
  decal.x=decal.x%dim.x;

  for(int yi=0;yi<dim.y;yi++)
    {

      int yi_decal=yi+decal.y;
      //si dépassement des dimensions
      if(yi+decal.y>dim.y-1)
	yi_decal=-dim.y+yi+decal.y;
      //si décage négatif
      if(yi+decal.y<0)
	yi_decal=dim.y+yi+decal.y;
      int hauteur_decal=yi_decal*dim.x;
      int hauteur=yi*dim.x;
      for(int xi=0;xi<dim.x;xi++)
	{
	  int xi_decal=xi+decal.x;

	  if(xi+decal.x>dim.x-1)
	    xi_decal=-dim.x+xi+decal.x;
	  if(xi+decal.x<0)
	    xi_decal=dim.x+xi+decal.x;

	  int pos1D_init=hauteur+xi;
	  int pos1D_shift=hauteur_decal+xi_decal;//nouvelle position après décalage
	  result[pos1D_shift]=entree[pos1D_init];

	}
    }
}



// *****************************************************************************
//
// *****************************************************************************


void TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y)
{
  fftw_plan_with_nthreads(6);

  int N=taille_x*taille_y;
  fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
  fftw_plan p;
  //Réservation memoire
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  //Récupération de l'image dans la partie reelle de l'entree
  for(int cpt=0;cpt<N;cpt++)
    {
      in[cpt][0]=entree_reelle[cpt];
      in[cpt][1]=entree_imag[cpt];
    }
  //calcul du plan, parametre servant a calculer et optimiser le FFT
  p=fftw_plan_dft_2d( taille_x,  taille_y, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */

  for(int cpt=0;cpt<(N);cpt++)
    {

      fft_reel[cpt]=out[cpt][0]/(N*N); //division par N^2 pour normaliser la fftw qui n'est pas normalisée
      fft_imag[cpt]=out[cpt][1]/(N*N);
    }
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}


// *****************************************************************************
//
// *****************************************************************************


void TF2D_INV(double entree_reelle[],double entree_imag[],double sortie_reelle[],double sortie_imag[],int taille_x,int taille_y)
{
  fftw_plan_with_nthreads(6);
  int N=taille_x*taille_y;
  fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
  fftw_plan p;
  //Réservation memoire
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  //Récupération de l'image dans la partie reelle de l'entree
  for(int cpt=0;cpt<N;cpt++)
    {
      in[cpt][0]=entree_reelle[cpt];
      in[cpt][1]=entree_imag[cpt];
    }
  //calcul du plan, parametre servant a calculer et optimiser le FFT
  p=fftw_plan_dft_2d(taille_x,  taille_y, in, out,FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */

  for(int cpt=0;cpt<(N);cpt++)
    {

      sortie_reelle[cpt]=out[cpt][0];
      sortie_imag[cpt]=out[cpt][1];
    }
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}


// *****************************************************************************
//
// *****************************************************************************


void Plan_ds_Vol(double *Vol3DRe, double *Vol3DIm, double *plan2DRe, double *plan2DIm, Var3D dimVol, Var2D dimPlan, int z3Di)
{
  int x3Di = 0, y3Di = 0;
  int taillePlan = dimPlan.x * dimPlan.y;
  int planVol3D = dimVol.x * dimVol.y;
  int ecartX = (dimVol.x - dimPlan.x) /2, ecartY = (dimVol.y - dimPlan.y) /2;
  int altitude_k = (z3Di +round(dimVol.z / 2)) * planVol3D;
  int cpt2D = 0, k = 0;

  for (int yi = 0; yi < dimPlan.y; yi++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
    for (int xi = 0; xi < dimPlan.x; xi++) { //on balaye l'image 2D en y, centre au milieu
      cpt2D = yi * dimPlan.x + xi;//calcul du cpt du tableau 1D de l'image 2D
      x3Di = xi + ecartX;
      y3Di = yi + ecartY;
      k = altitude_k + y3Di * dimVol.x + x3Di;
      Vol3DRe[k] = plan2DRe[cpt2D];
      Vol3DIm[k] = plan2DIm[cpt2D];
    }
  }
}


// *****************************************************************************
//
// *****************************************************************************


void read_dataM(char *inputfile,Mat* Data, int image_size, int IMAGE_NUM)
{
  //printf("Reading the Data Values form Binary File.............>");
  FILE *ifptr;
  float precision;
  Mat Data_tmp(image_size,image_size, CV_64F);
  ifptr = fopen(inputfile,"rb");
  if(ifptr == NULL) 
    { MSG_ASSERT(FALSE, "file not found"); }

  for(int z=0;z<IMAGE_NUM;z++)
    {
      for(int i=0;i<image_size;i++)
        {
	  for(int j=0; j< image_size; j++)
            {
	      fread(&precision,sizeof(precision),1,ifptr);
	      Data_tmp.at<double>(i,j)=precision;
            }
        }
      //cout<<"precision:"<<z<<endl;
      Data_tmp.copyTo(Data[z]);
    }
  fclose(ifptr);
  //printf(" Done.\n");
}
