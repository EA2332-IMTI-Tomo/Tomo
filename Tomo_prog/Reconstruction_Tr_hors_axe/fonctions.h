#ifndef __FONCTIONS__
#define __FONCTIONS__

/* --------------------------------------------------------------------------- */
// incl
/* --------------------------------------------------------------------------- */
#include "struct.h"
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <strstream>
//#include <Magick++.h>
#include <fftw3.h>

#include <cstring>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <tiffio.h>

//using namespace Magick;
using namespace std;

#include "projet.h"
typedef struct{
clock_t init, fin;
float   total;
}temps;

/* --------------------------------------------------------------------------- */
// Prototypes
/* --------------------------------------------------------------------------- */
void antigaussienne(double *tab, int Tx, int sigma, float A, int Exy);

void divCplx(nbCplx* imageA,nbCplx* imageB,nbCplx* resultat, int NbPix);
void multiplierCplx(nbCplx* image1,nbCplx* image2,nbCplx* resultat, int NbPix);
void decal2DGen(double* entree, double* result, Var2D dim,Var2D decal);
Var2D corr_crois(double* obj2D_A, double* obj2D_B, Var2D dim);

void methodeCarre(int NbPixROI2d, double *holo1,  double *holo2,  double *holo3,  double *holo4);
void changeDim2D(double* tab, double* tabFinal, Var2D dimInit,Var2D dimFin);
void changeDim2DCplx(nbCplx *tab, nbCplx* tabFinal, Var2D dimInit, Var2D dimFin);
//double *circshift(double *entree,int dimx,int dimy,int decal_x,int decal_y);
double bruit(int attenuation);
void calcPhase2pi(nbCplx* obj, Var2D taille,double* phaseMod2pi);
void calc_Uborn(nbCplx *TF_UBorn,nbCplx *UBorn,Var2D dim2DHA,Var2D PosSpec);
void Chrono(temps *t, string message);
void coupeCplx(nbCplx *src, nbCplx *dest, Var2D dim_src, Var2D dim_dest, Var2D coin);
void circshift2(double* entree, double* result, Var2D dim,Var2D decal);
void circshift3(double* entree, double* result, Var2D dim,Var2D decal);
void circshift2DCplx(nbCplx* entree, nbCplx* result, Var2D dim,Var2D decal);
void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D);
void circshift3DCplx(nbCplx *volume3D, nbCplx *volume3D_shift, Var3D dimFinal3D, Var3D decal3D);
void charger_image2D(unsigned char* phasei, string chemin,Var2D coin,Var2D taille);

int chargeBin(float *objet, string chemin,  int NbPix);

void circshift3D(double *volume3D, double *volume3D_shift,int taille_x,int taille_y,int taille_z);
void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D);
int CreerZoneFresnel(double *FresnelRe,double * FresnelIm, Var2D dim, Var2D centre, float d, float lambda);
string extract_string(std::string token,  std::string chemin_fic);
void genere_rectang3D(nbCplx *objet,Var3D posI_Coin,Var3D dimRect, Var3D dim);
void genere_rectang2D(double *objet,Var2D posI_Coin,Var2D dimRect,Var2D dim);
void genere_OTF_RB_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon);
void genere_OTF_T_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon);
void genere_OTF_RH_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon);
void holo2TF_UBorn(unsigned char *holo1, nbCplx *TF_UBorn,Var2D dimROI,Var2D decalROI, Var2D dim2DHA, Var2D coinHA);

void lire_bin(string chemin, double resultat[], short int precision, int NbParam);
//void rempli_tableau(unsigned char *finalArray, string path, Var2D coin, Var2D taille);
void TF2Dcplx(nbCplx *entree, nbCplx *fft, Var2D dim);
void TF2Dcplx_INV(nbCplx *fft_entree, nbCplx *objet, Var2D dim);
double *tukey2D(int dimx,int dimy, float alpha);

void genereCache(double masque[], int t_image, int t_mask, int centreX, int centreY);
void interp3D(double *volume_interp_3D, int taille_x,int taille_y,int taille_z);
void ecrire_rapport(int NXMAX,float rayon,float Rf, int DIMX_CCD2,int coin_x, int coin_y,short int precision_exportation,string chemin,int nb_proj,float n1,float NA,float Tp, int G);
//void fftw_plan_with_nthreads(int nthreads);
void multiplier_masque(double image[], unsigned char masque[], int t_image, int t_mask, int centreX, int centreY);
//void concatener(char chaine1[],char chaine2[],char resultat[]);

void multiplier_masque2(double image[], double masque[], int t_image, int t_mask, int centreX, int centreY);
void multiplier_masque2Cplx(nbCplx image[], double masque[], int t_image, int t_mask, Var2D posCentre);
void multiplier_masqueCplx2(nbCplx *image, nbCplx *masque, int t_image, int t_mask, Var2D CentreI);
void prepare_wisdom2D(Var2D dim, const char *chemin);
void prepare_wisdom3D(Var3D dim, char *chemin);
void SAV(double *var_sav, int NbPix2D, char* chemin, enum PRECISION,char options[]);
void SAV2(double *var_sav, int NbPix2D, string chemin, enum PRECISION precision, char options[]);

void SAV_Re2(nbCplx *var_sav, int taille, string chemin, enum PRECISION precision, char options[]);
void SAV_Im2(nbCplx *var_sav, int taille, string chemin, enum PRECISION precision, char options[]);
void SAV_Re(nbCplx *var_sav, int taille, char *chemin, enum PRECISION precision, char options[]);
void SAV_Im(nbCplx *var_sav, int taille, char *chemin, enum PRECISION precision, char options[]);

void SAV_Tiff3D(nbCplx *var_sav, string chemin_ind, string chemin_abs, int dim);
void SAV_Tiff2D(double *var_sav, string chemin, int dim);
void TF3DCplx(nbCplx *Objet3D_shift, nbCplx *Spect_shift,Var3D dimVol);
void TF3DCplx_INV(nbCplx *Spect3D_shift, nbCplx* Objet3D_shift,Var3D dimVol);

double max(double* entree, int dim);
void phase2pi(nbCplx* obj, Var2D taille,double* WrappedImage);
int retroPropag_Born(nbCplx *TF3D_PotObj, nbCplx *TF_Uborn_norm, double * sup_redon, int dim_final, Var2D posSpec, Var3D decal3D, Var2D NMAX, double rayon);
void retroPropagSA(int deltaZ, nbCplx *fft_shift_norm, nbCplx *planObjet, Var3D decal, Var2D NMAX, double rayon);
void decalCoupeCplx(nbCplx *fft, nbCplx *fft_tmp, Var2D NMAX,Var2D dimCCD);
void Plan_ds_VolCplx(nbCplx *Vol3D, nbCplx *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di);
void Plan_ds_Vol(double *Vol3D, double *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di);
void Vol_ds_Plan(double *Vol3D, double *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di);
int coordMaxMod2D(nbCplx *entree, int tailleTab);
void InitTabCplx(nbCplx *z,int taille);
void decal2DCplxGen(nbCplx* entree, nbCplx* result, Var2D dim,Var2D decal);
void recale(nbCplx* obj,nbCplx* objDecal,nbCplx *objRecal, Var3D dimVol);

///calcul sur cplx
void AXB_cplx(nbCplx *A, nbCplx *B, nbCplx *produit, int NbPix3D);///multiplication complexe
void conj_cplx(nbCplx *A, nbCplx *conj, int NbPix3D); ///conjugaison complexe
void mod_cplx(nbCplx *A, int NbPix3D);
void recal_obj(nbCplx *FFTobj, nbCplx *FFTobjDecal,nbCplx *objRecal, Var3D dimVol);
#endif

