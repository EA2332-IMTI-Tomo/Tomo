
/***************************************************************************
 *   Copyright (C) 2007 by mat   *
 *   mat@feynman   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


// #define CHEMIN "/home/bailleul/ProjetTOMO/session/session10032302-record"
#define CHEMIN "INPUT/session10032302-record"
#define MASKPATH "/home/bailleul/ProjetTOMO/MASQUES"
#define TEST_PATH "OUTPUT"



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <octave-2.9.9/octave/oct.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
//#include "/home/mat/Mat/Projet_C/ImageMagick-6.3.2/Magick++/lib/Magick++.h"
#include <Magick++.h>
//#include "/usr/include/ImageMagick/Magick++.h"
#include <fftw3.h>
#include <cstring>
#include <fstream>


using namespace std;
using namespace Magick;

double *circshift(double *entree,int dimx,int dimy,int decal_x,int decal_y);

unsigned char *rempli_tableau(string path, int coin_x, int coin_y,int taille_x,int taille_y);

void TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y);

unsigned char* charger_image2D(unsigned char* phasei,int numero, char* chemin,int cpt_fichier,int coin_x,int coin_y,int taille_x,int taille_y);

void circshift3D(double *volume3D, double *volume3D_shift,int taille_x,int taille_y,int taille_z);
//void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D)

void interp3D(double *volume_interp_3D, int taille_x,int taille_y,int taille_z);
void ecrire_rapport(int Nxmax,float rayon,float Rf, float K, int dimx_ccd2,int coin_x, int coin_y,double precision,char *chemin,int nb_proj,float n1,float NA,float Tp, int G);
//void fftw_plan_with_nthreads(int nthreads);

int main(int argc, char *argv[])
{
	/*pointeur pour ouverture de fichier
	FILE* fichier_TF_real = NULL;
	//ouverture de ce fichier en écriture binaire
	fichier_TF_real = fopen("/home/mat/tomo_test/test_resultat2.bin", "wb");*/
	//int dimx_ccd2=512,  dimy_ccd2=512;
	int dimx_ccd2=512,  dimy_ccd2=512;

	//valeur controlant l'exclusion de certains centres
	int xm0_limite=200; //centre maximum
	int ym0_limite=200;
	int rayon_inf=xm0_limite*xm0_limite+ym0_limite*ym0_limite;
	//valeur du coin pour la découpe
	int coin_x=50, coin_y=50;

if(coin_x+dimx_ccd2>740 || coin_y+dimy_ccd2>574)
printf("Attention, la zone decoupee sort de l'image\n");
//largeur : 740, hauteur : 574.

//short int precision;
float precision;
int jumeau_elimine=0;
int points_faux=0;
/*int fftwThreadInit;
 fftwThreadInit=fftw_init_threads();
int nthreads=2;
printf("fftwThreadInit: %i\n",fftwThreadInit);*/
//Récupération du chemin en tant que 1 er argument de la ligne de commande
/*	for (int i=0; i < argc; i++)
	{
		printf("Argument %li : %s \n", i+1, argv[i]);
	}
	 //getchar();//attendre un appui sur une touche
	long longueurChemin = 0;
	// On récupère la longueur de la chaîne dans longueur_chemin
	longueurChemin = strlen(argv[1]);
	char Chemin[longueurChemin];
	strcpy(Chemin, argv[1]);
	printf("%s\n",Chemin);

	if(argc<1)
	{
	char Chemin[]="/home/mat/Mat/Matlab/A_Faire/session07022601_diatom_tete2mort/session07022601-record";
	}
	char NomSupport[]="/home/mat/Mat/tranche_maper3D/support_redondance.bin";*/
//char Chemin[]="/opt/resultats2009/session09010701/session09010701-record";
//char Chemin[]="/opt/resultats2009/Reflexion/session09042704/session09042704-record";
//char Chemin[]="/opt2/resultats2009/session09010903/session09010903-record";
//char Chemin[]="/opt2/resultats2009/LYON/session09091101/session09091101-record";
//opt2/resultats2009/session09030402_virus_lyon2/session09030402-record";
//char Chemin[]="/home/mat/Mat/Matlab/A_Faire/session07100201-reglet14/session071002-record";
//char NomSupport[]="/home/mat/Mat/tranche_maper3D/support_redondance.bin";

 char Chemin[] = CHEMIN;

// ** char Chemin[]="/opt2/resultats2010/transmission/modulateur/session10032601_reglet_1.31v_800a/session10032601-record";// /opt2/resultats2009/lyon_aout_09/muscles/session09082401-record";
/*/opt2/resultats2009/LYON/session09111701
int coin_x_tmp;
int coin_y_tmp;*/
int premier_plan=1;
int Num_Angle_final=799;//
// int sans_signal=0;
int centres_exclus=0;
int NbAngle=Num_Angle_final-premier_plan+1;
int SautAngle=1;
//int cpt_angle=6;
	printf("coin x (default : %i) : ",coin_x);
	scanf("%i", &coin_x);
//if(coin_x_tmp==NULL)
//printf("coin_x_tmp : %i\n",coin_x_tmp);

	printf("coin y (default : %i): ",coin_y);

	scanf("%i", &coin_y);



//nombre de pixel (pour le tableau)
int N=dimx_ccd2*dimy_ccd2;

//---------coordonnees du coin haut gauche dans l'image originale (748x570) pour l'image recadrée --------



float n1=1.515;	//indice de l'huile
float NA=1.40;	//ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
float theta=asin(NA/n1);//theta_max
printf("theta : %f\n",theta);

float lambda=632*pow(10,-9); //longueur d'onde
int G=100;	//grossissement telan+objectif
float Tp=11.2*pow(10,-6); //Taille des pixels
float Rf=0.666;//1.2;
//float Rf=0.8;//1.7;//0.77;//1.3;// bonne valeur le 29/9; //rapport f1/f2


// Facteur d'echelle
float K=lambda*G/(2*NA*Tp*Rf);
printf("K : %f\n", K );
float Tps=Tp;//*K; //Taille pixel apres mise à l'echelle important si RF<>1

int Ntx=dimx_ccd2;//K; //Nombre total de pixel en x
/*int Nty=Ntx;*/ //(support carré)

float rayon=n1*Rf*Tps*Ntx/(G*lambda); //Rayon

float delta_zmax=rayon*cos(theta);//float?
//int rayon=round(n1*Rf*Tp*Ntx/(G*lambda));ici un commentaire archéologique
//int rayon=100;

//dimension du support dans Fourier
int Nxmax=round(n1*Rf*Tps*Ntx/(G*lambda)*NA/n1);//Nxmax=round(rayon*NA/n1)=round(rayon*sin theta_max);
int Nxmax_Rf=round(n1*Tps*Ntx/(G*lambda)*NA/n1);

float test_nxmax_Rf=n1*Tps*Ntx/(G*lambda)*NA/n1;
printf("testnxmax_Rf %f \n",test_nxmax_Rf);
printf("round(testnxmax_Rf) %f \n",round(test_nxmax_Rf));
int Nymax=Nxmax;
int nb_proj=0;




printf("Rayon %f \n",rayon);
//printf("Ntx %i \n",Ntx);
printf("Nxmax=%i,4*Nxmax=%i \n",Nxmax,4*Nxmax);
//Indice crée pour balayer le tableau final 1D (informatique) contenant le volume final de données en 3D
	int N_tab=64*Nxmax_Rf*Nxmax_Rf*Nxmax_Rf;
	//int N_tab=round(64*(n1*Tps*Ntx/(G*lambda)*NA/n1)*(n1*Tps*Ntx/(G*lambda)*NA/n1)*(n1*Tps*Ntx/(G*lambda)*NA/n1));//sans arrondis

printf("N_tab (indice espace 3D)= %i \n",N_tab);
	cout <<"Espace final : " << 4*round(n1*Tps*Ntx/(G*lambda)*NA/n1) << " pixels cube\n";

//printf("Précision d'exportation: %i bits\n chemin :%s \n",8*sizeof(precision),Chemin);

///---------------Calcul de quelques constantes, afin d'éviter de les calculer dans la boucle--------------
//dimx_arc
		const int dv0=4*Nxmax; //dimmension totale en x,y

		printf("dimx espace 3D, dv0: %i\n",dv0);
	//dimx_arc/2
// 		const int dv0s2=2*Nxmax;
// 		const int dv1s2=2*Nymax;

	//dimx_arc*dimy_arc
// 		const int dv1xdv0=16*Nxmax*Nymax;

	//dimz_arc
// 		const int dv2s2=2*Nxmax;
///////////////////////////////////////////////////////Rf

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

/// //////////////////////////////////////////////////fin calcul constante//////////////////////////


//reservation de 4 tableaux (image 2D)  pour le phase shifting
unsigned char* phase1=new unsigned char[N];
unsigned char* phase2=new unsigned char[N];
unsigned char* phase3=new unsigned char[N];
unsigned char* phase4=new unsigned char[N];

/// 4 tableaux pour transtypage des images 8 bits vers 64 bits (double).
double* phased1=new double[N];
double* phased2=new double[N];
double* phased3=new double[N];
double* phased4=new double[N];


//reservation reference et "noir" camera pour une "tentative" de correction de la ref
	unsigned char* reference=new unsigned char[N];
// 	unsigned char* noir_camera=new unsigned char[N];

//masque hamming
	unsigned char* masque=new unsigned char[N];

//cache objet jumeau
	unsigned char* cache_jumeau=new unsigned char[100];



//variables utilisées pour la TF2D

	//Avant Crop à NxMax;
	double *fft_reel_tmp=new double[N];
	double *fft_imag_tmp=new double[N];
	//Après Crop à NxMax;
	double *fft_reel=new double[4*Nxmax*Nymax];
	double *fft_imag=new double[4*Nxmax*Nymax];

//double *qualite_phase_shift=new double[N];
//espace3D reel, imaginaire et support de redondance

double *reel_arc=new double[N_tab];//partie reelle du papillon dans l'espace reciproque
double *imag_arc=new double[N_tab];//partie imaginaire du papillon dans l'espace reciproque
double *sup_redon=new double[N_tab];//pour normaliser par rapport au nombre de fois où l'on remplit la frequence
int *centre=new int[4*Nxmax*Nymax];//pour mettre la position des centres translatés, on crée une variable 2D de la taille d'un plan apres tomo

clock_t
temps_depart, /*temps de départ de programme */
temps_initial, /* temps initial en micro-secondes */
temps_final, /* temps d'arrivée du programme */
temps_arrivee;   /* temps final en micro-secondes */
float temps_cpu=0,     /* temps total en secondes */
temps_total=0;
temps_depart = clock();
temps_initial = clock ();

// int cpt_ref=0;
/*
reference=rempli_tableau("/home/mat/Mat/Matlab/A_Faire/session07101001/ref2.bmp", coin_x, coin_y,dimx_ccd2,dimy_ccd2);
noir_camera=rempli_tableau("/home/mat/Mat/Matlab/A_Faire/session07101001/fond_camera2.bmp", coin_x, coin_y,dimx_ccd2,dimy_ccd2);


//si un pixel est à zero dans reference, mettre 1
for(int pixel=0;pixel<N;pixel++)
{

	if(reference[pixel]<=noir_camera[pixel])
	{
	reference[pixel]=255;
	cpt_ref++;
	}
	else
	{
	reference[pixel]=reference[pixel]-noir_camera[pixel];
	}
}
double ref_tampon;

	FILE *fichier_reference;
		fichier_reference= fopen("/home/mat/tomo_test/ref_tampon.bin", "wb");
		for(int cpt=0;cpt<(dimx_ccd2*dimy_ccd2);cpt++)
			{
			//Ecriture sur le disque
			ref_tampon=(double)reference[cpt];
			fwrite(&ref_tampon,sizeof(double),1,fichier_reference);
			}
		fclose(fichier_reference);


printf("valeurs pixels dans fond camera> valeur ref dans ref : %i\n",cpt_ref);
*/
if(dimx_ccd2==512)
{
masque=rempli_tableau( MASKPATH "/kmasq_512x512_30.bmp", 0, 0,dimx_ccd2,dimy_ccd2);
printf("masque 512x512\n");
}
if(dimx_ccd2==256)
{
masque=rempli_tableau( MASKPATH "/kmasq_256x256.bmp", 0, 0,dimx_ccd2,dimy_ccd2);
printf("masque 256x256\n");
}

//taille totale du masque en pixel pour elimination objet jumeau;
const int t_mask=30;
cache_jumeau=rempli_tableau( MASKPATH "/k_30x30.bmp", 0, 0,t_mask,t_mask);
//printf("masque[128*128]= %d", masque[128*256+128]);
/*FILE * I_masque;
	I_masque = fopen("/home/mat/tomo_test/masque.bin", "wb");
	for(int cpt=0;cpt<N;cpt++)
		{
		//Ecriture sur le disque
		fwrite(&cache_jumeau[cpt],sizeof(unsigned char),1,I_masque);
		}
	fclose(I_masque);*/

printf("*******************************************\n");
printf("remplissage de l'espace réciproque\n");
printf("*******************************************\n");

printf("  \\,`//\n");
printf(" _).. `_\n");
printf("( __  -\\ \n");
printf("    '`.\n");
printf("   ( \\>\n");
printf("   _||_ \n");

FILE* test_existence;//nom bidon pour tester l'existence des fichiers

//début de la grosse boucle sur les angles/////////////////////////
//##############################################################################################################################

for(int cpt_angle=premier_plan;cpt_angle<premier_plan+NbAngle;cpt_angle=cpt_angle+SautAngle)//boucle sur tous les angles
{
//cout<<"cpt_angle : "<<cpt_angle<<endl;
	if((cpt_angle-100*(cpt_angle/100))==0)
	printf("cpt_angle=%i\n",cpt_angle);
	//*if (argc>1)
	//*tableau = rempli_tableau((string)*++argv);
	//*return 0;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////phase shifting////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//tester l'existence des fichiers au cas ou
	char test_Copie[100];
	char test_FinNom[15];
	//sauvegarder chemin car sprintf écrase
	strcpy(test_Copie,Chemin);
	sprintf(test_FinNom,"%i-001.bmp",cpt_angle);
	//créer le chemin final
	strcat(test_Copie,test_FinNom);
	//printf("test_Copie vaut : %s\n", test_Copie);
	test_existence = fopen(test_Copie, "rb");

if(test_existence!=NULL)
{
		fclose(test_existence);
		phase1=charger_image2D(phase1,1,Chemin,cpt_angle, coin_x, coin_y,dimx_ccd2,dimy_ccd2);
	/*	if(phase1==NULL)
		{
		printf("Le fichier n'existe pas\n");
		}
else{*/
		phase2=charger_image2D(phase2,2,Chemin,cpt_angle, coin_x, coin_y,dimx_ccd2,dimy_ccd2);
		phase3=charger_image2D(phase3,3,Chemin,cpt_angle, coin_x, coin_y,dimx_ccd2,dimy_ccd2);
		phase4=charger_image2D(phase4,4,Chemin,cpt_angle, coin_x, coin_y,dimx_ccd2,dimy_ccd2);



//Calcul du déphasage réel (méthode de Carré)
// double *tan_2delta_phi=new double[N];
// double somme_tan_2delta_phi=0;
// int nb_tan=0;
//
// //calcul de la valeur du saut de phase
//
// for(int pixel=0;pixel<N;pixel++)
// 		{
// 			phased1[pixel]=(double)phase1[pixel];
// 			phased2[pixel]=(double)phase2[pixel];
// 			phased3[pixel]=(double)phase3[pixel];
// 			phased4[pixel]=(double)phase4[pixel];
//
//  			tan_2delta_phi[pixel]=sqrt(fabs((3*(phased2[pixel]-phased3[pixel])-(phased1[pixel]-phased4[pixel]))/((phased2[pixel]-phased3[pixel])+(phased1[pixel]-phased4[pixel]))));
// 			if(tan_2delta_phi[pixel]<10)//la valeurde la tangente doit valoir environ tan(0.5*pi/2)=1, on élimine les valeurs aberrantes>10
// 			{
//  			somme_tan_2delta_phi=tan_2delta_phi[pixel]+somme_tan_2delta_phi;
// 			nb_tan++;
// 			}
// 		}
//
// double moyenne_tan_2delta_ph=somme_tan_2delta_phi/N;
// double somme_ecart_carre=0;
//
// //calcul écart type de la valeur du saut de phase
// 	for(int pixel=0;pixel<N;pixel++)
// 	{
// 		if(tan_2delta_phi[pixel]<10)//la valeur doit valoir environ 1, on élimine les valeurs aberrantes>10
// 		{
// 		somme_ecart_carre=somme_ecart_carre+(tan_2delta_phi[pixel]-moyenne_tan_2delta_ph)*(tan_2delta_phi[pixel]-moyenne_tan_2delta_ph);
// 		}
// 	}
//
// double ecart_type=sqrt(somme_ecart_carre/N);
// //affichage saut de phase
// printf("la moyenne du phase shifting sur cette image vaut: %f avec nb_tan :%f et ecart type : %f\n",moyenne_tan_2delta_ph,sqrt(nb_tan),ecart_type);
// 	//allocation résultat du phase shifting attention à ne pas désallouer hors de la boucle for


	double *plan_reel=new double[N];
	double *plan_imag=new double[N];

		for(int pixel=0;pixel<N;pixel++)
		{
			plan_reel[pixel]=((double)phase1[pixel]-(double)phase3[pixel])*(masque[pixel]);
			plan_imag[pixel]=((double)phase4[pixel]-(double)phase2[pixel])*(masque[pixel]);
	//plan_reel[pixel]=(((double)phase1[pixel]/(double)reference[pixel])-((double)phase3[pixel]/(double)reference[pixel]));//*(masque[pixel]);
	//plan_imag[pixel]=(((double)phase4[pixel]/(double)reference[pixel])-((double)phase2[pixel]/(double)reference[pixel]));//*(masque[pixel]);
			//qualite_phase_shift[pixel]=(double)phase4[pixel]-(double)phase3[pixel]-(double)phase1[pixel]-(double)phase2[pixel
		}

/*
					//Sauver les plans reels
		char Nom_PS[100]="/home/mat/tomo_test/plan_PS_2D/Plan_reel_2D_";
		char Num_PS[10];
		sprintf(Num_PS,"%i.bin",cpt_angle); //cpt_angle->Num_PS
		strcat(Nom_PS,Num_PS);
		FILE *fichier_ps;
		fichier_ps= fopen(Nom_PS, "wb");
		for(int cpt=0;cpt<(dimx_ccd2*dimy_ccd2);cpt++)
			{
			//Ecriture sur le disque
			fwrite(&plan_reel[cpt],sizeof(double),1,fichier_ps);
			}
		fclose(fichier_ps);

					//sauver les plan imag

	char Nom_PS_imag[100]="/home/mat/tomo_test/plan_PS_2D/Plan_imag_2D_";
		char Num_PS_imag[10];
		sprintf(Num_PS_imag,"%i.bin",cpt_angle); //cpt_angle->Num_PS
		strcat(Nom_PS_imag,Num_PS_imag);
		FILE *fichier_ps_imag;
		fichier_ps_imag= fopen(Nom_PS_imag, "wb");
		for(int cpt=0;cpt<(dimx_ccd2*dimy_ccd2);cpt++)
			{
			//Ecriture sur le disque
			fwrite(&plan_imag[cpt],sizeof(double),1,fichier_ps_imag);
			}
		fclose(fichier_ps_imag);
	*/
		/*FILE * I_plan_reel;
		I_plan_reel = fopen("/home/mat/tomo_test/plan_reel_C.bin", "wb");
		for(int cpt=0;cpt<N;cpt++)
			{
			//Ecriture sur le disque
			fwrite(&plan_reel[cpt],sizeof(double),1,I_plan_reel);
			}
		fclose(I_plan_reel);
		*/

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////// Circshift avant TF2D
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//int N=dimx_ccd2*dimy_ccd2;
		double  *plan_reel_shift=new double[N];
		double  *plan_imag_shift=new double[N];


		plan_reel_shift=circshift(plan_reel, dimx_ccd2,dimy_ccd2,dimx_ccd2/2,dimy_ccd2/2);
		plan_imag_shift=circshift(plan_imag,dimx_ccd2,dimy_ccd2,dimx_ccd2/2,dimy_ccd2/2);




		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////TF 2D du front d'onde: passage du plan image au plan réciproque////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		TF2D(plan_reel_shift,plan_imag_shift,fft_reel_tmp,fft_imag_tmp,dimx_ccd2,dimy_ccd2);



		delete[] plan_reel_shift;
		delete[] plan_imag_shift;
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////Découpage dans fourier à Nxmax (crop dans fourier = mise à l'échelle dans plan image)
		////////////////////////////la frequence max correspond maintenant a l'ouverture numerique experimentale//////////////////////////////////
		////////////////////////////fft_reel_tmp=entree, fft_reel=sortie//////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		int xi;
		int yi;
		///////////////////////////////////////////////premier demi-espace
		for (xi=0;xi<Nxmax;xi++)
		{
			for ( yi=0;yi<Nymax;yi++)
			{
				int cpt1=yi*dimx_ccd2+xi;
				int cpt2=yi*(Nxmax)*2+xi;

				fft_reel[cpt2]=fft_reel_tmp[cpt1];
				fft_imag[cpt2]=fft_imag_tmp[cpt1];
			}

			for (int yi=dimy_ccd2-Nymax;yi<dimy_ccd2;yi++)
			{
				int cpt1=yi*dimx_ccd2+xi;

				int cpt2 = xi+(yi-dimy_ccd2+2*Nymax)*2*(Nxmax);
				fft_reel[cpt2]=fft_reel_tmp[cpt1];
				fft_imag[cpt2]=fft_imag_tmp[cpt1];
			}
		}

		///////////////////////////////////////////////deuxieme demi-espace
		for(int xi=dimx_ccd2-Nxmax;xi<dimx_ccd2;xi++)
		{
				for (int yi=0;yi<Nymax;yi++)
			{
				int cpt1=yi*dimx_ccd2+xi;
				int cpt2=yi*(Nxmax)*2+(xi-dimx_ccd2+2*Nxmax);

				fft_reel[cpt2]=fft_reel_tmp[cpt1];
				fft_imag[cpt2]=fft_imag_tmp[cpt1];
			}

			for (int yi=dimy_ccd2-Nymax;yi<dimy_ccd2;yi++)
			{
				int cpt1=yi*dimx_ccd2+xi;

				int cpt2 = (xi-dimx_ccd2+2*Nxmax)+(yi-dimy_ccd2+2*Nymax)*2*(Nxmax);
				fft_reel[cpt2]=fft_reel_tmp[cpt1];
				fft_imag[cpt2]=fft_imag_tmp[cpt1];
			}
		}

		//printf("  Nxmax : %i, Nymax : %i\n",Nxmax,Nymax);

		/*FILE* TF_real = NULL;
		TF_real = fopen("/home/mat/tomo_test/TF_reel.bin", "wb");

		for(int cpt=0;cpt<(4*Nxmax*Nymax);cpt++)
			{
			//Ecriture sur le disque
			fwrite(&fft_reel[cpt],sizeof(double),1,TF_real);
			}
		fclose(TF_real);*/


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////circshift 2D apres TF
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		double  *fft_reel_shift=new double[4*Nxmax*Nymax];
		double  *fft_imag_shift=new double[4*Nxmax*Nymax];

		fft_reel_shift=circshift(fft_reel, 2*Nxmax,2*Nymax,Nxmax,Nymax);
		fft_imag_shift=circshift(fft_imag, 2*Nxmax,2*Nymax,Nxmax,Nymax);


		//delete(fft_reel_shift);//?????????
			//Sauver les Tf2D
		/*char Nom_TF[100]="/home/mat/tomo_test/TF2d/TF2D_";
		char Num_TF[10];
		sprintf(Num_TF,"%i.bin",cpt_angle);
		strcat(Nom_TF,Num_TF);
		FILE *fichier_tf2d;
		fichier_tf2d= fopen(Nom_TF, "wb");
		for(int cpt=0;cpt<(4*Nxmax*Nymax);cpt++)
			{
			//Ecriture sur le disque
			fwrite(&fft_reel_shift[cpt],sizeof(double),1,fichier_tf2d);
			}
		fclose(fichier_tf2d);*/

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////On cherche la position du maximum d'intensité					///
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//----------------------------Recherche du maximum non centré-----------------------------------------


		//gloubiboulga function /!\ achtung, dont touch /!!\ Caramba, Né tochas
		double *fft_module_shift=new double[4*Nxmax*Nymax];
// 		double fft_mod_max=0;
// 		double fft_reel_max=0;

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




		/*FILE* fichier_TF_mod = NULL;
		fichier_TF_mod = fopen("/home/mat/tomo_test/fft_module_shift.bin", "wb");
			for(int cpt=0;cpt<(4*Nxmax*Nymax);cpt++)
				{
				//Ecriture sur le disque
				fwrite(&fft_module_shift[cpt],sizeof(double),1,fichier_TF_mod);
				}

			fclose(fichier_TF_mod);
		*/


		double max_part_reel = fft_reel_shift[cpt_max];
		double max_part_imag = fft_imag_shift[cpt_max];
		double max_module = fft_module_shift[cpt_max];

		//printf("max_part_reel : %f, max_part_imag: %f,  max_module %f \n ",max_part_reel, max_part_imag, max_module);


		/*FILE * TF_reel_shift;
		TF_reel_shift = fopen("/home/mat/tomo_test/test_resultat_shift.bin", "wb");
		for(int cpt=0;cpt<(4*Nxmax*Nymax);cpt++)
			{
			//Ecriture sur le disque
			fwrite(&fft_reel_shift[cpt],sizeof(double),1,TF_reel_shift);
			}
		fclose(TF_reel_shift);*/

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////normalisation par le pic central
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		double  *fft_reel_shift_norm=new double[4*Nxmax*Nymax];
		double  *fft_imag_shift_norm=new double[4*Nxmax*Nymax];


		for(int cpt=0;cpt<(4*Nxmax*Nymax);cpt++)
			{
			fft_reel_shift_norm[cpt]=(fft_reel_shift[cpt]*max_part_reel+fft_imag_shift[cpt]*max_part_imag)/max_module;
			//fft_reel[cpt]=(fft_reel[cpt]*fft_reel[cpt_max]+fft_imag[cpt]*fft_imag[cpt_max])/fft_module[cpt_max];danger!!!!
			fft_imag_shift_norm[cpt]=(fft_imag_shift[cpt]*max_part_reel-fft_reel_shift[cpt]*max_part_imag)/max_module;
			//fft_imag[cpt]=(fft_imag[cpt]*fft_reel[cpt_max]-fft_reel[cpt]*fft_imag[cpt_max])/fft_module[cpt_max];danger!!!!
			}

	delete[] fft_reel_shift;
	delete[] fft_imag_shift;
	delete[] fft_module_shift;

		//--------------------------------------


		//Coordonnées dans l'espace 2D à partir de l'indice 1D: (xc,yc)=(0,0)=en haut à gauche de l'image
		int xmi=cpt_max%(2*Nxmax);
		int ymi=cpt_max/(2*Nymax);

		int xc=xmi-Nxmax;
		int yc=ymi-Nymax;
		//coordonnée objet jumeaux
		int xmij=2*Nxmax-xmi;
		int ymij=2*Nymax-ymi;
	////////////////////////////////////Virer l'objet jumeaux
		if((xc*xc+yc*yc)>900)//35*35=1225;objet jumeau pas trop près de l'objet
		{ jumeau_elimine++;

			for(int x=xmij-t_mask/2;x<xmij+t_mask/2;x++)
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




		/*	FILE* fichier_fft_reel = NULL;
			//ouverture de ce fichier en écriture binaire
			char nom_fichier[100];
			sprintf(nom_fichier,"/home/mat/tomo_test/TF2d/fft_reel_shift-%ld.bin",cpt_angle);
			fichier_fft_reel = fopen(nom_fichier, "wb");
			for(int cpt=0;cpt<4*Nxmax*Nymax;cpt++)
			{
				precision=fft_reel_shift_norm[cpt];
				fwrite(&precision,sizeof(precision),1,fichier_fft_reel);
			}
			fclose(fichier_fft_reel);
		*/




		//printf("xmi: %i, ymi : %i\n",xmi,ymi);



		//printf("cpt_max vaut : %i ; xc : %i , yc : %i \n", cpt_max, xc, yc);

		//Conversion vers le centre du tableau
	//printf("coucou2\n");
		//delete(fft_module_shift);//une fois le max trouve, on peut liberer fft_module
		//delete[] fft_module_shift;// modif vl du 4/10/2007
		/*FILE * TF_imag_shift_norm;
		TF_imag_shift_norm = fopen("/home/mat/tomo_test/test_resultat_imag_shift_norm_C.bin", "wb");
			for(int cpt=0;cpt<(4*Nxmax*Nymax);cpt++)
				{
				//Ecriture sur le disque
				fwrite(&fft_imag_shift[cpt],sizeof(double),1,TF_imag_shift_norm);
				}
		fclose(TF_imag_shift_norm);*/


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////Mapping 3D
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		//coordonnée dans l'image2D centrée (xm0,ym0)=(0,0)=au centre de l'image
		int xm0=(xmi-Nxmax);//
//printf("xm0 %i \n", xm0);
		int ym0=(ymi-Nymax);

			if(xm0==0 && ym0==0)
			printf("(xm0,ym0)=(0,0) dans le plan %i\n", cpt_angle);
	//if((xm0*xm0+ym0*ym0)>rayon_inf || fabs(xm0)>xm0_limite || fabs(ym0)>ym0_limite)
	if((xm0*xm0+ym0*ym0)>rayon_inf)
	{
	//printf("xm0 : %i, ym0:%i\n",xm0,ym0);
	//printf("\nFichier %i sans signal\n", cpt_angle);
	centres_exclus++;
	}
	else
	{
		centre[xmi*2*Nxmax+ymi]=cpt_angle;//sauvegarde des centres en coordonnées non centrée; on met le numero d'angle
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

		double zm0_carre = rayon*rayon-xm0*xm0-ym0*ym0;
		if(round(zm0_carre)>-1)
		{

		zm0=sqrt(zm0_carre);
		//printf("zm0:%i\n",zm0);
		nb_proj++;

		/* Essai en calculant x et y a partir de cpt
			for(int cpt=0;cpt<4*Nxmax*Nymax;cpt++)
			{
			int x = (cpt%(2*Nxmax));
			int y = (cpt/(2*Nymax));
		}*/
	temps_initial = clock ();
			for (int x = -Nxmax; x < Nxmax; x++)//on balaye l'image 2D en x , origine (0,0) de l'image au milieu
			{
				for (int y = -Nymax; y < Nymax; y++)//on balaye l'image 2D en y, centre au milieu
				{

					int cpt=(y+Nymax)*2*Nxmax+x+Nxmax;//calcul du cpt du tableau 1D de l'image 2D
					//argument sous la racine calcul dans ARC_centre
					arg_z_arc=r2-x*x-y*y;

					if(arg_z_arc>delta_zmax*delta_zmax)//ne pas depasser l'ouverture numérique pour 1 hologramme
					{

						z_arc=round(sqrt(arg_z_arc)-zm0);//altitude au centre du volume
						k=round((-xm0+x+dv0s2rf)+(-ym0+y+dv1s2rf)*dv0rf+(z_arc+dv2s2rf)*dv1xdv0rf);//indice du tableau 1D du volume 3D
						sup_redon[k]+=1;//pour calculer le support

						reel_arc[k]+=fft_reel_shift_norm[cpt];//pour calculer l'image
						imag_arc[k]+=fft_imag_shift_norm[cpt];//pour calculer l'image

					}
					else
					points_faux++;
				} //fin for y
			}

		}//fin if zm0>-1
	}//fin else xm0_limite

			temps_final = clock ();
		temps_cpu = (temps_final - temps_initial) * 1e-6;
		//printf("temps apres lecture : %f\n",temps_cpu);
	temps_total=temps_total+temps_cpu;
		temps_initial = clock();

	delete[] fft_reel_shift_norm;
	delete[] fft_imag_shift_norm;
// désallocations qui cause un segfault!
  	delete[] plan_reel;
 	delete[] plan_imag;

} //fin de test de validite du nom de fichier
else{
	printf("fichier %i inexistant\n",cpt_angle);
	//fclose(test_existence);

	}
}//fin de boucle for sur tous les angles  on peut désallouer les variables définit hors de la boucle

delete[] masque;
delete[] cache_jumeau;

delete[] phased1;
delete[] phased2;
delete[] phased3;
delete[] phased4;

////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

printf("points_faux : %i,points_faux\n",points_faux);




//////////////////////////////////////////////////////exportation de centre

FILE* fichier_centre = NULL;
	//ouverture de ce fichier en écriture binaire
	fichier_centre = fopen( TEST_PATH "/centre.bin", "wb");
	//printf("reel_arc[100] : %ld\n",reel_arc[0]);
	for(int cpt=0;cpt<(4*Nxmax*Nymax);cpt++)
		{
		//Ecriture sur le disque
		fwrite(&centre[cpt],sizeof(int),1,fichier_centre);// /!\sizeofint vaut 32 bits
		}
	printf("ecriture de centre en 32 bits\nDimension centre : %d x%d",2*Nxmax,2*Nymax);
	fclose(fichier_centre);
// désallocations qui cause un segfault!
delete[] centre;

printf("jumeaux re-elimines : %i\n",jumeau_elimine);

/*int cpt_max=0;
//Recherche du MAX dans reel_arc
	for(int cpt=1;cpt<(64*Nxmax*Nxmax*Nymax);cpt++)
		{

			if(reel_arc[cpt]>reel_arc[cpt_max])
			{
				cpt_max=cpt;
			}
		}

printf("max de reel_arc %f \n",reel_arc[cpt_max]);*/

//////////////////////////////ecriture de sup_redon
printf("temps_total proj : %f \n",temps_total);
printf("coucou \n");
FILE* fichier_sup_redon = NULL;
//ouverture de ce fichier en écriture binaire
fichier_sup_redon = fopen( TEST_PATH "/sup_redon_C.bin", "wb");
//printf("reel_arc[100] : %ld\n",reel_arc[0]);
for(int cpt=0;cpt<(N_tab);cpt++)
	{
	//Ecriture sur le disque
	precision=sup_redon[cpt];
	fwrite(&precision,sizeof(precision),1,fichier_sup_redon);
	}
printf("ecriture de sup_redon.bin terminée,  \n");
fclose(fichier_sup_redon);

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////mesure du temps de calcul
/////////////////////////////////////////////////////////////////////////////////////

	printf("nb_proj: %i\n",nb_proj);
	printf("Nombre de centre exclus : %i\n",centres_exclus);
	temps_final = clock ();
	temps_cpu = (temps_final - temps_initial) * 1e-6/nb_proj;
	printf("temps moyen pour 1 angle: %f\n",temps_cpu);
	temps_cpu = (temps_final - temps_initial) * 1e-6;
	printf("temps total pour %i angle(s): %f\n", nb_proj, temps_cpu);
	temps_initial=clock();//enclenchement du chronometre




/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////normalisation par sup_redon
/////////////////////////////////////////////////////////////////////////////////////


	for(int cpt=0;cpt<N_tab;cpt++)
	{
		if (sup_redon[cpt]==0) ////////////////remplace les 0 de sup_redon par des 1---> evite la division par 0
		{
			//printf("sup_redon= %f \n" , sup_redon[cpt]);
			sup_redon[cpt]=1;
		}
	/*if (sup_redon[cpt]>30.0)
		{
			sup_redon[cpt]=30.0;
		}*/
		reel_arc[cpt] = reel_arc[cpt]/sup_redon[cpt];//normalisation par sup_redon
		imag_arc[cpt] = imag_arc[cpt]/sup_redon[cpt];
	}

	temps_final = clock ();
	temps_cpu = (temps_final - temps_initial) * 1e-6;
	printf("temps apres normalisation : %lf\n",temps_cpu);

	////////////////////////////on libere la memoire de sup_redon

	delete[] sup_redon;

//////////////interpolation dans Fourier


	/*FILE* fichier_reel_arc_norm_avant = NULL;
	//ouverture de ce fichier en écriture binaire
	fichier_reel_arc_norm_avant = fopen("/home/mat/tomo_test/reel_arc_norm_avant.bin", "wb");
	//printf("reel_arc[100] : %ld\n",reel_arc[0]);
	for(int cpt=0;cpt<(64*Nxmax*Nymax*Nxmax);cpt++)
		{
		//Ecriture sur le disque
		precision=reel_arc[cpt];
		fwrite(&reel_arc[cpt],8,1,fichier_reel_arc_norm_avant);
		}
printf("precision : %i",sizeof(precision));

	printf("ecriture de reel_arc_norm_avant.bin terminée \n");
	fclose(fichier_reel_arc_norm_avant);*/
printf("coucou_1\n");
/*
for (int x = 0; x < dv0rf; x++)//on balaye l'image, référentiel avec centre au milieu
			{
				for (int y = 0; y < dv0rf; y++)//on balaye l'image,
				{
					int z_min=0;
					int z_max=0;
					int cpt=(x+y*dv0rf+z_min*dv1xdv0rf);

					while(reel_arc[cpt]==0 && imag_arc[cpt]==0 && z_min<dv0rf-1)
					{
						z_min++;
						cpt=(x+y*dv0rf+z_min*dv1xdv0rf);
					}
					//printf("coucou1, cpt : %i, x : %i, y : %i, z_min : %i, z max : %i\n",cpt,x,y,z_min,z_max);
					while(z_min<dv0rf-2)
					{
						z_max=z_min+1;
						cpt=x+y*dv0rf+z_max*dv1xdv0rf;
						while(reel_arc[cpt]==0 && imag_arc[cpt]==0 && z_max<dv0rf-1)
						{
							z_max++;
							cpt=(x+y*dv0+z_max*dv1xdv0rf);
						}
						//printf("coucou2, cpt : %i, x : %i, y : %i, z_min : %i, z max:%i,dv0  : %i\n",cpt,x,y,z_min,z_max,dv0);
						if(z_max<dv0-1)
						{
							for(int z_interp=z_min+1;z_interp<z_max;z_interp++)//on y met la moyenne
							{

							cpt=(x+y*dv0rf+z_interp*dv1xdv0rf);
	//reel_arc[cpt]=((z_interp-z_min)*reel_arc[x+y*dv0rf+z_max*dv1xdv0rf]+(z_max-z_interp)*reel_arc[x+y*dv0rf+z_min*dv1xdv0rf])/(z_max-z_min);
	//imag_arc[cpt]=((z_interp-z_min)*imag_arc[x+y*dv0rf+z_max*dv1xdv0rf]+(z_max-z_interp)*imag_arc[x+y*dv0rf+z_min*dv1xdv0rf])/(z_max-z_min);
							}
						}
						z_min=z_max;


					}

				}

			}

*/


//interp3D(reel_arc,  dv0rf,dv0rf, dv0rf);
//interp3D(imag_arc,  dv0rf,dv0rf, dv0rf);

printf("coucou_2");
//////////////////////////////ecriture de papillon binarisé
double *papillon_masque=new double[N_tab];
for(int compteur=0;compteur<N_tab;compteur++)
	{
	if(reel_arc[compteur]==0)
	papillon_masque[compteur]=0;
	else
	papillon_masque[compteur]=1;
	}


FILE* fichier_papillon_masque = NULL;
		//ouverture de ce fichier en écriture binaire
		fichier_papillon_masque = fopen( TEST_PATH "/papillon_masque.bin", "wb");
		//printf("reel_arc[100] : %ld\n",reel_arc[0]);

		for(int cptb=0;cptb<N_tab;cptb++)
			{
			//Ecriture sur le disque
			precision=papillon_masque[cptb];
			fwrite(&precision,sizeof(precision),1,fichier_papillon_masque);
			//cout<<"\npapillon_masque(cpt)="<<papillon_masque[cptb];
			}

fclose(fichier_papillon_masque);
delete[] papillon_masque;
cout<<"\npapillon_masque(cpt) termine\n";

//////////////////////////////ecriture de reel_arc_norm
/*	FILE* fichier_reel_arc_norm = NULL;
	//ouverture de ce fichier en écriture binaire
	fichier_reel_arc_norm = fopen("/home/mat/reel_arc_norm.bin", "wb");
	//printf("reel_arc[100] : %ld\n",reel_arc[0]);
	for(int cpt=0;cpt<(N_tab);cpt++)
		{
		//Ecriture sur le disque
		precision=reel_arc[cpt];
		fwrite(&precision,sizeof(precision),1,fichier_reel_arc_norm);
		}
	printf("ecriture de reel_arc_norm.bin terminée \n");
fclose(fichier_reel_arc_norm);
*/
//////////////////////////////ecriture de imag_arc_norm

/*	FILE* fichier_imag_arc_norm = NULL;
	//ouverture de ce fichier en écriture binaire
	fichier_imag_arc_norm = fopen("/home/mat/imag_arc_norm_C.bin", "wb");
	//printf("reel_arc[100] : %ld\n",reel_arc[0]);
	for(int cpt=0;cpt<(N_tab);cpt++)
		{
		//Ecriture sur le disque
		precision=imag_arc[cpt];
		fwrite(&precision,sizeof(precision),1,fichier_imag_arc_norm);
		}
	printf("ecriture de imag_arc_norm.bin terminée \n");
	fclose(fichier_imag_arc_norm);
*/


/////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////circshift avant TF3D
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
printf("*******************************************\n");
printf("circshift avant TF3D \n");
printf("*******************************************\n");


/////////////////////creation des tableaux qui recoivent le circshift
double *reel_arc_shift=new double[N_tab];//partie reelle du papillon dans l'espace reciproque
double *imag_arc_shift=new double[N_tab];//partie imaginaire du papillon dans l'espace reciproque

temps_initial = clock();//enclenchement chronometre

circshift3D(reel_arc, reel_arc_shift, dv0rf, dv0rf, dv0rf);//ancienne valeur, causant des problèmes d'arrondis
circshift3D(imag_arc, imag_arc_shift, dv0rf, dv0rf, dv0rf);

//circshift3D(reel_arc, reel_arc_shift, 4*Nxmax/Rf, 4*Nymax/Rf, 4*Nxmax/Rf);ancienne valeur, causant des problèmes d'arrondis
//circshift3D(imag_arc, imag_arc_shift, 4*Nxmax/Rf, 4*Nymax/Rf, 4*Nxmax/Rf);


temps_final = clock ();
temps_cpu = (temps_final - temps_initial) * 1e-6;
printf("temps apres circshift 3D: %f\n",temps_cpu);
printf("*******************************************\n");
printf("TF3D...\n");
printf("*******************************************\n");
printf("2s calcul de la TF3D\n");
printf("     .\"\".    .\"\".\n");
printf("     |  |   /  /\n");
printf("     |  |  /  /\n");
printf("     |  | /  /\n");
printf("     |  |/  ;-.__ \n");
printf("     |  ` _/  /  /\n");
printf("     |  /` ) /  /\n");
printf("     | /  /_/\\_/ \n");
printf("     |/  /      |\n");
printf("     (  ' \\ '-  |\n");
printf("      \\    `.  /\n");
printf("       |      |\n");
printf("       |      |\n");

temps_initial = clock();//enclenchement chronometre
//////////////////////////////suppression des tableaux non circshiftés
delete[] reel_arc;
delete[] imag_arc;


//////////////////////////// TF3D
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//TF3D(reel_arc_shift,imag_arc_shift,res_reel, res_imag);

//fftw_plan_with_nthreads(2);
	int N3D=N_tab;
	//Déclaration des variables pour la FFT : entre,sortie et "fftplan"
	fftw_complex *in3D, *out3D;
	fftw_plan p3D;
	//Réservation memoire
	in3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3D);
	//out3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3D);
	// p = fftw_plan_dft_2d(nx,ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	//Récupération de l'image dans la partie reelle de l'entree
	for(int cpt=0;cpt<(N3D);cpt++)
	{
		in3D[cpt][0]=reel_arc_shift[cpt];
		in3D[cpt][1]=imag_arc_shift[cpt];
	}
	//libertation de la memoire d'entree
	delete[] reel_arc_shift;
	delete[] imag_arc_shift;
	//Réservation memoire
	//in3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3D);
	out3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3D);

	//calcul du plan, parametre servant a calculer et optimiser le FFT

	p3D=fftw_plan_dft_3d(dv0rf, dv0rf, dv0rf, in3D, out3D,FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p3D); // repeat as needed
	fftw_destroy_plan(p3D);
	fftw_free(in3D);



temps_final = clock ();
temps_cpu = (temps_final - temps_initial) * 1e-6;
printf("temps apres TF 3D: %f\n",temps_cpu);
printf("*******************************************\n");
printf("circshift et calcul du module\n");
printf("*******************************************\n");
temps_initial = clock();
double *final_reel=new double[N_tab];//partie reelle du resultat final
double *final_imag=new double[N_tab];//partie imaginaire du resultat final

	for(int cpt=0;cpt<(N3D);cpt++)
		{
		final_reel[cpt]=out3D[cpt][0];
		final_imag[cpt]=out3D[cpt][1];
		}
	 fftw_free(out3D);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////circshift apres TF3D
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////creation des tableaux qui recoivent le circshift
	double *final_reel_shift=new double[N_tab];//partie reelle du resultat final shifté (centré)
	double *final_imag_shift=new double[N_tab];//partie imaginaire du resultat final shifté (centré)

	circshift3D(final_reel, final_reel_shift, dv0rf, dv0rf, dv0rf);
	circshift3D(final_imag, final_imag_shift, dv0rf, dv0rf, dv0rf);
//////////////////////////////suppression des tableaux non circshiftés
	delete[] final_reel;
	delete[] final_imag;




/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////calcul du module carre
/////////////////////////////////////////////////////////////////////////////////////
	double *valeur_module_shift=new double[N3D];

	for(int cpt=0;cpt<(N3D);cpt++)
		{
		valeur_module_shift[cpt] = final_reel_shift[cpt]*final_reel_shift[cpt]+final_imag_shift[cpt]*final_imag_shift[cpt];//calcul du module
		}

	temps_final = clock ();
	temps_cpu = (temps_final - temps_initial) * 1e-6;
	printf("temps apres circshift 3D final et calcul du module: %f\n",temps_cpu);
	temps_initial = clock();//enclenchement chronometre


printf("coucou-1");
/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////ecriture des resultats
/////////////////////////////////////////////////////////////////////////////////////
	printf("*******************************************\n");
	printf("ecriture des résultats\n");
	printf("*******************************************\n");
 	printf("coucou0\n");
	FILE* fichier_final_reel_shift = NULL;
 	FILE* fichier_final_imag_shift = NULL;
 	FILE* fichier_final_modul_shift = NULL;

 printf("coucou1\n");

 	fichier_final_reel_shift = fopen( TEST_PATH "/final_reel_shift.bin", "wb");
 	fichier_final_imag_shift = fopen( TEST_PATH "/final_imag_shift.bin", "wb");
 	fichier_final_modul_shift = fopen( TEST_PATH "/final_modul_shift.bin", "wb");

	//ecrire_rapport(Nxmax,rayon,Rf,K,dimx_ccd2,coin_x,coin_y,precision,Chemin,nb_proj,n1,NA,Tp,G);
//
	for(int cpt=0;cpt<(N3D);cpt++)
		{
		precision=final_reel_shift[cpt];
		fwrite(&precision,sizeof(precision),1,fichier_final_reel_shift);

		precision=final_imag_shift[cpt];
		fwrite(&precision,sizeof(precision),1,fichier_final_imag_shift);

		precision=valeur_module_shift[cpt];
		fwrite(&precision,sizeof(precision),1,fichier_final_modul_shift);//ecriture du module
		}
//
printf("coucou2");
//Libérer la mémoire des 4 images expérimentales.
	delete[] phase1;
	delete[] phase2;
	delete[] phase3;
	delete[] phase4;
	delete[] reference;
//
printf("coucou3");
//on libere la TF 2D du front d'onde et son cropé
	delete[] fft_reel;
	delete[] fft_imag;
	delete[] fft_reel_tmp;
	delete[] fft_imag_tmp;
//libération memoire allouée pour les threads
 void fftw_cleanup_threads(void);

printf("ecriture de final_reel_shift.bin,  final_imag_shift.bin, final_modul_shift : OK \n");
fclose(fichier_final_reel_shift);
fclose(fichier_final_imag_shift);
fclose(fichier_final_modul_shift);

//liberation des variables "resultat final"
delete[] final_reel_shift;
delete[] final_imag_shift;
delete[] valeur_module_shift;

temps_arrivee = clock ();
temps_cpu = (temps_arrivee-temps_depart ) * 1e-6;
printf("temps total: %f\n",temps_cpu);
return 0;

}



/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//Interpolation3D : attend un volume "troué" et les taille en x,y,et z

void interp3D(double *volume_interp_3D, int taille_x,int taille_y,int taille_z)
{

int z=0;
// int cpt=x+y*dv0rf+z*taille_y*dv0rf;
int z_min=taille_z;
int z_max=0;

// cout<< "valeur papillon : "<< volume_interp_3D[cpt] << endl;
for (int x=0; x < taille_x; x++)//on balaye l'image, référentiel avec centre au milieu
	{
	for (int y=0; y<taille_y; y++)//on balaye l'image,
		{//printf("y:%i  x: : %i ",y,x);
				//cpt=x+y*taille_x+z*taille_y*taille_x;
// 				for (int z = 0; z < taille_z; z++)//on balaye l'image,
// 				{cpt=x+y*taille_x+z*taille_y*taille_x;
// // 					if(sup_redon[cpt]==599)
// // 					{
// // 					cout << "x : " << x <<"y : "<< y<< "z : "<<z<<"cpt : "<<cpt<<"sup_redon : " << sup_redon[cpt] << endl;
// // 					}
//
// 				}
		z=0;
		//printf("cpt : %i, sup redon : %d\n",x+y*taille_x+z*taille_y*taille_x,volume_interp_3D[x+y*taille_x+z*taille_y*taille_x]);
		while(z<taille_z) // pas de boucle for car l'indice est variable
			{float valeur_test=volume_interp_3D[x+y*taille_x+z*taille_y*taille_x];
			if(valeur_test!=0 && z<taille_z) //si valeur !=0, alors borne inférieure mais z_min ne peut valoir 19.
				{	//printf("coucou1 sup_redon : %f,x,y,z : %i,%i,%i\n",volume_interp_3D[x+y*taille_x+z*taille_y*taille_x],x,y,z);
				z_min=z; //on a trouvé z_min
				//printf("x,y: %i,%i z : %i, sup_redon : %f\n",x,y,z, volume_interp_3D[x+y*taille_x+z*taille_y*taille_x]);
				z_max=z_min+1; //initialiser z_max au point suivant (sinon while jamais verifie)

					while( z_max<taille_z && volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]==0) //soit trouver prochain z_max, soit fin de colonne
						{	//printf("coucou2\n,%i,%i,%i,sup_redon : %f\n",x,y,z_max,volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]);
						z_max=z_max+1;
						z=z_max;

						if(z_max==taille_z-1)
							{
							//printf("fin de colonne\n");
							}
						}


					if(z_max!=taille_z && z_max-z_min>1 && volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]!=0) //il faut au moins un trou pour interpoler
						{
						//printf("interpolation x,y,zmin, z_max: %i,%i,%i, %i\n", x,y,z_min,z_max);
						//y=ax+b-> interpolation linéaire
						double a=(volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]-volume_interp_3D[x+y*taille_x+z_min*taille_y*taille_x])/(z_max-z_min);
						double b=volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]-a*z_max;
						//printf("hello: z_min: %i, zmax:%i\n",z_min,z_max);
						for(int cpt_z=z_min+1; cpt_z<z_max; cpt_z++)

						{
						volume_interp_3D[x+y*taille_x+cpt_z*taille_y*taille_x]=a*cpt_z+b;
						//printf("interpolaion finie\n");
						}
// 						x;
// 						volume_interp_3D[x+y*taille_x+z_min*taille_y*taille_x];
// 						z_min;
// 						volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x];
// 						z_max;

						z=z_max-1; // nouveau compteur=ancienne borne sup (-1 car z++)

						}




								// redémarrer le compteur de ligne à z_max
				}
			z++;
			}

		}
	}
}
//##############§FIN INTERP3D#####################################################


void ecrire_rapport(int Nxmax,float rayon,float Rf, float K, int dimx_ccd2,int coin_x, int coin_y,double precision,char* chemin,int nb_proj,float n1,float NA,float Tp, int G)
{
//////////////////////////////////////////////Ecriture du rapport de calcul
time_t date;
time(&date);

FILE *fichier_rapport ;    /*  pointeur vers un type fichier     */
/*  ouverture du fichier 'mon_fichier.txt' pour ecriture (w) en mode texte (t)  */
fichier_rapport = fopen ( TEST_PATH "/rapport_calcul.txt", "wt") ;
/*  en cas d'échec de l'ouverture, le pointeur est NULL: intercepter ce cas  */
if (fichier_rapport == NULL){
    /*  message d'alerte et fin du programme  */
    printf ("impossible de créer le fichier rapport_calcul.txt\n") ;
    exit (0) ;}
fprintf(fichier_rapport,"Date     : %s\nNxmax    : %i\n,Rayon : %f\nRf       : %f\nK : %f\ndimx_ccd : %d\ncoin_x   : %d\ncoin_y   : %d\nPrecision: %lu bits\nSession  : %s\nnb_proj  : %i\nindice n1: %f\nNA       : %f\nT_pixel  : %e\nG        : %i",ctime(&date),Nxmax,rayon,Rf,K,dimx_ccd2,coin_x,coin_y,8*sizeof(precision),chemin,nb_proj,n1,NA,Tp,G);
/*  fermeture du fichier  */
fclose (fichier_rapport) ;
}

/////////////////////////////////////////////////////////////////////////////////////
//fonction chargeant une image 2D. Retourne un pointeur sur le tableau 1D qui la contient
//  Nécessite le nom du pointeur à remplir,
//le numéro du pas de phase shifting (1,2 3 ou 4), le chemin, et cpt_fichier
//pour boucler sur les 1000 projections
/////////////////////////////////////////////////////////////////////////////////////

unsigned char* charger_image2D(unsigned char* phasei,int numero,char *Chemin,int cpt_fichier,int coin_x,int coin_y,int taille_x,int taille_y)
{
	char CopieChemin[100];
	char FinNom[15];

	/*FILE* fichier_TF_real2 = NULL;*///nom bidon pour tester charger_image2D

	//sauvegarder chemin car sprintf écrase
	strcpy(CopieChemin,Chemin);
	sprintf(FinNom,"%ld-00%ld.bmp", (long int) cpt_fichier, (long int) numero);
	//créer le chemin final
	strcat(CopieChemin, FinNom);
	//printf("Copiechemin vaut : %s\n", CopieChemin);
	//charger l'image dans le tableau
///tester la validité du chemin "copiechemin"
/*FILE *test_chemin = fopen(CopieChemin,"r");

if (test_chemin == NULL) {
  fprintf( stderr, "attention, Ficher %i absent.\n",cpt_fichier); //fichier invalide? alors on retourne une erreur
  phasei=NULL;
return phasei;
fclose(test_chemin);
}
else //sinon on peut essayer de remplir le tableau
	{
		phasei = rempli_tableau(CopieChemin, coin_x, coin_y,taille_x,taille_y);
		return phasei;
fclose(test_chemin);
	}*/
phasei = rempli_tableau(CopieChemin, coin_x, coin_y,taille_x,taille_y);
		return phasei;
}

/////////////////////////////////////////////////////////////////////////////////////
//Fonction transferant les pixels de l'image 2D de taille (nx,ny) dans un tableau 1D de taille (nx*ny)
// nécessite la bibliotheque imageMagick
/////////////////////////////////////////////////////////////////////////////////////
unsigned char *rempli_tableau(string path, int coinx, int coiny,int taille_x, int taille_y)
{
	/////// declaration des variables

	Image Monimage;
	//int i, currentImageWidth, currentImageHeight;
	unsigned char *finalArray;

	////// chargement en memoire de l'image

	Monimage.read(path);

	////// Mon image est N&B
	Monimage.type( GrayscaleType );
	//Monimage.display();
	//Monimage.channelDepth();
	Monimage.compressType(NoCompression);
 	//Monimage.crop(Geometry(dimx_ccd,dimy_ccd, 0, 0) );
	//Monimage.write("/home/mat/tomo_test/test.bmp");
	////// extraction des caractéristiques de l'image
	//currentImageWidth = Monimage.columns();
	//currentImageHeight = Monimage.rows();
	////// reservation de la taille memoire du tableau final
	finalArray = new unsigned char[taille_x*taille_y];
	////// lecture de l'image
	Monimage.getPixels(coinx,coiny,taille_x,taille_y);
	////// ecriture de l'image dans un tableau
	Monimage.writePixels(GrayQuantum, finalArray);
	//affichage pour contrôle
	//printf("valeur de la case 256 %d \n", finalArray[256]);
	return finalArray;
	// ** // delete[] finalArray;
}





/////////////////////////////////////////////////////////////////////////////////////
//fonction de calcul TF2D
/////////////////////////////////////////////////////////////////////////////////////
void TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y)
{
	int N=taille_x*taille_y;
	//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
	fftw_complex *in, *out;
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

		fft_reel[cpt]=out[cpt][0];
		fft_imag[cpt]=out[cpt][1];
		}


	fftw_destroy_plan(p);
	fftw_free(in); fftw_free(out);
}

/////////////////////////////////////////////////////////////////////////////////////
//fonction de circshift "generale"
/////////////////////////////////////////////////////////////////////////////////////

double* circshift(double* entree, int dimx,int dimy,int decal_x,int decal_y)
{
	double * entree_shift=new double[dimx*dimy];
	for(int xi=0;xi<decal_x;xi++)
	{
		for(int yi=0;yi<decal_y;yi++)
		{
		int pixel=yi*dimx+xi;
		int pixel_shift=(yi+decal_y)*dimx+xi+decal_x;
	//1er quadrant vers 4 eme
		entree_shift[pixel_shift]=entree[pixel];
	//4 eme quadrant vers 1er
		entree_shift[pixel]=entree[pixel_shift];
	//2eme vers 3eme
	entree_shift[(yi+decal_y)*dimx+xi]=entree[pixel+decal_x];
	//3eme vers 2eme
	entree_shift[pixel+decal_x]=entree[(yi+decal_y)*dimx+xi];
		}
	}

return entree_shift;
delete[] entree_shift; //est ce utile dans une fonction?
}

/*void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D)
{

  for(int xi=0;xi<dimFinal3D.x;xi++)
  {
    for(int yi=0;yi<dimFinal3D.y;yi++)
    {
      for(int zi=0;zi<dimFinal3D.z;zi++)
      {
	///variables décalés
        short int x2=xi+decal3D.x;
        short int y2=yi+decal3D.y;
        short int z2=zi+decal3D.z;
	///compteur zi*taille_x*taille_y+yi*taille_x+xi
        if(xi+decal3D.x>dimFinal3D.x)
          x2=(xi+decal3D.x)-dimFinal3D.x;
        else
        {
          if(xi+decal3D.x<0)
            x2=dimFinal3D.x+(xi+decal3D.x);
          //cout<<"coucou"<<endl;
        }
        if(yi+decal3D.y>dimFinal3D.y)
          y2=(yi+decal3D.y)-dimFinal3D.y;
        else
        {
          if(yi+decal3D.y<0)
            y2=dimFinal3D.y+(yi+decal3D.y);
        }
        if(zi+decal3D.z>dimFinal3D.z)
        {
          z2=(zi+decal3D.z)-dimFinal3D.z;
        }
        else
        {
          if(zi+decal3D.z<0)
            z2=dimFinal3D.z+(zi+decal3D.z);
        }
        int pixel=zi*dimFinal3D.x*dimFinal3D.y+yi*dimFinal3D.x+xi;
        int pixel_decal=z2*dimFinal3D.x*dimFinal3D.y+y2*dimFinal3D.x+x2;
//         if(volume3D[pixel]==344444)
//         {
//           printf("xi+decal3D.x : %i,yi+decal3D.y : %i,zi+decal3D.z : %i\n",xi+decal3D.x,yi+decal3D.y,zi+decal3D.z);
//           printf("xi %i: ,yi : %i,zi %i: \n",xi,yi,zi);
//           printf("x2 %i: ,y2 : %i,z2 %i: \n",x2,y2,z2);
//           printf("volume3D_shift[pixel]: %f",volume3D_shift[pixel]);
//         }
        volume3D_shift[pixel_decal]=volume3D[pixel];
      }
    }
  }



}*/

//ancienne fonction circshift3D uniquement pour decalage d'une demi dimension
void circshift3D(double *volume3D, double *volume3D_shift,int taille_x,int taille_y,int taille_z)
{
//////attention il y a 8 cubes:
////// face superieur:  3 4     010 110
//                      1 2     000 100

//face inferieur :     7  8     011 111
//		       5  6     001 101

// 1   <-----> 8
// 2   <-----> 7
// 3   <-----> 6
// 4   <-----> 5

	for(int xi=0;xi<taille_x/2;xi++)
	{
			for(int yi=0;yi<taille_y/2;yi++)
			{
				for(int zi=0;zi<taille_z/2;zi++)
				{

			//cube 1 vers cube 8  //  000--->111
	//test volume3D_shift[zi*taille_x*taille_y+yi*taille_x+xi]=1;
		volume3D_shift[(zi+taille_x/2)*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi+taille_x/2]=volume3D[zi*taille_x*taille_y+yi*taille_x+xi];


			//cube 8 vers cube 1
		volume3D_shift[zi*taille_x*taille_y+yi*taille_x+xi]=volume3D[(zi+taille_x/2)*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi+taille_x/2];


			//cube 2 vers cube 7  //  100 ---->  011
		volume3D_shift[(zi+taille_x/2)*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi]=volume3D[zi*taille_x*taille_y+yi*taille_x+xi+taille_x/2];

			//cube 7 vers cube 2
		volume3D_shift[zi*taille_x*taille_y+yi*taille_x+xi+taille_x/2]=volume3D[(zi+taille_x/2)*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi];


			//cube 3 vers cube 6 // 010 ---> 101
		volume3D_shift[(zi+taille_x/2)*taille_x*taille_y+yi*taille_x+xi+taille_x/2]=volume3D[zi*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi];

			//cube 6 vers cube 3
		volume3D_shift[zi*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi]=volume3D[(zi+taille_x/2)*taille_x*taille_y+yi*taille_x+xi+taille_x/2];


			//cube 4 vers cube 5 // 110 --->  001
		volume3D_shift[(zi+taille_x/2)*taille_x*taille_y+yi*taille_x+xi]=volume3D[zi*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi+taille_x/2];

			//cube 5 vers cube 4
		volume3D_shift[zi*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi+taille_x/2]=volume3D[(zi+taille_x/2)*taille_x*taille_y+yi*taille_x+xi];



				}


			}
	}

}
