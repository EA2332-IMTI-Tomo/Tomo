// résolution des images d'entrée
#define IMAGE_X 1312
#define IMAGE_Y 1082
// taille de la fenêtre utile
#define WINDOW_X 1024
#define WINDOW_Y 1024
#define LAST_ANGLE_NUM 10

//taille de la fenêtre de sélection du cercle
#define CIRCLE_DIAMETER 242
#define FICHIER_MATTHIEU1 "/home/mat/tomo_test/TF2d/Tf2d.bin"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <octave-2.9.9/octave/oct.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>

#include <Magick++.h>
#include <fftw3.h>

#include <cstring>
#include <fstream>
#include <sstream>
#include <assert.h>

using namespace std;
using namespace Magick;


/* --------------------------------------------------------------------------- */
// Types
/* --------------------------------------------------------------------------- */

typedef struct {
        double Re,Im;
} nb_complexe;

typedef struct {
        int x,y;
} Var2D;

typedef struct {
        int x,y,z;
} Var3D;

/* --------------------------------------------------------------------------- */
// Prototypes
/* --------------------------------------------------------------------------- */

void changeDim2D(double* tab, double* tabFinal, Var2D dimInit,Var2D dimFin);
double *circshift(double *entree,int dimx,int dimy,int decal_x,int decal_y);

void circshift2(double* entree, double* result, Var2D dim,Var2D decal);
void circshift3(double* entree, double* result, Var2D dim,Var2D decal);
void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D);

void rempli_tableau(unsigned char *finalArray, string path, Var2D coin, Var2D taille);
void TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y);
void TF2D_INV(double entree_reelle[],double entree_imag[], double sortie_reelle[],double sortie_imag[],int taille_x,int taille_y);
double *tukey2D(int dimx,int dimy, float alpha);
void charger_image2D(unsigned char* phasei, char* chemin,int coin_x,int coin_y,int taille_x,int taille_y);

void circshift3D(double *volume3D, double *volume3D_shift,int taille_x,int taille_y,int taille_z);
void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D);
void genereCache(double masque[], int t_image, int t_mask, int centreX, int centreY);
void interp3D(double *volume_interp_3D, int taille_x,int taille_y,int taille_z);
void ecrire_rapport(int NXMAX,float rayon,float Rf, float K, int DIMX_CCD2,int coin_x, int coin_y,short int precision_exportation,char *chemin,int nb_proj,float n1,float NA,float Tp, int G);
//void fftw_plan_with_nthreads(int nthreads);
void multiplier_masque(double image[], unsigned char masque[], int t_image, int t_mask, int centreX, int centreY);
void concatener(char chaine1[],char chaine2[],char resultat[]);
void antigaussienne(double *tab, int Tx, int sigma, float A, int Exy);
void multiplier_masque2(double image[], double masque[], int t_image, int t_mask, int centreX, int centreY);

void crop_window(double* img_src, size_t src_dimx, size_t src_dimy, double* img_dst, size_t dst_dimx, size_t dst_dimy, size_t ul_corner_x, size_t ul_corner_y);


/* -------------------------------------------------------------------------- */
// Usage
/* -------------------------------------------------------------------------- */


static void
usage(int argc, char **argv)
{
        if ((argc - 1) == 0) {
                printf("Programme de reconstruction en tomographie hors-axe\n");
                printf("Usage: %s <paramètres obligatoires> <paramètres optionnels>\n", argv[0]);
                printf("Paramètres obligatoires: \n");
                printf("-i <répertoire>: répertoire des images acquises à traiter \n");
                printf("-c <cx> <cy>: centre du cercle sur l'image 1024x1024 (orig°: top-left corner) \n");

                exit(EXIT_FAILURE);
        }

}


/* --------------------------------------------------------------------------- */
// Main
/* --------------------------------------------------------------------------- */


int main(int argc, char *argv[])
{

        // ----------------------------------------------------------------------
        // parsing arguments
        // ----------------------------------------------------------------------

        usage(argc, argv);

        char* input_dir = NULL;
        size_t circle_cx = 0, circle_cy = 0;

        while (argc > 0) {
                if (!strcmp(argv[0], "-i") && (argc > 1)) {
                        input_dir = argv[1];
                        argc--;
                        argv++;
                }
                if (!strcmp(argv[0], "-c") && (argc > 2)) {
                        circle_cx = atoi(argv[1]);
                        circle_cy = atoi(argv[2]);
                        argc-=2;
                        argv+=2;
                }

                argc--;
                argv++;
        }

        assert(input_dir);
        assert(circle_cx && circle_cy);

        cout << endl << "input dir: " << input_dir;
        cout << endl << "circle centered on: (" << circle_cx << " , " << circle_cy << ")";
        cout.flush();


        // ----------------------------------------------------------------------
        //


        //VECTRA1024 2
        Var2D dimCCD= {WINDOW_X, WINDOW_Y},decalCCD= {dimCCD.x/2,dimCCD.y/2};
        int DIMX_CCD2=dimCCD.x,  DIMY_CCD2=dimCCD.y;
        //valeur controlant l'exclusion de certains centres

        int xm0_limite=200; //centre maximum
        int ym0_limite=200;
        int rayon_inf=xm0_limite*xm0_limite+ym0_limite*ym0_limite;

        //valeur du coin pour la découpe
        //découpe centrée - la découpe sera réalisée à la lecture du fichier image
        Var2D coin { (IMAGE_X - WINDOW_X)/2 , (IMAGE_Y - WINDOW_Y)/2 };
        int coin_x=coin.x,coin_y=coin.y;

        //   if(coin_x+DIMX_CCD2>IMAGE_X || coin_y+DIMY_CCD2> IMAGE_Y)//
        //     printf("Attention, la zone decoupee sort de l'image\n");

        // Définition de la fenêtre de sélection centrée sur (cx, cy)
        size_t upperleft_x = circle_cx - (CIRCLE_DIAMETER / 2);
        size_t upperleft_y = circle_cy - (CIRCLE_DIAMETER / 2);
        assert (upperleft_x < WINDOW_X);
        assert (upperleft_y < WINDOW_Y);

        // création d'une coupe 2D centrée sur le cercle
        double *decoupe=new double[CIRCLE_DIAMETER * CIRCLE_DIAMETER];

        // ----------------------------------------------------------------------
        //


        float precision;//pour conversion de type lors de l'écriture finale
        const short int precision_exportation=sizeof(precision);
        cout<<"precision_export="<<precision_exportation<<endl;
        int jumeau_elimine=0;
        int points_faux=0;
        int fftwThreadInit;

        fftwThreadInit=fftw_init_threads();
        int nthreads=3;
        printf("fftwThreadInit: %i\n",fftwThreadInit);


        //##VECTRA
        char* Dossier_acquiz= input_dir; //dossier des données
        char NumSession[]="i0";//partie fixe du nom du fichier

        char NomHoloPart1[strlen(Dossier_acquiz)+strlen("i0xxx")]; //réserver la place pour un nom montant au max à l'hologramme 999
        concatener(Dossier_acquiz,NumSession,NomHoloPart1);//concaténer le dossier avec la date, partie fixe au début du nom

        char Fichier_holo[256];//réserver la place pour le nom final

        //chemin pour référence d'amplitude
        char CheminModRef[strlen(Dossier_acquiz)+strlen("mod_ref.pgm")+1];
        concatener(Dossier_acquiz,"/mod_ref.pgm",CheminModRef);

        int coin_x_tmp;
        int coin_y_tmp;
        const int premier_plan=1, Num_Angle_final= LAST_ANGLE_NUM,//
                               NbAngle=Num_Angle_final-premier_plan+1, SautAngle=1, N=DIMX_CCD2*DIMY_CCD2;//nombre de pixel (pour le tableau 2D)
        /*printf("coin x (default : %i) : ",coin_x);//coordonnees du coin haut gauche dans l'image originale (748x570) pour l'image recadrée --
          scanf("%i", &coin_x);
          //if(coin_x_tmp==NULL)
          //printf("coin_x_tmp : %i\n",coin_x_tmp);
          printf("coin y (default : %i): ",coin_y);
          scanf("%i", &coin_y);*/
        //K; //Nombre total de pixel en x
        /*int Nty=DIMX_CCD2;*/ //(support carré)

        const float
        n1=1.515,	//indice de l'huile
           NA=1.40,	//ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
            theta=asin(NA/n1),//theta_max
            lambda=632*pow(10,-9), //longueur d'onde
            G=100,	//grossissement telan+objectif
                             //Tp=11.2*pow(10,-6), //Taille des pixels
            Tp=8*pow(10,-6),
            Rf=1/1.5,//1.2;=1/facteur grossissement
            K=lambda*G/(2*NA*Tp*Rf),// Facteur d'echelle
            Tps=Tp,//*K; //Taille pixel apres mise à l'echelle important si RF<>1
            tailleTheoPixelHolo=Tp/(G/Rf)*pow(10,9);//Pour info


        float
        rayon=n1*Rf*Tps*DIMX_CCD2/(G*lambda), //Rayon
              delta_zmax=rayon*cos(theta);//float?

        printf("theta : %f\n",theta);
        printf("K : %f\n", K );
        printf("Rayon %f \n",rayon);
        cout<<"###############################################"<<endl;
        //cout<<"Taille théorique des pixels espace holo:"<<tailleTheoPixelHolo<<" nm"<<endl;
        //printf("DIMX_CCD2 %i \n",DIMX_CCD2);

        //int rayon=round(n1*Rf*Tp*DIMX_CCD2/(G*lambda));ici un commentaire archéologique
        //int rayon=100;

        //dimension du support dans Fourier

        const int
        NXMAX=round(n1*Rf*Tps*DIMX_CCD2/(G*lambda)*NA/n1),//NXMAX=round(rayon*NA/n1)=round(rayon*sin theta_max);
              //Variable pour (éventuellement) Forcer la dimension de l'espace à une certaine taille, attention, ceci changera la taille des pixels tomo
              dim_final=4*NXMAX;
        //
        double zoom=double(dim_final)/double(4*NXMAX);
        cout<<"zoom="<<zoom<<endl;

        //NXMAX_Rf=round(zoom*n1*Tps*DIMX_CCD2/(G*lambda)*NA/n1); //dimension espace final
        //cout<<"NXMAX_Rf="<<NXMAX_Rf<<endl;

        const int NYMAX=NXMAX,
                  N_tab=64*NXMAX*NXMAX*NXMAX;//Nombre de pixel dans l'esapce final (indice max du tableau 3D+1)

        int centres_exclus=0, nb_proj=0;

        //printf("testNXMAX_Rf %f \n",NXMAX_Rf);
        //printf("round(testNXMAX_Rf) %f \n",round(NXMAX_Rf));
        //printf("NXMAX=%i,4*NXMAX=%i \n",NXMAX,4*NXMAX);

// float
        // test_NXMAX_Rf=n1*Tps*DIMX_CCD2/(G*lambda)*NA/n1;

        float tailleTheoPixelTomo=tailleTheoPixelHolo*DIMX_CCD2/(4*NXMAX);

        cout<<"Taille théorique pixel tomo:"<<tailleTheoPixelTomo<<endl;

        if(zoom!=1) {
                cout<<"dim_final forcée à "<<dim_final<<" au lieu de "<<4*NXMAX<<"-->facteur de zoom="<<zoom<<endl;
                cout<<"nouvelle taille pixel tomo imposée="<<tailleTheoPixelTomo/zoom<<"nm"<<endl;
                cout<<"###########################################"<<endl;
        }

        //Indice crée pour balayer le tableau final 1D (informatique) contenant le volume final de données en 3D
        //int N_tab=round(64*(n1*Tps*DIMX_CCD2/(G*lambda)*NA/n1)*(n1*Tps*DIMX_CCD2/(G*lambda)*NA/n1)*(n1*Tps*DIMX_CCD2/(G*lambda)*NA/n1));//sans arrondis

        printf("N_tab (indice espace 3D)= %i \n",N_tab);
        //cout <<"Espace final : " << round(4*zoom*n1*Tps*DIMX_CCD2/(G*lambda)*NA/n1) << " pixels cube\n";

        //printf("Précision d'exportation: %i bits\n chemin :%s \n",8*sizeof(precision),Chemin);

        cerr << endl << "dv0:" << round(4 * NXMAX) << " NXMAX:" << NXMAX;

        ///---------------Calcul de quelques constantes, afin d'éviter de les calculer dans la boucle--------------

        const float dv0=round(4*NXMAX),
                        dv0rf=round(4*NXMAX),
                        dv0s2rf=round(2*NXMAX),
                        dv1s2rf=round(2*NXMAX),
                        dv1xdv0rf=round(16*NXMAX*NXMAX),
                        dv2s2rf=round(2*NXMAX);


        printf("dimx espace 3D, dv0: %f\n",dv0);
        fflush(stdout);
        cout.flush();

        /// //////////////////////////////////////////////////fin calcul constante//////////////////////////


        //reservation de 4 tableaux (image 2D)  pour le phase shifting
        unsigned char* holo1=new unsigned char[DIMX_CCD2*DIMY_CCD2];
        double *holo_Re=new double[DIMX_CCD2*DIMY_CCD2];
        double *holo_Im=new double[DIMX_CCD2*DIMY_CCD2];

        double* ampli_ref=new double[DIMX_CCD2*DIMY_CCD2]; //reservation reference et "noir" camera pour une "tentative" de correction de la ref

        // 	unsigned char* noir_camera=new unsigned char[N];
        double* masque=new double[DIMX_CCD2*DIMY_CCD2]; //masque hamming
        unsigned char* cache_jumeau=new unsigned char[100];//cache objet jumeau

        //variables utilisées pour la TF2D
        //Avant Crop à NXMAX;
        double *fft_reel_tmp=new double[DIMX_CCD2*DIMY_CCD2];
        double *fft_imag_tmp=new double[DIMX_CCD2*DIMY_CCD2];


        //Après Crop à NXMAX;
        double *fft_reel=new double[4*NXMAX*NYMAX];
        double *fft_imag=new double[4*NXMAX*NYMAX];
        double *fft_reel_norm=new double[4*NXMAX*NYMAX];
        double *fft_imag_norm=new double[4*NXMAX*NYMAX];
        double *fft_module=new double[4*NXMAX*NYMAX];

        double  *holo_Re_shift=new double[N];
        double  *holo_Im_shift=new double[N];

        double *fft_module_shift=new double[4*NXMAX*NYMAX];

        double *fft_reel_shift=new double[4*NXMAX*NYMAX];
        double *fft_imag_shift=new double[4*NXMAX*NYMAX];
        ///-TF normalisée en amplitude cplx
        double  *fft_reel_shift_norm=new double[4*NXMAX*NYMAX];
        double  *fft_imag_shift_norm=new double[4*NXMAX*NYMAX];

        //double *qualite_phase_shift=new double[N];
        //espace3D reel, imaginaire et support de redondance

        double *reel_arc=new double[N_tab];//partie reelle du papillon dans l'espace reciproque
        double *imag_arc=new double[N_tab];//partie imaginaire du papillon dans l'espace reciproque
        double *sup_redon=new double[N_tab];//pour normaliser par rapport au nombre de fois où l'on remplit la frequence
        int *centre=new int[4*NXMAX*NYMAX];//pour mettre la position des centres translatés, on crée une variable 2D de la taille d'un plan apres tomo

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
        ///Chargement de la réference en module et calcul amplitude
        unsigned char* mod_ref=new unsigned char[N]; //reservation reference et "noir" camera pour une "tentative" de correction de la ref
        rempli_tableau(mod_ref,CheminModRef,coin,dimCCD);
        for(int pixel=0; pixel<N; pixel++) {
                ampli_ref[pixel]=sqrt((double)mod_ref[pixel]);
        }

        ///Variable de masquage :  fenetre de Tukey, écrasement jumeau
        float alpha=0.1;//coeff pour le masque de tuckey
        masque=tukey2D(DIMX_CCD2,DIMY_CCD2,alpha);
        //taille totale du masque en pixel pour elimination objet jumeau;
        const int t_mask=31;
        //cache_jumeau=rempli_tableau("../masques/k_20x20.bmp", 0, 0,t_mask,t_mask);
        double* cache_jumeau2=new double[t_mask*t_mask];
        antigaussienne(cache_jumeau2,t_mask,round(t_mask*2),1,0);
        /*  FILE * pEcriture1;
            pEcriture1=fopen("/home/mat/tomo_test/gaussienne.bin","w+");
            if(pEcriture1==0)
            cout<<"problème lors de l'écriture"<<endl;
            fwrite (cache_jumeau2,1,sizeof(double)*t_mask*t_mask,pEcriture1);
            fclose (pEcriture1);*/

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



        Var2D dim2D= {2*NXMAX,2*NYMAX},decal2D= {NXMAX,NYMAX};
        // ----------------------------------------------------------------------
        //
        // ----------------------------------------------------------------------
        //début de la grosse boucle sur les angles/////////////////////////
        //##############################################################################################################################


        char TamponChemin[ strlen(input_dir) + 256 ];
        for(int cpt_angle=premier_plan; cpt_angle<premier_plan+NbAngle; cpt_angle=cpt_angle+SautAngle) { //boucle sur tous les angles
                //cout<<"cpt_angle : "<<cpt_angle<<endl;
                if((cpt_angle-100*(cpt_angle/100))==0)
                        printf("cpt_angle=%i\n",cpt_angle);


                /*
                char TamponChemin[sizeof(Fichier_holo)];
                char Num_angle[4];
                sprintf(Num_angle,"%i",cpt_angle);
                //créer le chemin final
                concatener(NomHoloPart1,Num_angle,Fichier_holo);
                strcpy(TamponChemin,Fichier_holo);
                concatener(TamponChemin,".pgm",TamponChemin);
                */
                ///Concaténer le numéro d'angle dans le nom
                sprintf(TamponChemin, "%s/i%04i_ha.pgm", input_dir, cpt_angle);
                test_existence = fopen(TamponChemin, "rb");
                cout<<"";

                if(test_existence!=NULL) {
                        fclose(test_existence);
                        fprintf(stdout, "%s\n", TamponChemin);
                        //continue;
                        ///charger l'hologramme hors axe
                        charger_image2D(holo1,TamponChemin, coin_x, coin_y,DIMX_CCD2,DIMY_CCD2);

                        for(int pixel=0; pixel<N; pixel++) {
                                ampli_ref[pixel]=sqrt((double)mod_ref[pixel]);
                                holo_Re[pixel]=(double)holo1[pixel];
                        }

                        //////////////////////////// Circshift avant TF2D
                        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        circshift2(holo_Re,holo_Re_shift, dimCCD,decalCCD);
                        memset(holo_Im_shift, 0, (dimCCD.x*dimCCD.y)*sizeof(double));	  //pas besoin de décaler la partie imag car pleine de zéro.
                        //////////////////////////TF 2D du front d'onde: passage du plan image au plan réciproque////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        FILE *fichier_holo2d;
                        fichier_holo2d= fopen("/home/mat/tomo_test/TF2d/holo2d", "a+b");
                        if(fichier_holo2d==0)
                                cout<<"Erreur d'ouverture du fichier TF2d"<<endl;
                        for(int cpt=0; cpt<(dimCCD.x*dimCCD.y); cpt++) {
                                //Ecriture sur le disque
                                fwrite(&holo_Re_shift[cpt],sizeof(double),1,fichier_holo2d);
                        }
                        fclose(fichier_holo2d);

                        TF2D(holo_Re_shift,holo_Im_shift,fft_reel_tmp,fft_imag_tmp,dimCCD.x,dimCCD.y);

                        //Décalage  sur fft_reel_tmp, pour recnetrer le spectre avant découpe (pas obligatoire mais plus clair)
                        circshift2(fft_reel_tmp, fft_reel_shift, dimCCD,decalCCD);
                        circshift2(fft_imag_tmp, fft_imag_shift, dimCCD,decalCCD);

                        cout<<"dimCCD.x: "<<dimCCD.x<<endl;
                        cout<<"decalCCD: "<<decalCCD.x<<endl;

                        FILE *fichier_TF2d;
                        fichier_TF2d= fopen(FICHIER_MATTHIEU1, "a+b");
                        if(fichier_TF2d==0)
                                cout<<"Erreur d'ouverture du fichier TF2d"<<endl;
                        for(int cpt=0; cpt<(dimCCD.x*dimCCD.y); cpt++) {
                                //Ecriture sur le disque
                                fwrite(&fft_reel_shift[cpt],sizeof(double),1,fichier_TF2d);
                        }
                        fclose(fichier_TF2d);




                        crop_window(fft_reel_shift, 1024, 1024, decoupe, 2*NXMAX, 2*NYMAX, upperleft_x, upperleft_y);
                        // copie de la couche là où il faudrait
                        for (size_t s = 0; s < 2*NXMAX * 2*NXMAX; s++)
                                fft_reel[s] = decoupe[s];

                        crop_window(fft_imag_shift, 1024, 1024, decoupe, 2*NXMAX, 2*NYMAX, upperleft_x, upperleft_y);
                        // copie de la couche là où il faudrait
                        for (size_t s = 0; s < 2*NXMAX*2*NXMAX; s++)
                                fft_imag[s] = decoupe[s];


                        /*
                              FILE *fichier_TF2d_crop;
                        	  fichier_TF2d_crop= fopen("/home/mat/tomo_test/TF2d/tf_crop", "a+b");
                        	  if(fichier_TF2d_crop==0)
                                    cout<<"Erreur d'ouverture du fichier tf_crop"<<endl;
                        	  for(int cpt=0;cpt<(2*NXMAX*2*NXMAX);cpt++)
                        	    {
                        	      //Ecriture sur le disque
                        	      fwrite(&fft_reel[cpt],sizeof(double),1,fichier_TF2d_crop);
                        	    }
                        	  fclose(fichier_TF2d_crop);*/


                        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////circshift 2D apres TF
                        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        //Sauver les Tf2D partie réelle
                        /*char Nom_TF[100]="/home/mat/tomo_test/TF2d/TF2D_";
                          char Num_TF[10];
                          sprintf(Num_TF,"%i.bin",cpt_angle);
                          strcat(Nom_TF,Num_TF);*/

                        /*
                        FILE* fichier_tf2d_shift;
                        fichier_tf2d_shift= fopen(FICHIER_MATTHIEU2, "a+b");
                        if(fichier_tf2d_shift==0)
                              cout<<"Erreur d'ouverture du fichier TF2D"<<endl;
                        for(int cpt=0;cpt<(4*NXMAX*NYMAX);cpt++)
                          {
                            //Ecriture sur le disque
                            fwrite(&fft_reel_shift[cpt],sizeof(double),1,fichier_tf2d_shift);
                          }
                        fclose(fichier_tf2d_shift);
                        */

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////On cherche la position du maximum d'intensité					///
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        //----------------------------Recherche du maximum dans ref non centré-----------------------------------------

                        int cpt_max=0;
                        fft_module[0]=pow(fft_reel[0],2)+pow(fft_imag[0],2);

                        //Recherche du MAX dansle module
                        for(int cpt=1; cpt<(4*NXMAX*NYMAX); cpt++) {
                                fft_module[cpt]=pow(fft_reel[cpt],2)+pow(fft_imag[cpt],2);
                                if(fft_module[cpt]>fft_module[cpt_max]) {
                                        cpt_max=cpt;
                                }
                        }

                        FILE *fichier_TF2d_crop;
                        fichier_TF2d_crop= fopen("/home/mat/tomo_test/TF2d/tf_crop", "a+b");
                        if(fichier_TF2d_crop==0)
                                cout<<"Erreur d'ouverture du fichier tf_crop"<<endl;
                        for(int cpt=0; cpt<(2*NXMAX*2*NXMAX); cpt++) {
                                //Ecriture sur le disque
                                fwrite(&fft_reel[cpt],sizeof(double),1,fichier_TF2d_crop);
                        }
                        fclose(fichier_TF2d_crop);

                        /////////////////////////////
                        //on met 1 dans l'image centre en cpt_max
                        /////////////////////////////

                        //centre[cpt_max]++;

                        /*FILE* fichier_TF_mod = NULL;
                          fichier_TF_mod = fopen("/home/mat/tomo_test/fft_module_shift.bin", "a+b");
                          for(int cpt=0;cpt<(4*NXMAX*NYMAX);cpt++)
                          {
                          //Ecriture sur le disque
                          fwrite(&fft_module_shift[cpt],sizeof(double),1,fichier_TF_mod);
                          }
                          fclose(fichier_TF_mod);
                        */


                        double max_part_reel = fft_reel[cpt_max];
                        double max_part_imag = fft_imag[cpt_max];
                        double max_module = fft_module[cpt_max];

                        //printf("max_part_reel : %f, max_part_imag: %f,  max_module %f \n ",max_part_reel, max_part_imag, max_module);


                        /*FILE * TF_reel_shift;
                          TF_reel_shift = fopen("/home/mat/tomo_test/test_resultat_shift.bin", "wb");
                          for(int cpt=0;cpt<(4*NXMAX*NYMAX);cpt++)
                          {
                          //Ecriture sur le disque
                          fwrite(&fft_reel_shift[cpt],sizeof(double),1,TF_reel_shift);
                          }
                          fclose(TF_reel_shift);*/

                        ///Coordonnées dans l'espace 2D à partir de l'indice 1D: (xc,yc)=(0,0)=en haut à gauche de l'image
                        int xmi=cpt_max%(2*NXMAX);
                        int ymi=cpt_max/(2*NYMAX);
                        cout<<"xmi="<<xmi<<"//ymi="<<ymi<<endl;
                        int xc=xmi-NXMAX;
                        int yc=ymi-NYMAX;




                        ///sauver tf  avant normalisation
                        FILE* fichier_fft_reel;
                        fichier_fft_reel = fopen("/home/mat/tomo_test/TF2d_apres_masquage/fft_reel.bin", "a+b");//a+b=mode incrémental=tout écrire dans le même fichier
                        for(int cpt=0; cpt<4*NXMAX*NYMAX; cpt++) {
                                precision=fft_reel[cpt];
                                fwrite(&precision,sizeof(precision),1,fichier_fft_reel);
                        }
                        fclose(fichier_fft_reel);
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////normalisation par le pic central
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                        ///Normalisation
                        for(int cpt=0; cpt<(4*NXMAX*NYMAX); cpt++) {
                                fft_reel_norm[cpt]=(fft_reel[cpt]*max_part_reel+fft_imag[cpt]*max_part_imag)/max_module;
                                fft_imag_norm[cpt]=(fft_imag[cpt]*max_part_reel-fft_reel[cpt]*max_part_imag)/max_module;
                        }


                        /*FILE* TF_imag_shift_norm=NULL;
                          TF_imag_shift_norm = fopen("/home/mat/tomo_test/test_resultat_imag_shift_norm_C.bin", "a+b");
                          for(int cpt=0;cpt<(4*NXMAX*NYMAX);cpt++)
                          {
                          //Ecriture sur le disque
                          fwrite(&fft_imag_shift[cpt],sizeof(double),1,TF_imag_shift_norm);
                          }
                          fclose(TF_imag_shift_norm);*/


                        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////Mapping 3D
                        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                        //coordonnée dans l'image2D centrée (xm0,ym0)=(0,0)=au centre de l'image
                        int xm0=(xmi-NXMAX);//
                        printf("xm0 %i \n", xm0);
                        int ym0=(ymi-NYMAX);
                             printf("ym0 %i \n", ym0);
                        if(xm0==0 && ym0==0)
                                printf("(xm0,ym0)=(0,0) dans le plan %i\n", cpt_angle);
                        //if((xm0*xm0+ym0*ym0)>rayon_inf || fabs(xm0)>xm0_limite || fabs(ym0)>ym0_limite)
                        if((xm0*xm0+ym0*ym0)>rayon_inf) { //inutile?
                                //printf("xm0 : %i, ym0:%i\n",xm0,ym0);
                                //printf("\nFichier %i sans signal\n", cpt_angle);
                                centres_exclus++;
                        } else {
                                centre[xmi*2*NXMAX+ymi]=cpt_angle;//sauvegarde des centres en coordonnées non centrée; on met le numero d'angle
                                //pour vérifier les pb. A ouvrir en sizeof(int)=32 bits mais non signé
                                //printf("cpt_angle:%i, xmi: %i, ymi : %i\n",cpt_angle,xmi,ymi);

                                //création de variable pour éviter N calculs dans la boucle sur le volume 3D
                                double r2=rayon*rayon;

                                //indice 3D
                                int k=0;
                                double arg_z_arc=0;
                                double z_arc=0;
                                double zm0;


                                //printf("round(rayon*rayon-(xm0)^2-(ym0)^2: %i\n",round(rayon*rayon-(xm0)^2-(ym0)^2));

                                double zm0_carre = rayon*rayon-xm0*xm0-ym0*ym0;
                                if(round(zm0_carre)>-1) {

                                        zm0=sqrt(zm0_carre);
                                        //printf("zm0:%i\n",zm0);
                                        nb_proj++;

                                        /* Essai en calculant x et y a partir de cpt
                                           for(int cpt=0;cpt<4*NXMAX*NYMAX;cpt++)
                                           {
                                           int x = (cpt%(2*NXMAX));
                                           int y = (cpt/(2*NYMAX));
                                           }*/
                                        temps_initial = clock ();
                                        int NXMAX_CARRE=NXMAX*NXMAX;

                                        for (int y = -NYMAX; y < NYMAX; y++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
                                                int y_carre=y*y;
                                                for (int x = -NXMAX; x < NXMAX; x++) { //on balaye l'image 2D en y, centre au milieu

                                                        int cpt=(y+NYMAX)*2*NXMAX+x+NXMAX;//calcul du cpt du tableau 1D de l'image 2D


                                                        if(x*x+y_carre<NXMAX_CARRE) { //ne pas depasser l'ouverture numérique pour 1 hologramme
                                                                double z_carre=r2-x*x-y_carre; //altitude au carré des données

                                                                double z=round(sqrt(z_carre)-zm0);
                                                                double altitude=(z+dv2s2rf)*dv1xdv0rf; //donne n'importequoi sans l'arrondi sur z!!

                                                                k=(-xm0+x+dv0s2rf)+(-ym0+y+dv1s2rf)*dv0rf+round(altitude);//indice du tableau 1D du volume 3D

                                                                //cout<<"k"<<k<<endl;

                                                                reel_arc[k]+=fft_reel_norm[cpt];//pour calculer l'image
                                                                imag_arc[k]+=fft_imag_norm[cpt];//pour calculer l'image
                                                                sup_redon[k]+=1;//pour calculer le support
                                                        } else
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


                } //fin de test de validite du nom de fichier
                else {
                        printf("fichier %i inexistant\n",cpt_angle);
                        //fclose(test_existence);

                }
        }//fin de boucle for sur tous les angles  on peut désallouer les variables définit hors de la boucle




        delete[] fft_reel_shift;
       delete[] fft_imag_shift;
        delete[] fft_module_shift;

         delete[] fft_reel_shift_norm;
        delete[] fft_imag_shift_norm;

        delete[] fft_reel_norm;
        delete[] fft_imag_norm;
        delete[] fft_module;
   // delete[] holo_Im_shift;
       //  delete[] holo_Re_shift;
            delete[] masque;
       delete[] cache_jumeau;
        delete[] holo1;
         delete[] holo_Re;
 //delete[] fft_reel_shift;
        ////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////

        printf("points_faux!! : %i,points_faux\n",points_faux);




        //////////////////////////////////////////////////////exportation de centre

        FILE* fichier_centre = NULL;
        //ouverture de ce fichier en écriture binaire
        fichier_centre = fopen("/home/mat/tomo_test/centre.bin", "wb");
        //printf("reel_arc[100] : %ld\n",reel_arc[0]);
        for(int cpt=0; cpt<(4*NXMAX*NYMAX); cpt++) {
                //Ecriture sur le disque
                fwrite(&centre[cpt],sizeof(int),1,fichier_centre);// /!\sizeofint vaut 32 bits
        }
        printf("ecriture de centre en 32 bits\nDimension centre : %d x%d\n",2*NXMAX,2*NYMAX);
        fclose(fichier_centre);
        // désallocations qui cause un segfault!
        delete[] centre;

        printf("jumeaux re-elimines : %i\n",jumeau_elimine);

        /*int cpt_max=0;
        //Recherche du MAX dans reel_arc
        for(int cpt=1;cpt<(64*NXMAX*NXMAX*NYMAX);cpt++)
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
        fichier_sup_redon = fopen("/home/mat/tomo_test/sup_redon_C.bin", "wb");
        //printf("reel_arc[100] : %ld\n",reel_arc[0]);
        for(int cpt=0; cpt<(N_tab); cpt++) {
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


        for(int cpt=0; cpt<N_tab; cpt++) {
                if (sup_redon[cpt]==0) { ////////////////remplace les 0 de sup_redon par des 1---> evite la division par 0
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




        /*FILE* fichier_reel_arc_norm_avant = NULL;
        //ouverture de ce fichier en écriture binaire
        fichier_reel_arc_norm_avant = fopen("/home/mat/tomo_test/reel_arc_norm_avant.bin", "wb");
        //printf("reel_arc[100] : %ld\n",reel_arc[0]);
        for(int cpt=0;cpt<(64*NXMAX*NYMAX*NXMAX);cpt++)
        {
        //Ecriture sur le disque
        precision=reel_arc[cpt];
        fwrite(&reel_arc[cpt],8,1,fichier_reel_arc_norm_avant);
        }
        printf("precision : %i",sizeof(precision));

        printf("ecriture de reel_arc_norm_avant.bin terminée \n");
        fclose(fichier_reel_arc_norm_avant);*/

        //////////////interpolation dans Fourier
        //interp3D(reel_arc,  dv0rf,dv0rf, dv0rf);
        //interp3D(imag_arc,  dv0rf,dv0rf, dv0rf);

        printf("coucou_2");
        //////////////////////////////ecriture de papillon binarisé
        double *papillon_masque=new double[N_tab];
        for(int compteur=0; compteur<N_tab; compteur++) {
                if(reel_arc[compteur]==0)
                        papillon_masque[compteur]=0;
                else
                        papillon_masque[compteur]=1;
        }


        FILE* fichier_papillon_masque = NULL;
        //ouverture de ce fichier en écriture binaire
        fichier_papillon_masque = fopen("/home/mat/papillon_masque.bin", "wb");
        //printf("reel_arc[100] : %ld\n",reel_arc[0]);

        for(int cptb=0; cptb<N_tab; cptb++) {
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
        Var3D
        decal3D= {dv0rf/2,dv0rf/2,dv0rf/2},
                 dimFinal= {dv0rf,dv0rf,dv0rf};
        circshift3D2(reel_arc, reel_arc_shift, dimFinal,decal3D);//ancienne valeur, causant des problèmes d'arrondis
        circshift3D2(imag_arc, imag_arc_shift, dimFinal,decal3D);

        //circshift3D(reel_arc, reel_arc_shift, dv0rf, dv0rf, dv0rf);//ancienne valeur, causant des problèmes d'arrondis
        //circshift3D(imag_arc, imag_arc_shift, dv0rf, dv0rf, dv0rf);

        //circshift3D(reel_arc, reel_arc_shift, 4*NXMAX/Rf, 4*NYMAX/Rf, 4*NXMAX/Rf);ancienne valeur, causant des problèmes d'arrondis
        //circshift3D(imag_arc, imag_arc_shift, 4*NXMAX/Rf, 4*NYMAX/Rf, 4*NXMAX/Rf);


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

        fftw_plan_with_nthreads(nthreads);
        int N3D=N_tab;
        //Déclaration des variables pour la FFT : entre,sortie et "fftplan"
        fftw_complex *in3D, *out3D;
        fftw_plan p3D;
        //Réservation memoire
        in3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3D);
        //out3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3D);
        // p = fftw_plan_dft_2d(nx,ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

        //Récupération de l'image dans la partie reelle de l'entree
        for(int cpt=0; cpt<(N3D); cpt++) {
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
        void fftw_cleanup_threads(void);



        temps_final = clock ();
        temps_cpu = (temps_final - temps_initial) * 1e-6;
        printf("temps apres TF 3D: %f\n",temps_cpu);
        printf("*******************************************\n");
        printf("circshift et calcul du module\n");
        printf("*******************************************\n");
        temps_initial = clock();
        double *final_reel=new double[N_tab];//partie reelle du resultat final
        double *final_imag=new double[N_tab];//partie imaginaire du resultat final

        for(int cpt=0; cpt<(N3D); cpt++) {
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

        //circshift3D(final_reel, final_reel_shift, dv0rf, dv0rf, dv0rf);
        //circshift3D(final_imag, final_imag_shift, dv0rf, dv0rf, dv0rf);
        circshift3D2(final_reel, final_reel_shift, dimFinal, decal3D);
        circshift3D2(final_imag, final_imag_shift, dimFinal, decal3D);
        //////////////////////////////suppression des tableaux non circshiftés
        delete[] final_reel;
        delete[] final_imag;
        /////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////calcul du module carre
        /////////////////////////////////////////////////////////////////////////////////////
        double *valeur_module_shift=new double[N3D];

        for(int cpt=0; cpt<(N3D); cpt++) {
                valeur_module_shift[cpt] = final_reel_shift[cpt]*final_reel_shift[cpt]+final_imag_shift[cpt]*final_imag_shift[cpt];//calcul du module
        }
        temps_final = clock ();
        temps_cpu = (temps_final - temps_initial) * 1e-6;
        printf("temps apres circshift 3D final et calcul du module: %f\n",temps_cpu);
        temps_initial = clock();//enclenchement chronometre

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

        fichier_final_reel_shift = fopen("/home/mat/tomo_test/final_reel_shift.bin", "wb");
        fichier_final_imag_shift = fopen("/home/mat/tomo_test/final_imag_shift.bin", "wb");
        fichier_final_modul_shift = fopen("/home/mat/tomo_test/final_modul_shift.bin", "wb");

        ecrire_rapport(NXMAX,rayon,Rf,K,DIMX_CCD2,coin_x,coin_y,precision_exportation,Fichier_holo,nb_proj,n1,NA,Tp,G);
        //
        for(int cpt=0; cpt<(N3D); cpt++) {
                precision=final_reel_shift[cpt];
                fwrite(&precision,sizeof(precision),1,fichier_final_reel_shift);

                precision=final_imag_shift[cpt];
                fwrite(&precision,sizeof(precision),1,fichier_final_imag_shift);

                precision=valeur_module_shift[cpt];
                fwrite(&precision,sizeof(precision),1,fichier_final_modul_shift);//ecriture du module
        }
        //
        printf("coucou2");


        delete[] ampli_ref;
        delete[] mod_ref;
        //
        printf("coucou3");
        //on libere la TF 2D du front d'onde et son cropé
        delete[] fft_reel;
        delete[] fft_imag;
        delete[] fft_reel_tmp;
        delete[] fft_imag_tmp;
        //libération memoire allouée pour les threads
        void fftw_cleanup_threads(void);


        fclose(fichier_final_reel_shift);
        fclose(fichier_final_imag_shift);
        fclose(fichier_final_modul_shift);
        printf("ecriture de final_reel_shift.bin,  final_imag_shift.bin, final_modul_shift : OK \n");
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

///####################fonction################"

void changeDim2D(double* tab, double* tabFinal, Var2D dimInit, Var2D dimFin)
{
        int diffX=round(dimFin.x-dimInit.x)/2;
        int diffY=round(dimFin.y-dimInit.y)/2;

        for(int x=0; x<dimInit.x; x++) {
                for( int y=0; y<dimInit.y; y++) {
                        int cpt_init=y*dimInit.x+x;
                        int cpt_final=(y+diffY)*dimFin.x+x+diffX;
                        tabFinal[cpt_final]=tab[cpt_init];
                }
        }
}



void concatener(char chaine1[], char chaine2[], char resultat[])//fonction limitée mais intéret car copie locale (pas besoin de tampon)
{

        //récupérer la longueur du resultat
        size_t LongrChaine1=strlen(chaine1);
        size_t LongrChaine2=strlen(chaine2);
        short unsigned int LongrChaine=LongrChaine1+LongrChaine2;
        char tampon[LongrChaine];//+1 pour zero final?

        strcpy(tampon,chaine1);
        strcat(tampon,chaine2);

        for(short unsigned int cpt=0; cpt<LongrChaine+1; cpt++)
                resultat[cpt]=tampon[cpt];

}

//Interpolation3D : attend un volume "troué" et les dimensions en x,y,et z

void interp3D(double *volume_interp_3D, int taille_x,int taille_y,int taille_z)
{
        int z=0;
        // int cpt=x+y*dv0rf+z*taille_y*dv0rf;
        int z_min=taille_z;
        int z_max=0;

        // cout<< "valeur papillon : "<< volume_interp_3D[cpt] << endl;
        for (int x=0; x < taille_x; x++) { //on balaye l'image, référentiel avec centre au milieu
                for (int y=0; y<taille_y; y++) { //on balaye l'image,
                        //printf("y:%i  x: : %i ",y,x);
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
                        while(z<taille_z) { // pas de boucle for car l'indice est variable
                                float valeur_test=volume_interp_3D[x+y*taille_x+z*taille_y*taille_x];
                                if(valeur_test!=0 && z<taille_z) { //si valeur !=0, alors borne inférieure mais z_min ne peut valoir 19.
                                        //printf("coucou1 sup_redon : %f,x,y,z : %i,%i,%i\n",volume_interp_3D[x+y*taille_x+z*taille_y*taille_x],x,y,z);
                                        z_min=z; //on a trouvé z_min
                                        //printf("x,y: %i,%i z : %i, sup_redon : %f\n",x,y,z, volume_interp_3D[x+y*taille_x+z*taille_y*taille_x]);
                                        z_max=z_min+1; //initialiser z_max au point suivant (sinon while jamais verifie)

                                        while( z_max<taille_z && volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]==0) { //soit trouver prochain z_max, soit fin de colonne
                                                //printf("coucou2\n,%i,%i,%i,sup_redon : %f\n",x,y,z_max,volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]);
                                                z_max=z_max+1;
                                                z=z_max;

                                                if(z_max==taille_z-1) {
                                                        //printf("fin de colonne\n");
                                                }
                                        }
                                        if(z_max!=taille_z && z_max-z_min>1 && volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]!=0) { //il faut au moins un trou pour interpoler
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
///#######masque pour &écraser jumeau############""
void antigaussienne(double *tab, int Tx, int sigma, float A, int Exy)
{

        int corr_paire=0;
        if(Tx%2==0) {
                cout<<"taille Tx="<<Tx<<" paire : le masque ne serait pas centré!"<<endl;
                corr_paire=1;
        }

        if(sigma==0)
                sigma=1;
        short  int x,y, Tinf=-round(Tx/2),Tsup=round(Tx/2);
        short unsigned int cptx,cpty;
        float Ex,Ey,sigmax,sigmay;
        Ex=Exy;
        Ey=Exy;
        sigmax=sigma;
        sigmay=sigma;
        if(Tx==1)
                tab[0]=0;
        else {
                for(x=Tinf; x<Tsup+1-corr_paire; x++) {
                        cptx=x+Tsup;
                        for( y=Tinf; y<Tsup+1-corr_paire; y++) {
                                cpty=y+Tsup;
                                tab[cpty*Tx+cptx]=1-A*exp(-(pow((x-Ex),4)/(2*sigmax*sigmax)+pow((y-Ey),4)/(2*sigmay*sigmay)));
                        }
                }
        }

}

//#############################################################"
double *tukey2D(int dimx,int dimy, float alpha)
{
        int N=dimx;
        double *  tuk2D=new double[dimx*dimy];
        double tuk1Dx [dimx];
        double tuk1Dy [dimy];

        int borne1=round(alpha*(N-1)/2);
        int borne2=round((N-1)*(1-alpha/2));

        //memset(tuk2D, 0, dim_entree.x*dim_entree.y*8);
        //memset(tuk1Dx, 0, dim_entree.x*8);
        for(int cpt=0; cpt<borne1+1; cpt++)
                tuk1Dx[cpt]=0.5*(1+cos(3.1415*(2*cpt/(alpha*(N-1))-1)));
        for(int cpt=borne1+1; cpt<borne2+1; cpt++)
                tuk1Dx[cpt]=1;
        for(int cpt=borne2+1; cpt<N; cpt++)
                tuk1Dx[cpt]=0.5*(1+cos(3.1415*(2*cpt/(alpha*(N-1))-2/alpha+1)));


        for(int cpt=0; cpt<N*N; cpt++) {
                int cptx=cpt%(N);
                int cpty=cpt/(N);
                tuk2D[cpt]=tuk1Dx[cptx]*tuk1Dx[cpty];
        }
        return tuk2D;
}

void ecrire_rapport(int NXMAX,float rayon,float Rf, float K, int DIMX_CCD2,int coin_x, int coin_y,short int precision_exportation,char* chemin,int nb_proj,float n1,float NA,float Tp, int G)
{
        time_t date;
        time(&date);
        FILE *fichier_rapport ;
        /*  ouverture pour ecriture (w) en mode texte (t) */
        fichier_rapport = fopen ("/home/mat/tomo_test/rapport_calcul.txt", "wt") ;
        if (fichier_rapport == NULL) {
                printf ("impossible de créer le fichier rapport_calcul.txt\n") ;
                exit (0) ;
        }
        //fprintf(fichier_rapport,"Date     : %s\nNXMAX    : %i\n,Rayon : %i\nRf       : %f\nK : %f\ndimx_ccd : %d\ncoin_x   : %d\ncoin_y   : %d\nPrecision: %i bits\nSession  : %s\nnb_proj  : %i\nindice n1: %f\nNA       : %f\nT_pixel  : %e\nG        : %i",ctime(&date),NXMAX,rayon,Rf,K,DIMX_CCD2,coin_x,coin_y,8*sizeof(precision),chemin,nb_proj,n1,NA,Tp,G);
        fprintf(fichier_rapport,"Date : %s\n NXMAX=%i\n Rf=%f\n Precision: %i bits\n Session  : %s\n Nombre de projections : %i",ctime(&date),NXMAX,Rf,8*precision_exportation,chemin,nb_proj);
        fclose (fichier_rapport);
}


void genereCache(double masque[], int t_image, int t_mask, int centreX, int centreY)
{
        int t_imageX=t_image;
        int t_imageY=t_image;
        if((t_mask%2)!=0)
                cout<<"fonction genere_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;
        for(int pixel=0; pixel<2*t_image*2*t_image; pixel++) {
                masque[pixel]=0;
        }
        for(int x=centreX-t_mask/2; x<centreX+t_mask/2; x++) {
                for(int y=centreY-t_mask/2; y<centreY+t_mask/2; y++) {
                        //le masque repasse du coté gauche lorsque le cache touche le bord droit! A corriger (?)
                        int x_2=x;
                        int y_2=y;
                        if(x>2*t_imageX)
                                x_2=x-2*t_imageX;
                        if(y>2*t_imageY)
                                y_2=y-2*t_imageY;
                        if(x<0)
                                x_2=2*t_imageX+x;
                        if(y<0)
                                y_2=2*t_imageY+y;
                        ///coordonnées 1D du centre
                        int cptj=2*t_imageY*y_2+x_2;
                        masque[cptj]=1;
                }
        }
}

/////////////////////////////////////////////////////////////////////////////////////
//fonction chargeant une image 2D. Retourne un pointeur sur le tableau 1D qui la contient
//  Nécessite le nom du pointeur à remplir,
//le numéro du pas de phase shifting (1,2 3 ou 4), le chemin, et cpt_fichier
//pour boucler sur les 1000 projections
/////////////////////////////////////////////////////////////////////////////////////


void charger_image2D(unsigned char* phasei, char *imgFile,int coin_x,int coin_y,int taille_x,int taille_y)
{
        Var2D taille {taille_x,taille_y},coin {coin_x,coin_y};
        //char CopieChemin[100];
        //char FinNom[15];
        //strcpy(CopieChemin,Chemin);//sauvegarder chemin car sprintf écrase
        //sprintf(FinNom,"-00%i.bmp",numero);
        //créer le chemin final
        //strcat(CopieChemin, FinNom);
        //printf("Copiechemin vaut : %s\n", CopieChemin);
        //charger l'image dans le tableau
        ///tester la validité du chemin "copiechemin"
        /*FILE *test_chemin = fopen(CopieChemin,"r");

          if (test_chemin == NULL)
          {
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
        rempli_tableau(phasei, imgFile, coin,taille);
        //return phasei;
}

/////////////////////////////////////////////////////////////////////////////////////
//Fonction transferant les pixels de l'image 2D de taille (nx,ny) dans un tableau 1D de taille (nx*ny)
// nécessite la bibliotheque imageMagick
/////////////////////////////////////////////////////////////////////////////////////

void rempli_tableau(unsigned char *finalArray, string path, Var2D coin, Var2D taille)
{
        Image Monimage;
        //int i, currentImageWidth, currentImageHeight;
        Monimage.read(path);////// chargement en memoire de l'image
        Monimage.type( GrayscaleType );	////// Mon image est N&B
        //Monimage.display();
        //Monimage.channelDepth();
        Monimage.compressType(NoCompression);
        //Monimage.crop(Geometry(dimx_ccd,dimy_ccd, 0, 0) );
        //Monimage.write("/home/mat/tomo_test/test.bmp");
        //currentImageWidth = Monimage.columns();////// extraction des caractéristiques de l'image
        //currentImageHeight = Monimage.rows();
        //	finalArray = new unsigned char[taille_x*taille_y];////// reservation de la taille memoire du tableau final
        ////// lecture de l'image
        Monimage.getPixels(coin.x,coin.y,taille.x,taille.y);
        ////// ecriture de l'image dans un tableau
        Monimage.writePixels(GrayQuantum, finalArray);
}

/////////////////////////////////////////////////////////////////////////////////////
//fonction de calcul TF2D
/////////////////////////////////////////////////////////////////////////////////////
void TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y)
{
        fftw_plan_with_nthreads(3);
        int N=taille_x*taille_y;
        fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
        fftw_plan p;
        //Réservation memoire
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        //Récupération de l'image dans la partie reelle de l'entree
        for(int cpt=0; cpt<N; cpt++) {
                in[cpt][0]=entree_reelle[cpt];
                in[cpt][1]=entree_imag[cpt];
        }
        //calcul du plan, parametre servant a calculer et optimiser le FFT
        p=fftw_plan_dft_2d( taille_x,  taille_y, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p); /* repeat as needed */

        for(int cpt=0; cpt<(N); cpt++) {

                fft_reel[cpt]=out[cpt][0]; //division par N^2 pour normaliser la fftw qui n'est pas normalisée
                fft_imag[cpt]=out[cpt][1];
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
}

void TF2D_INV(double entree_reelle[],double entree_imag[],double sortie_reelle[],double sortie_imag[],int taille_x,int taille_y)
{
        fftw_plan_with_nthreads(3);
        int N=taille_x*taille_y;
        fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
        fftw_plan p;
        //Réservation memoire
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        //Récupération de l'image dans la partie reelle de l'entree
        for(int cpt=0; cpt<N; cpt++) {
                in[cpt][0]=entree_reelle[cpt];
                in[cpt][1]=entree_imag[cpt];
        }
        //calcul du plan, parametre servant a calculer et optimiser le FFT
        p=fftw_plan_dft_2d( taille_x,  taille_y, in, out,FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(p); /* repeat as needed */

        for(int cpt=0; cpt<(N); cpt++) {

                sortie_reelle[cpt]=out[cpt][0];
                sortie_imag[cpt]=out[cpt][1];
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
}

/////////////////////////////////////////////////////////////////////////////////////
//fonction de circshift

/////////////////////////////////////////////////////////////////////////////////////
void circshift2(double* entree, double* result, Var2D dim,Var2D decal)
{
        //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
        decal.y=decal.y%dim.y;
        decal.x=decal.x%dim.x;

        for(int yi=0; yi<dim.y; yi++) {
                int yi_decal=yi+decal.y;
                //si dépassement des dimensions
                if(yi+decal.y>dim.y-1)
                        yi_decal=-dim.y+yi+decal.y;
                //si décage négatif
                if(yi+decal.y<0)
                        yi_decal=dim.y+yi+decal.y;
                int hauteur_decal=yi_decal*dim.x;
                int hauteur=yi*dim.x;
                for(int xi=0; xi<dim.x; xi++) {
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
void circshift3(double* entree, double* result, Var2D dim,Var2D decal)
{
        //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
        decal.y=decal.y%dim.y;
        decal.x=decal.x%dim.x;

        for(int yi=0; yi<decal.y; yi++) {
                for(int xi=0; xi<decal.x; xi++)

                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
                        //1er quadrant vers 4 eme
                        result[pixel_shift]=entree[pixel];
                        //4 eme quadrant vers 1er
                        result[pixel]=entree[pixel_shift];
                        //2eme vers 3eme
                        result[(yi+decal.y)*dim.x+xi]=entree[pixel+decal.x];
                        //3eme vers 2eme
                        result[pixel+decal.x]=entree[(yi+decal.y)*dim.x+xi];
                }
        }
}


double* circshift(double* entree, int dimx,int dimy,int decal_x,int decal_y)
{
        double * entree_shift=new double[dimx*dimy];//fuite memoire!
        for(int yi=0; yi<decal_y; yi++) {
                for(int xi=0; xi<decal_x; xi++)

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

}
void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D)
{
        decal3D.x=decal3D.x%dimFinal3D.x;//élmiiner les "modulos"
        decal3D.y=decal3D.y%dimFinal3D.y;
        decal3D.z=decal3D.z%dimFinal3D.y;

        unsigned short int xi,yi,zi=0;
        short int x2,y2,z2=0; //signé car une fois décalé, peuvent être négatifs!
        const unsigned int taille_plan =dimFinal3D.x*dimFinal3D.y;

        for(zi=0; zi<dimFinal3D.z; zi++) {
                if(zi+decal3D.z>dimFinal3D.z-1) { //dépassement à droite
                        z2=zi+decal3D.z-dimFinal3D.z;
                } else {
                        if(zi+decal3D.z<0) { //dépassement à gauche
                                z2=dimFinal3D.z+(decal3D.z+zi);
                        } else {
                                z2=zi+decal3D.z;
                        }
                }
                int nb_pixelz_decal=z2*taille_plan;
                unsigned int nb_pixelz=zi*taille_plan;
                for(yi=0; yi<dimFinal3D.y; yi++) {
                        if(yi+decal3D.y>dimFinal3D.y-1) { //dépassement à droite
                                y2=yi+decal3D.y-dimFinal3D.y;
                        } else {
                                if(yi+decal3D.y<0) { //dépassement à gauche
                                        y2=dimFinal3D.y+(decal3D.y+yi);
                                } else {
                                        y2=yi+decal3D.y;
                                }
                        }
                        int nb_lignes=yi*dimFinal3D.x;
                        int nb_lignes_decal=y2*dimFinal3D.x;

                        for(xi=0; xi<dimFinal3D.x; xi++) {
                                if(xi+decal3D.x>dimFinal3D.x-1) { //dépassement à droite
                                        x2=xi+decal3D.x-dimFinal3D.x;
                                } else {
                                        if(xi+decal3D.x<0) { //dépassement à gauche
                                                x2=dimFinal3D.x+(decal3D.x+xi);
                                        } else {
                                                x2=xi+decal3D.x;
                                        }
                                }

                                // pixel=zi*dimFinal3D.x*dimFinal3D.y+yi*dimFinal3D.x+xi;
                                //pixel_decal=z2*dimFinal3D.x*dimFinal3D.y+y2*dimFinal3D.x+x2;
                                //volume3D_shift[pixel_decal]=volume3D[pixel];
                                //volume3D_shift[z2*dimFinal3D.x*dimFinal3D.y+y2*dimFinal3D.x+x2]=volume3D[zi*dimFinal3D.x*dimFinal3D.y+yi*dimFinal3D.x+xi];
                                volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2]=volume3D[nb_pixelz+nb_lignes+xi];
                                //memcpy(volume3D_shift+(nb_pixelz_decal+nb_lignes_decal+x2), volume3D+(nb_pixelz+nb_lignes+xi), 8); commande ok mais plus lent !?

                        }
                }
        }
}
///fonction de masquage
void multiplier_masque(double image[], unsigned char masque[], int t_image, int t_mask, int centreX, int centreY)
{
        int t_imageX=t_image;
        int t_imageY=t_image;
        //if((t_mask%2)!=0)
        //cout<<"fonction multiplier_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;

        for(int x=centreX-t_mask/2; x<centreX+t_mask/2; x++) {
                for(int y=centreY-t_mask/2; y<centreY+t_mask/2; y++) {
                        //si le masque déborde de l'image (attention le masque repasse du coté gauche lorsque le jumeau touche le bord droit! A corriger
                        int x_2=x;
                        int y_2=y;
                        if(x>2*t_imageX)
                                x_2=x-2*t_imageX;
                        if(y>2*t_imageY)
                                y_2=y-2*t_imageY;
                        if(x<0)
                                x_2=2*t_imageX+x;
                        if(y<0)
                                y_2=2*t_imageY+y;
                        ///coordonnées 1D du centre
                        int cptj=2*t_imageY*y_2+x_2;
                        //if((double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)!=1)
                        ///cout<<"masque:"<<(double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)<<endl;
                        image[cptj]=image[cptj]*masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)]/255;
                }
        }
}

///fonction de masquage
void multiplier_masque2(double image[], double masque[], int t_image, int t_mask, int centreX, int centreY)
{
        int t_imageX=t_image;
        int t_imageY=t_image;
        //if((t_mask%2)!=0)
        //cout<<"fonction multiplier_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;

        for(int x=centreX-t_mask/2; x<centreX+t_mask/2; x++) {
                for(int y=centreY-t_mask/2; y<centreY+t_mask/2; y++) {
                        //si le masque déborde de l'image (attention le masque repasse du coté gauche lorsque le jumeau touche le bord droit! A corriger
                        int x_2=x;
                        int y_2=y;
                        if(x>2*t_imageX)
                                x_2=x-2*t_imageX;
                        if(y>2*t_imageY)
                                y_2=y-2*t_imageY;
                        if(x<0)
                                x_2=2*t_imageX+x;
                        if(y<0)
                                y_2=2*t_imageY+y;
                        ///coordonnées 1D du centre
                        int cptj=2*t_imageY*y_2+x_2;
                        //if((double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)!=1)
                        ///cout<<"masque:"<<(double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)<<endl;
                        image[cptj]=image[cptj]*masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)];
                }
        }
}



// soient deux images allouées en row-major
// découpe une fenêtre de img_src et l'écrit dans img_dst
// il faut préciser la taille des deux images et la position du coin supérieur gauche de la fenêtre

void crop_window(double* img_src, size_t src_dimx, size_t src_dimy, double* img_dst, size_t dst_dimx, size_t dst_dimy, size_t ul_corner_x, size_t ul_corner_y)
{
        assert(dst_dimx < src_dimx);
        assert(dst_dimy < src_dimy);

        assert((ul_corner_x + dst_dimx) <= src_dimx);
        assert((ul_corner_y + dst_dimy) <= src_dimy);

        double *start_src = img_src + (ul_corner_y * src_dimx);

        double *curs_src = start_src;
        size_t pos_dst = 0;
        size_t i, j;

        for (i = 0; i < dst_dimy; i++) {
                for (j = ul_corner_x; j < ul_corner_x + dst_dimx; j++) {
                        img_dst[pos_dst] = curs_src[j];
                        pos_dst++;
                }
                curs_src += src_dimx;
        }
}


