#include "projet.h"
#include "fonctions.h"
//#include "CImg.h"


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

//using namespace cimg_library;
using namespace std;
using namespace Magick;




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
{        // ----------------------------------------------------------------------
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
        float precision;//pour conversion de type en flottant (32 bits) lors de l'écriture finale
        const short int precision_exportation=sizeof(precision);
        cout<<"precision_export="<<precision_exportation<<endl;

        int points_faux=0;

        ///constantes multithread FFTW
        int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        int nthreads=3;
        printf("fftwThreadInit: %i\n",fftwThreadInit);
        ///Chemin fichiers-------------------------------------
        char* Dossier_acquiz= input_dir; //dossier des données
        char NumSession[]="i0";//partie fixe du nom du fichier
        char NomHoloPart1[strlen(Dossier_acquiz)+strlen("i0xxx")]; //réserver la place pour un nom montant au max à l'hologramme 999
        concatener(Dossier_acquiz,NumSession,NomHoloPart1);//concaténer le dossier avec la date, partie fixe au début du nom
        //chemin pour référence d'amplitude
        char CheminModRef[strlen(Dossier_acquiz)+strlen("mod_ref.pgm")+1];
        concatener(Dossier_acquiz,"/mod_ref.pgm",CheminModRef);

        const int premier_plan=1, Num_Angle_final= LAST_ANGLE_NUM,//
        NbAngle=Num_Angle_final-premier_plan+1, SautAngle=1, NbPixROI2d=DIMX_CCD2*DIMY_CCD2;//nombre de pixel (pour le tableau 2D)

        ///-------CONSTANTES EXPERIMENTALES->PREPROC?
        const float
        n1=1.515,	//indice de l'huile
        NA=1.40,	//ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
        theta=asin(NA/n1),//theta_max
        lambda=632*pow(10,-9), //longueur d'onde
        G=100,	//grossissement telan+objectif

        Tp=8*pow(10,-6),//cam. photon focus
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

        const int
        NXMAX=round(n1*Rf*Tps*DIMX_CCD2/(G*lambda)*NA/n1),//NXMAX=round(rayon*NA/n1)=round(rayon*sin theta_max);
              //Variable pour (éventuellement) Forcer la dimension de l'espace à une certaine taille, attention, ceci changera la taille des pixels tomo
        dim_final=512;//peut être différent de 4*NXMAX, mais l'image final sera zoomée;
        //
        double zoom=double(dim_final)/double(4*NXMAX);
        cout<<"zoom="<<zoom<<endl;

        //NXMAX_Rf=round(zoom*n1*Tps*DIMX_CCD2/(G*lambda)*NA/n1); //dimension espace final
        //cout<<"NXMAX_Rf="<<NXMAX_Rf<<endl;

        const int NYMAX=NXMAX;
        // Définition de la fenêtre de sélection centrée sur (cx, cy)
        size_t upperleft_x = circle_cx - (2*NXMAX / 2);
        size_t upperleft_y = circle_cy - (2*NYMAX / 2);
        assert (upperleft_x < WINDOW_X);
        assert (upperleft_y < WINDOW_Y);
        // création d'une coupe 2D centrée sur le cercle
        double *decoupe=new double[2*NXMAX*2*NYMAX];
        int centres_exclus=0, nb_proj=0;
        float tailleTheoPixelTomo=tailleTheoPixelHolo*DIMX_CCD2/(4*NXMAX);

        cout<<"Taille théorique pixel tomo:"<<tailleTheoPixelTomo<<endl;

        if(zoom!=1) {
                cout<<"dim_final forcée à "<<dim_final<<" au lieu de "<<4*NXMAX<<"-->facteur de zoom="<<zoom<<endl;
                cout<<"nouvelle taille pixel tomo imposée="<<tailleTheoPixelTomo/zoom<<"nm"<<endl;
                cout<<"###########################################"<<endl;
        }

//        printf("Précision d'exportation: %i bits\n chemin :%s \n",8*sizeof(precision),Chemin);
        //cerr << endl << "dv0:" << round(4 * NXMAX) << " NXMAX:" << NXMAX;

        ///---------------Calcul de quelques constantes, afin d'éviter de les calculer dans la boucle--------------

        const int dimVolX=round(dim_final),
                        dv0s2rf=round(dim_final/2),
                        dv1s2rf=round(dim_final/2),
                        dv1xdimVolX=round(dim_final*dim_final),
                        dv2s2rf=round(dim_final/2);
        printf("dimx espace 3D, dimVolX: %i\n",dimVolX);
        fflush(stdout);
        cout.flush();

        /// //////////////////////////////////////////////////fin calcul constante//////////////////////////

        ///------------réservation des plan images/Tf, ROI découpée sur la caméra------------------------------------------
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
        double  *holo_Re_shift=new double[DIMX_CCD2*DIMY_CCD2];
        double  *holo_Im_shift=new double[DIMX_CCD2*DIMY_CCD2];

        double *fft_reel_tmp_centre=new double[DIMX_CCD2*DIMY_CCD2];
        double *fft_imag_tmp_centre=new double[DIMX_CCD2*DIMY_CCD2];

        ///---------------Espace 2D Après Découpe  à NXMAX;-----------------------------------------------------
        double *fft_reel=new double[4*NXMAX*NYMAX];
        double *fft_imag=new double[4*NXMAX*NYMAX];
        double *fft_reel_norm=new double[4*NXMAX*NYMAX];
        double *fft_imag_norm=new double[4*NXMAX*NYMAX];
        double *fft_module=new double[4*NXMAX*NYMAX];
        double *fft_module_shift=new double[4*NXMAX*NYMAX];
        double *fft_reel_shift=new double[4*NXMAX*NYMAX];
        double *fft_imag_shift=new double[4*NXMAX*NYMAX];

        //TF normalisée en amplitude cplx
        double  *fft_reel_shift_norm=new double[4*NXMAX*NYMAX];
        double  *fft_imag_shift_norm=new double[4*NXMAX*NYMAX];
        int *centre=new int[4*NXMAX*NYMAX];//pour mettre la position des centres translatés, on crée une variable 2D de la taille d'un plan apres tomo

        ///---------espace3D final --------------------------------------------------
        const int N_tab=dim_final*dim_final*dim_final;//Nombre de pixel dans l'esapce final (indice max du tableau 3D+1)
        printf("N_tab (indice espace 3D)= %i \n",N_tab);
        double *reel_arc=new double[N_tab];//partie reelle du papillon dans l'espace reciproque
        double *imag_arc=new double[N_tab];//partie imaginaire du papillon dans l'espace reciproque
        double *sup_redon=new double[N_tab];//pour normaliser par rapport au nombre de fois où l'on remplit la frequence

        ///-----Réservation chrono----------------
        clock_t
        temps_depart, /*temps de départ de programme */
        temps_initial, /* temps initial en micro-secondes */
        temps_final, /* temps d'arrivée du programme */
        temps_arrivee;   /* temps final en micro-secondes */
        float temps_cpu=0,     /* temps total en secondes */
        temps_total=0;
        temps_depart = clock();
        temps_initial = clock ();

        ///Chargement de la réference en module et calcul amplitude
        unsigned char* mod_ref=new unsigned char[NbPixROI2d]; //reservation reference et "noir" camera pour une "tentative" de correction de la ref
        rempli_tableau(mod_ref,CheminModRef,coin,dimCCD);
        for(int pixel=0; pixel<NbPixROI2d; pixel++) {
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
//        Var2D dim2D= {2*NXMAX,2*NYMAX},decal2D= {NXMAX,NYMAX};

        ///début de la grosse boucle sur les angles/////////////////////////
        //##############################################################################################################################

        char TamponChemin[ strlen(input_dir) + 256 ];
        for(int cpt_angle=premier_plan; cpt_angle<premier_plan+NbAngle; cpt_angle=cpt_angle+SautAngle) { //boucle sur tous les angles
                //cout<<"cpt_angle : "<<cpt_angle<<endl;
                if((cpt_angle-100*(cpt_angle/100))==0)
                        printf("cpt_angle=%i\n",cpt_angle);

                ///Concaténer le numéro d'angle dans le nom
                sprintf(TamponChemin, "%s/i%04i_ha.pgm", input_dir, cpt_angle);
                test_existence = fopen(TamponChemin, "rb");

                if(test_existence!=NULL) {
                        fclose(test_existence);
                       // fprintf(stdout, "%s\n", TamponChemin);
                        //continue;
                        ///charger l'hologramme hors axe
                        charger_image2D(holo1,TamponChemin, coin_x, coin_y,DIMX_CCD2,DIMY_CCD2);
                        for(int pixel=0; pixel<NbPixROI2d; pixel++) {
                                ampli_ref[pixel]=sqrt((double)mod_ref[pixel]);
                                holo_Re[pixel]=(double)holo1[pixel];
                        }
                        //////////////////////////// Circshift avant TF2D
                        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        circshift2(holo_Re,holo_Re_shift, dimCCD,decalCCD);
                        memset(holo_Im_shift, 0, (dimCCD.x*dimCCD.y)*sizeof(double));	  //pas besoin de décaler la partie imag car pleine de zéro.
                        //////////////////////////TF 2D du front d'onde: passage du plan image au plan réciproque////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                       /* FILE *fichier_holo2d;
                        fichier_holo2d= fopen("/home/mat/tomo_test/TF2d/holo2d", "a+b");
                        if(fichier_holo2d==0)
                                cout<<"Erreur d'ouverture du fichier TF2d"<<endl;
                        for(int cpt=0; cpt<(dimCCD.x*dimCCD.y); cpt++) {
                                //Ecriture sur le disque
                                fwrite(&holo_Re_shift[cpt],sizeof(double),1,fichier_holo2d);
                        }
                        fclose(fichier_holo2d);*/

                        TF2D(holo_Re_shift,holo_Im_shift,fft_reel_tmp,fft_imag_tmp,dimCCD.x,dimCCD.y);

                        //Décalage  sur fft_reel_tmp, pour recentrer le spectre avant découpe (pas obligatoire mais plus clair)
                        circshift2(fft_reel_tmp, fft_reel_tmp_centre, dimCCD,decalCCD);
                        circshift2(fft_imag_tmp, fft_imag_tmp_centre, dimCCD,decalCCD);

                        crop_window(fft_reel_tmp_centre, 1024, 1024, decoupe, 2*NXMAX, 2*NYMAX, upperleft_x, upperleft_y);
                        // copie de la couche là où il faudrait
                        for (size_t s = 0; s < 2*NXMAX * 2*NXMAX; s++)
                                fft_reel[s] = decoupe[s];

                        crop_window(fft_imag_tmp_centre, 1024, 1024, decoupe, 2*NXMAX, 2*NYMAX, upperleft_x, upperleft_y);
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

                        ///----------------------------Recherche du maximum dans ref non centré-----------------------------------------

                        int cpt_max=0;
                        fft_module[0]=pow(fft_reel[0],2)+pow(fft_imag[0],2);

                        //Recherche du MAX dansle module
                        for(int cpt=1; cpt<(4*NXMAX*NYMAX); cpt++) {
                                fft_module[cpt]=pow(fft_reel[cpt],2)+pow(fft_imag[cpt],2);
                                if(fft_module[cpt]>fft_module[cpt_max]) {
                                        cpt_max=cpt;
                                }
                        }
                      /*  FILE *fichier_TF2d_crop;
                        fichier_TF2d_crop= fopen("/home/mat/tomo_test/TF2d/tf_crop", "a+b");
                        if(fichier_TF2d_crop==0)
                                cout<<"Erreur d'ouverture du fichier tf_crop"<<endl;
                        for(int cpt=0; cpt<(2*NXMAX*2*NXMAX); cpt++) {
                                //Ecriture sur le disque
                                fwrite(&fft_reel[cpt],sizeof(double),1,fichier_TF2d_crop);
                        }
                        fclose(fichier_TF2d_crop);*/

                        double max_part_reel = fft_reel[cpt_max];
                        double max_part_imag = fft_imag[cpt_max];
                        double max_module = fft_module[cpt_max];

                        //printf("max_part_reel : %f, max_part_imag: %f,  max_module %f \n ",max_part_reel, max_part_imag, max_module);

                        ///Coordonnées dans l'espace 2D à partir de l'indice 1D: (xc,yc)=(0,0)=en haut à gauche de l'image
                        int xmi=cpt_max%(2*NXMAX);
                        int ymi=cpt_max/(2*NYMAX);
                      //  cout<<"xmi="<<xmi<<"//ymi="<<ymi<<endl;
                        int xc=xmi-NXMAX;
                        int yc=ymi-NYMAX;
                        ///sauver tf  avant normalisation
                      /*  FILE* fichier_fft_reel;
                        fichier_fft_reel = fopen("/home/mat/tomo_test/TF2d_apres_masquage/fft_reel.bin", "a+b");//a+b=mode incrémental=tout écrire dans le même fichier
                        for(int cpt=0; cpt<4*NXMAX*NYMAX; cpt++) {
                                precision=fft_reel[cpt];
                                fwrite(&precision,sizeof(precision),1,fichier_fft_reel);
                        }
                        fclose(fichier_fft_reel);*/
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
                       // printf("xm0 %i \n", xm0);
                        int ym0=(ymi-NYMAX);
                       //      printf("ym0 %i \n", ym0);
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


                                double zm0;

                                //printf("round(rayon*rayon-(xm0)^2-(ym0)^2: %i\n",round(rayon*rayon-(xm0)^2-(ym0)^2));

                                double zm0_carre = rayon*rayon-xm0*xm0-ym0*ym0;
                                if(round(zm0_carre)>-1) {

                                        zm0=sqrt(zm0_carre);
                                        nb_proj++;
                                        temps_initial = clock ();
                                        int NXMAX_CARRE=NXMAX*NXMAX;

                                        for (int y = -NYMAX; y < NYMAX; y++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
                                                int y_carre=y*y;
                                                for (int x = -NXMAX; x < NXMAX; x++) { //on balaye l'image 2D en y, centre au milieu
                                                        int cpt=(y+NYMAX)*2*NXMAX+x+NXMAX;//calcul du cpt du tableau 1D de l'image 2D
                                                        if(x*x+y_carre<NXMAX_CARRE) { //ne pas depasser l'ouverture numérique pour 1 hologramme
                                                                double z_carre=r2-x*x-y_carre; //altitude au carré des données
                                                                double z=round(sqrt(z_carre)-zm0);
                                                                double altitude=(z+dv2s2rf)*dv1xdimVolX; //donne n'importequoi sans l'arrondi sur z!!

                                                                k=(-xm0+x+dv0s2rf)+(-ym0+y+dv1s2rf)*dimVolX+round(altitude);//indice du tableau 1D du volume 3D
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

        delete[] fft_reel_tmp_centre;
        delete[] fft_imag_tmp_centre;
        delete[] fft_reel_shift;
        delete[] fft_imag_shift;
        delete[] fft_module_shift;
        delete[] fft_reel_shift_norm;
        delete[] fft_imag_shift_norm;
        delete[] fft_reel_norm;
        delete[] fft_imag_norm;
        delete[] fft_module;
        delete[] holo_Im_shift;
        delete[] holo_Re_shift;
        delete[] masque;
        delete[] cache_jumeau;
        delete[] holo1;
        delete[] holo_Re;

        ////////////////////////////////////////////////////////////////////
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

        delete[] centre;

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


        /// ///////////////////////////mesure du temps de calcul

        printf("nb_proj: %i\n",nb_proj);
        printf("Nombre de centre exclus : %i\n",centres_exclus);
        temps_final = clock ();
        temps_cpu = (temps_final - temps_initial) * 1e-6/nb_proj;
        printf("temps moyen pour 1 angle: %f\n",temps_cpu);
        temps_cpu = (temps_final - temps_initial) * 1e-6;
        printf("temps total pour %i angle(s): %f\n", nb_proj, temps_cpu);
        temps_initial=clock();//enclenchement du chronometre


        ///------------------normalisation par sup_redon------------------------------------------------
        for(int cpt=0; cpt<N_tab; cpt++) {
                if (sup_redon[cpt]==0) { ////////////////remplace les 0 de sup_redon par des 1---> evite la division par 0
                        //printf("sup_redon= %f \n" , sup_redon[cpt]);
                        sup_redon[cpt]=1;
                }
                reel_arc[cpt] = reel_arc[cpt]/sup_redon[cpt];//moyennage par sup_redon
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
        //interp3D(reel_arc,  dimVolX,dimVolX, dimVolX);
        //interp3D(imag_arc,  dimVolX,dimVolX, dimVolX);

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
        decal3D={dimVolX/2,dimVolX/2,dimVolX/2},

                 dimFinal= {dimVolX,dimVolX,dimVolX};
                        cout<<"decal3D.x"<<decal3D.x<<endl;
        circshift3D2(reel_arc, reel_arc_shift, dimFinal,decal3D);//ancienne valeur, causant des problèmes d'arrondis
        circshift3D2(imag_arc, imag_arc_shift, dimFinal,decal3D);

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

        p3D=fftw_plan_dft_3d(dimVolX, dimVolX, dimVolX, in3D, out3D,FFTW_BACKWARD, FFTW_ESTIMATE);

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

        //circshift3D(final_reel, final_reel_shift, dimVolX, dimVolX, dimVolX);
        //circshift3D(final_imag, final_imag_shift, dimVolX, dimVolX, dimVolX);
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

        //ecrire_rapport(NXMAX,rayon,Rf,K,DIMX_CCD2,coin_x,coin_y,precision_exportation,Fichier_holo,nb_proj,n1,NA,Tp,G);
        ///---Ecriture des 3 volumes.
        for(int cpt=0; cpt<(N3D); cpt++) {
                precision=final_reel_shift[cpt];
                fwrite(&precision,sizeof(precision),1,fichier_final_reel_shift);

                precision=final_imag_shift[cpt];
                fwrite(&precision,sizeof(precision),1,fichier_final_imag_shift);

                precision=valeur_module_shift[cpt];
                fwrite(&precision,sizeof(precision),1,fichier_final_modul_shift);//ecriture du module
        }

        delete[] ampli_ref;
        delete[] mod_ref;
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

