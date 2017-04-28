#include "projet.h"
#include "fonctions.h"
#include "deroulement.h"
#include "Correction_aberration.h"
//#include "CImg.h"
#include <chrono>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <octave-2.9.9/octave/oct.h>
#include <time.h>
#include <ctime>
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
                printf("-o repertoire de sortie");

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
     /*   usage(argc, argv);


        char* input_dir = NULL;
        char* output_dir= NULL;
        size_t circle_cx = 0, circle_cy = 0;

       while (argc > 0) {//chemin acquiz
                if (!strcmp(argv[0], "-i") && (argc > 1)) {
                        input_dir = argv[1];
                        argc--;
                        argv++;
                }
                if (!strcmp(argv[0], "-o") && (argc > 1)) {
                        output_dir = argv[1];
                        argc--;
                        argv++;
                }
              /*  if (!strcmp(argv[0], "-c") && (argc > 2)) {
                        circle_cx = atoi(argv[1]);
                        circle_cy = atoi(argv[2]);

                        argc-=2;
                        argv+=2;
                }

                argc--;
                argv++;
        }*/
        string home=getenv("HOME");
        string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
        string chemin_config_GUI=getenv("HOME")+fin_chemin_gui_tomo;
        string repertoire_config=extract_string("CHEMIN_CONFIG",home+fin_chemin_gui_tomo);
        string chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);
        string chemin_acquis=extract_string("CHEMIN_ACQUIS",home+fin_chemin_gui_tomo);
//        cout<<"output_dir="<<output_dir<<endl;

        string tampon="UBornfinal_Im.raw";
        string result=chemin_result+tampon;

        //Effacer précédent résultats.
        //cout<<"bazar appppend"<<chemin_result.append("Ubornfinal_Im.raw").c_str()<<endl;
        if( remove(result.c_str())!= 0 )
        perror( "Fichier Ubornfinal_Im inexistant" );

        tampon="UBornfinal_Re.raw";
        result=chemin_result+tampon;
        if( remove(result.c_str()) != 0 )
        perror( "Fichier UBornfinal_Re inexistant" );

       /* char* pPath_home;
        pPath_home = getenv ("HOME");
        if (pPath_home==NULL)
        printf ("le repertoire %s n'existe pas",pPath_home);
        return 0;*/

        ///attention, écrase les arguments d'entrée
        // input_dir="/opt/resultats2015/Ludo6/ACQUIS/";
        //input_dir="/opt/resultats2015/ACQUIS_bbr_pollen/";
        //input_dir="/opt/resultats2015/bille_4V_avec_mod_ref/ACQUIS/";
       // string fic_cfg_recon="/opt/Acquiz_tomo/2016/pollen_marguerite/recon.txt";
       // string fic_cfg_recon="/opt/Acquiz_tomo/2016/pollen_marguerite/deg_0/recon.txt";
        // input_dir="/opt/Acquiz_tomo/2016/pollen_marguerite/deg_0/";
        //string fic_cfg_recon="/opt/resultat2015/bertrand_zeolithe_20151110/recon.txt";
        //input_dir="/opt/Acquiz_tomo/resultats2015/bille_4_2_avec_blanc_cplx/bille2_4.2V/ACQUIS/";
       // input_dir="/opt/Acquiz_tomo/2016/bille_62micron_3112016/blanc/";
        //string home=pPath_home;
       // string repertoire=input_dir;
       // string rep_sortie=output_dir;
        string fic_cfg_recon=chemin_acquis+"recon.txt";
        cout<<"fichier recon"<<fic_cfg_recon<<endl;
       // ecrire_val("BORN",1975,fic_cfg_recon);
        string fic_cfg_manip=repertoire_config+"config_manip.txt";
        cout<<"chemin config"<<fic_cfg_manip<<endl;



        size_t circle_cx=extract_val("CIRCLE_CX",fic_cfg_manip);
        size_t circle_cy=extract_val("CIRCLE_CY",fic_cfg_manip);


        assert(chemin_acquis.c_str());
       // assert(circle_cx && circle_cy);
        cout << endl << "input dir: " << chemin_acquis;
        cout << endl << "Hors-axe centre sur : (" << circle_cx << " , " << circle_cy << ")";
        cout.flush();

        //ROI sur camera, a priori : 1024
        Var2D dimROI= {WINDOW_X, WINDOW_Y},decalROI= {dimROI.x/2,dimROI.y/2};
        //valeur du coin pour la découpe
        //découpe centrée - la découpe sera réalisée à la lecture du fichier image
        Var2D coinROI { (IMAGE_DIMX - WINDOW_X)/2 , (IMAGE_DIMY - WINDOW_Y)/2 }; //ROi centré dans l'image CCD.

        //   if(coin_x+DIMX_CCD2>IMAGE_X || coin_y+DIMY_CCD2> IMAGE_Y)//
        //     printf("Attention, la zone TF_UBorne sort de l'image\n");
        float precision;//pour conversion de type en flottant (32 bits) lors de l'écriture finale
        int points_faux=0;
        int centres_exclus=0;  //centre maximum
        ///--------------constantes multithread FFTW
            int fftwThreadInit;
            fftwThreadInit=fftw_init_threads();
            int nthreads=4;
            printf("fftwThreadInit bool: %i\n",fftwThreadInit);
        ///---Chemin fichiers-------------------------------------
            string Dossier_acquiz= chemin_acquis,Numsession="i0",NomHoloPart1=Dossier_acquiz+Numsession;
            string CheminModRef=Dossier_acquiz+"mod_ref.bmp";
            string Chemin_mask=Dossier_acquiz+"Image_mask.pgm";
            cout<<"CheminModRef : "<<CheminModRef<<endl;


        ///--selections des hologrammes par numéro
            const int   premier_plan=extract_val("FIRST_LANE",fic_cfg_recon),
                        Num_Angle_final=extract_val("FINAL_ANGLE",fic_cfg_recon),//
            NbAngle=Num_Angle_final-premier_plan, SautAngle=1, NbPixROI2d=dimROI.x*dimROI.y;//nombre de pixel (pour le tableau 2D)
            cout<<"Nb Angle="<<NbAngle<<endl;

            ///-------CONSTANTES EXPERIMENTALES->PREPROC?
            cout<<"##################### INFO MANIP ##################"<<endl;
            const float
            n0=extract_val("N0",fic_cfg_manip),	//indice de l'huile
            NA=extract_val("NA",fic_cfg_manip),	//ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
            theta=asin(NA/n0),//theta_max
            lambda=extract_val("LAMBDA",fic_cfg_manip), //longueur d'onde
            kv=2*3.1416/lambda,
            f_tube=extract_val("F_TUBE",fic_cfg_manip), ///focale lentille tube
            fobj=extract_val("F_OBJ",fic_cfg_manip),///focale objectif
            G=f_tube/fobj,	//grossissement telan+objectif

            TpCam=extract_val("TPCAM",fic_cfg_manip),//cam. photon focus
            Rf=extract_val("RF",fic_cfg_manip),//1.2;=1/facteur grossissement
            Gt=G/Rf,
            K=lambda*G/(2*NA*TpCam*Rf),// Facteur d'echelle
            Tps=TpCam,//*K; //Taille pixel apres mise à l'echelle important si RF<>1
            tailleTheoPixelHolo=TpCam/Gt*pow(10,9);//Pour info, taille des pixels sur un hologramme=Tpcam/GT
            cout<<"n0="<<n0<<endl;

            //cout<<"taille théorique pixel holo : "<< tailleTheoPixelHolo<<" nm"<<endl;
            ///##########Calcul Rayon théo, fréquence max+ facteur zoom, taille pixel,

            float
            //rayon=n0*Rf*TpCam*DIMX_CCD2/(G*lambda), //Rayon
            Rth=round(dimROI.x*TpCam*n0/(lambda*Gt)); //Rayon théorique

            const int
            NXMAX=extract_val("NXMAX",fic_cfg_manip),

            //NXMAX=round(n0*Rf*TpCam*dimROI.x/(Gt*lambda)*NA/n0),NYMAX=NXMAX,//NXMAX=round(rayon*NA/n0)=round(rayon*sin theta_max);

            NXMAX_Th=round(n0*TpCam*dimROI.x/(Gt*lambda)*NA/n0),
            NYMAX=NXMAX;
            cout<<"NXMAX theorique="<<NXMAX_Th<<endl;
            if(NXMAX==0)
                cout<<"----- ERREUR : NXMAX absent du fichier de configuration" <<fic_cfg_manip<<endl;

            float R_exp=round(NXMAX*n0/NA);//calcul du rayon à partir de la corde NXMAX

            float rayon=R_exp;

            float tailleTheoPixelUborn=tailleTheoPixelHolo*(float)dimROI.x/(2*NXMAX),
            delta_zmax=rayon*cos(theta);//float?
            cout<<"##################### INFO RECONSTRUCTION ##################"<<endl;
            //cout<<"test: "<<(float)dimROI.x/(2*NXMAX)<<endl;
            cout<<"Taille en pixel de UBorn="<<2*NXMAX<<endl;
            cout<<"Rayon EWALD utilisé="<<rayon<<endl;
            cout<<"Rayon theorique="<<Rth<<" pixels"<<endl;
            printf("theta : %f rad\n: %f, rayon théorique Ewald : %f \n",theta,Rth);

            cout<<"Rayon d'EWALD calculé à partir du fichier de config="<<rayon<<endl;
            cout<<"taille théorique Uborn en nm, avant tomo : "<<tailleTheoPixelUborn<<"nm"<<endl;
            cout<<"##########################################################"<<endl;
          //  const int dim_final=4*NXMAX;//512;//peut être différent de 4*NXMAX, mais l'image final sera (dé)zoomée;
           // float K_tomo_final=(float)dim_final/(4*NXMAX)/2;//facteur de zoom final/tomo/facteur 2 pour synthèse ouverture.

          //  float tailleTheoPixelTomo=tailleTheoPixelUborn*K_tomo_final;
          //  double zoom=double(dim_final)/double(4*n0*TpCam*dimROI.x/(Gt*lambda)*NA/n0);

         /*   if(zoom!=1) {
                    cout<<"dim_final forcée à "<<dim_final<<" au lieu de "<<4*NXMAX<<"-->facteur de zoom="<<zoom<<endl;
                    cout<<"nouvelle taille pixel tomo imposée par le zoom numérique= "<<tailleTheoPixelTomo/zoom<<" nm"<<endl;
                    cout<<endl<<"###########################################"<<endl;
            }*/

            ///--------------VARIABLE DE TAILLE : holo, tomo, nbpixel, decalage
            Var2D dimSpctHolo={2*NXMAX, 2*NYMAX}, dimSpctTomo={4*NXMAX,4*NYMAX}, NMAX={NXMAX,NYMAX}, DecalU_Born=NMAX,fftDecalTomo=dimSpctHolo;
            const int NbPixHolo=dimSpctHolo.x*dimSpctHolo.y;
            // cout <<"Dim finale : " << round(4*zoom*n0*TpCam*dimROI.x/(G*lambda)*NA/n0) << " pixels cube\n";



            //NXMAX_Rf=round(zoom*n0*Tps*DIMX_CCD2/(G*lambda)*NA/n0); //dimension espace final
            //cout<<"NXMAX_Rf="<<NXMAX_Rf<<endl;

            ///------ Définition de la fenêtre de sélection du Hors-axe centrée sur (cx, cy)
            size_t upperleft_x = circle_cx - (2*NXMAX / 2);
            size_t upperleft_y = circle_cy - (2*NYMAX / 2);
            assert (upperleft_x < WINDOW_X);
            assert (upperleft_y < WINDOW_Y);
            // création d'une coupe 2D centrée sur le cercle support du champ hors axe (=spectre du champ complexe, Uborr ou UD)
            nbCplx *TF_UBorn=NULL;
            int nb_proj=0;

            fflush(stdout);
            cout.flush();

            /// //////////////////////////////////////////////////fin calcul constante//////////////////////////

            ///------------Variable 2D avant découpe, taille=ROI découpée sur la caméra : HOLOGRAMME------------------------------------------
            unsigned char* holo1=new unsigned char[NbPixROI2d];
            nbCplx *holo=new nbCplx[NbPixROI2d];

            // 	unsigned char* noir_camera=new unsigned char[N];
            double* masque=new double[NbPixROI2d]; //masque hamming

            unsigned char* cache_jumeau=new unsigned char[100];//cache objet jumeau

            nbCplx *TF_Holo=new nbCplx[NbPixROI2d];  //Avant Crop à NXMAX (hologramme);
            nbCplx *holo_shift=new nbCplx[NbPixROI2d];
            nbCplx *TF_Holo_centre=new nbCplx[NbPixROI2d];


        ///---------------Variables 2D Après Découpe  à 2NXMAX : variables liées à U_56-----------------------------------------------------
           // nbCplx *fft=new nbCplx[4*NXMAX*NYMAX];
            Var2D dim2DHA= {2*NXMAX,2*NYMAX},decal2D= {NXMAX,NYMAX}, coinHA={circle_cx-NXMAX,circle_cy-NYMAX};///coordonnée du coin haut gauche de découpe du hors axe.
            const unsigned int NbPixU_Born=dim2DHA.x*dim2DHA.y;

            double* ampli_ref=new double[NbPixU_Born]; //reservation reference pour "tentative" de correction de la ref
            nbCplx *TF_UBornTot=new nbCplx[NbPixU_Born*NbAngle];///variable stockant les N champs complexes decoupés depuis la zone 1024 (pour utilsier wisdom en 1024)
            nbCplx *TF_UBorn_norm=new nbCplx[NbPixU_Born], *TF_UBorn_I=new nbCplx[4*NXMAX*NXMAX];
            nbCplx *TF_UBorn_normC=new nbCplx[NbPixU_Born], *TF_UBorn_normI=new nbCplx[NbPixU_Born];
            nbCplx *UBorn_I=new nbCplx[NbPixU_Born];
            nbCplx *UBorn=new nbCplx[NbPixU_Born];
            double *UBornAmp=new double[NbPixU_Born];
            double* tukeyHA=new double[NbPixU_Born]; //masque hamming
            nbCplx *UBornFinal=new nbCplx[4*NXMAX*NXMAX];
            nbCplx *UBornFinalDecal=new nbCplx[4*NXMAX*NXMAX];
            double *phase_stack=new double[NbPixU_Born*NbAngle];///pile 3D pour les 600 holgrammes de taille NbPIxUborn
            double *Amp_stack=new double[NbPixU_Born*NbAngle];///pile 3D pour les 600 holgrammes de taille NbPIxUborn
            double *phase2Pi=new double[4*NXMAX*NXMAX];
            double *tab_posSpec=new double[NbAngle*2];  ///stockage des speculaires pour exportations
            double *UnwrappedPhase=new double[4*NXMAX*NXMAX];
            double *PhaseFinal=new double[4*NXMAX*NXMAX];
            double *UBornAmpFinal=new double[NbPixU_Born];
            double *poly_aber=new double[4*NXMAX*NXMAX];
            double *TF_champMod=new double[NbPixU_Born];
            double *centre=new double[NbPixU_Born];//pour mettre la position des centres translatés, on crée une variable 2D de la taille d'un plan apres tomo
            double *mediane_phase=new double[NbPixU_Born];
            double *mediane_Amp=new double[NbPixU_Born];
            double *col_phase=new double[NbAngle];///vecteur colonne pour calcul mediane
            double *col_Amp=new double[NbAngle];//idem

         //pour mettre la position des centres translatés, on crée une variable 2D de la taille d'un plan apres tomo

        ///-----Réservation chrono----------------
            clock_t t_init,t_fin,t_init_deroul,t_fin_deroul, t_init_aber,t_fin_aber,t_init_TF,t_fin_TF,t_init_decal,t_fin_decal;
            clock_t t_init_2pi, t_fin_2pi, t_ini_born,t_fin_born;
            double t_total,t_total_deroul,t_total_aber,t_total_TF,t_total_decal,t_total_2pi,t_total_born;

        ///Chargement de la réference en module, decoupe ROI, mise à l'echelle à 2*NXMAX et calcul amplitude
             //charge_ref(ampli_ref, CheminModRef, NbPixU_Born, dimROI, coinROI,dim2DHA);
            charge_refOpenCV(ampli_ref, CheminModRef,NbPixU_Born,  dimROI, coinROI, dim2DHA);

//            SAV2(ampli_ref,NbPixU_Born,output_dir+"Tomo_Images/ampli_ref.bin",FLOAT,"wb");
            ///Charger masque aberration
             Mat mask = imread(Chemin_mask, 0);
              if(! mask.data ){
                cout <<  "Masque non trouvé, création masque unité" << std::endl ;
                 mask=255*Mat::ones(dim2DHA.x,dim2DHA.y, CV_8UC1);
                }
                else{
                cout<<"chargement du masque"<<Chemin_mask<<endl;
                }
             if(mask.rows!=2*NXMAX){
                cout<<"Problème aberrations : masque "<<Chemin_mask<< " de largeur "<<mask.rows<<", mais image de largeur "<<2*NXMAX<<endl;             }


            //SAV2(ampli_ref,NbPixU_Born,chemin_result+"ampli_ref",FLOAT,"wb");
            //Mat src0=Mat(dim2DHA.x, dim2DHA.y,CV_64F, ampli_ref);

            ///Chargement du masque pour correction aberration/ampli

            mask.convertTo(mask, CV_8U);

            bool b_CorrAber=false;
            bool b_Deroul=false;
            bool b_Born=true;

            b_CorrAber=extract_val("C_ABER",fic_cfg_recon);///corriger les aberrations?
            b_Deroul=extract_val("DEROUL",fic_cfg_recon);///Dérouler la phase?
            b_Born=extract_val("BORN",fic_cfg_recon);///Born vrai ? Sinon Rytov
            if(b_Born==true){
            cout<<"BORN=1"<<endl;
            }
            else{
                cout<<"RYTOV=1"<<endl;
            }


            int NbPtOk=countM(mask);
            cout<<"NbPtOk="<<NbPtOk<<endl;
            if(NbPtOk==0)
               b_CorrAber=false;

            ///corriger les aberrations sur l'amplitude de la ref--------------------------
            Mat src=Mat(dim2DHA.x, dim2DHA.y, CV_64F, ampli_ref);
            Mat AmpRef_corr(dim2DHA.x, dim2DHA.y, CV_64F);

            if(b_CorrAber==true){
                    cout<<"Corr_Aber=1"<<endl;
                AmpRef_corr=ampliCorr(src, mask,poly_aber,3,NbPtOk);

                for(size_t x=0;x<dim2DHA.x;x++){
                    for(size_t y=0;y<dim2DHA.y;y++){
                        size_t cpt=x+y*dim2DHA.x;
                        ampli_ref[cpt]=AmpRef_corr.at<double>(y,x);
                    }
                }
            }

            //SAV2(ampli_ref,NbPixU_Born,chemin_config_GUI+"truc.txt",FLOAT,"wb");
            ///normaliser l'énergie (L2) dans l'amplitude de la ref--------------------------
            double moy_ref=moyenne(ampli_ref,NbPixU_Born);
            for(size_t cpt=0;cpt<NbPixU_Born;cpt++)
                ampli_ref[cpt]=ampli_ref[cpt]/moy_ref;

            cout<<"valeur moyenne="<<moy_ref<<endl;

            ///Variable de masquage :  fenetre de Tukey,
            float alpha=0.1;//coeff pour le masque de tuckey
            masque=tukey2D(dimROI.x,dimROI.y,0.05);
            tukeyHA=tukey2D(dim2DHA.x,dim2DHA.y,alpha);
            cout<<"Repertoire de sortie : "<<chemin_result<<endl;
            SAV2(masque,dimROI.x*dimROI.y,chemin_result+"/tukey.bin",FLOAT,"wb");

            Var2D posSpec={0,0};
            FILE* test_existence;//nom bidon pour tester l'existence des fichiers

        ///début de la boucle sur les angles/////////////////////////

       /* int bool_wisdom1024=fftw_import_wisdom_from_filename("/home/mat/Projet_tomo/Tomo_Images/wisdom/test1024.wisdom");//charger ou calculer le fichier wisdom
                if(bool_wisdom1024==0){
                        cout<<"################-----Calcul wisdom 1024 (~10 min)"<<endl;
                        prepare_wisdom(dimROI,"/home/mat/Projet_tomo/Tomo_Images/wisdom/test1024.wisdom");
                }*/

        ///##################################################début de la boucle sur les angles/////////////////////////
        char charAngle[4+1];
        t_init=clock();
        size_t  NbAngleOk=0,NumAngle=0;
        for(int cpt_angle=premier_plan; cpt_angle<premier_plan+NbAngle; cpt_angle=cpt_angle+SautAngle) //boucle sur tous les angles : TF1024+decoupe
        {
                if((cpt_angle-100*(cpt_angle/100))==0)
                printf("cpt_angle=%i\n",cpt_angle);
                ///Concaténer le numéro d'angle dans le nom
                sprintf(charAngle,"%03i",cpt_angle);
                //string TamponChemin2=Dossier_acquiz+"i"+charAngle+"-001.pgm";
                string TamponChemin2=Dossier_acquiz+"i"+charAngle+".pgm";
                //cout<<TamponChemin2<<endl;
                test_existence = fopen(TamponChemin2.c_str(), "rb");
                if(test_existence!=NULL) {
                    fclose(test_existence);
                    ///charger l'hologramme
                    charger_image2D(holo1,TamponChemin2, coinROI, dimROI);
                    holo2TF_UBorn(holo1, TF_UBornTot+NbAngleOk*NbPixU_Born,dimROI,decalROI, dim2DHA, coinHA,NbAngleOk,masque);
                    NbAngleOk++;
                    //cout<<NbAngleok<<endl;
                }
                else {
                        printf("fichier %i inexistant\n",cpt_angle);
                        //fclose(test_existence);
                }
        }
        //SAV_Re(TF_UBornTot,NbPixU_Born*NbAngle,"/home/mat/Projet_tomo/Tomo_Images/TF_Uborn3d.bin",FLOAT,"wb");
        //recharger wisdom 2D mini
       //NumAngle=0;
       fftw_forget_wisdom();

     /*  int bool_wisdom2D=fftw_import_wisdom_from_filename("/home/mat/Projet_tomo/Tomo_config/wisdom/test2D.wisdom");//charger ou calculer le fichier wisdom
                if(bool_wisdom2D==0)
                    prepare_wisdom(dimSpctHolo,"/home/mat/Projet_tomo/Tomo_config/wisdom/test2D.wisdom");*/

        for(int cpt_angle=premier_plan; cpt_angle<premier_plan+NbAngleOk; cpt_angle=cpt_angle+SautAngle) //boucle sur tous les angles : correction aberrations
        {
            ///Recherche de la valeur maximum du module dans ref non centré-----------------------------------------
            TF_UBorn=TF_UBornTot+NumAngle*NbPixU_Born;//récupérer la Tf2 dans le volume des N hologrammes
            // SAV_Re(TF_UBorn, NbPixU_Born, "/home/mat/Projet_tomo/Tomo_Images/TF2d/TF_UBorn", FLOAT,"a+b");
            SAV_Re2(TF_UBorn, NbPixU_Born, chemin_result+"/TF2d/TF_UBorn", FLOAT,"a+b");

            NumAngle++;
            int cpt_max=coordSpec(TF_UBorn, TF_champMod,NMAX);
            double  max_part_reel = TF_UBorn[cpt_max].Re,///sauvegarde de la valeur cplx des  spéculaires
                    max_part_imag = TF_UBorn[cpt_max].Im,
                    max_module = TF_champMod[cpt_max];

            int kxmi=cpt_max%(2*NXMAX), kymi=cpt_max/(2*NYMAX);
            posSpec={kxmi,kymi};///coord informatique speculaire
            tab_posSpec[cpt_angle]=(double)posSpec.x;
            tab_posSpec[cpt_angle+NbAngle]=(double)posSpec.y;
            centre[kxmi*2*NMAX.x+kymi]=cpt_angle;

            if(b_CorrAber==true){
                calc_Uborn(TF_UBorn,UBorn,dim2DHA,posSpec);///--/!\ recale le spectre dans Uborn!
                  //SAV_Im(UBorn,NbPixU_Born,"/home/mat/Projet_tomo/Tomo_Images/Uborn_im_norm_avant.raw",FLOAT,"a+b");
                //t_init_2pi=clock();
                calcPhase2pi(UBorn, dim2DHA,phase2Pi);//atan
                //phase2pi(UBorn, dim2DHA,phase2Pi);//asin

                //t_fin_2pi=clock();
                //t_total_2pi=(double)(t_fin_2pi-t_init_2pi)/CLOCKS_PER_SEC+t_total_2pi;
                SAV2(phase2Pi,NbPixU_Born,chemin_result+"/phase.raw",FLOAT,"a+b");
               // t_init_deroul=clock();
               if(b_Deroul==true){
                phaseUnwrapping_Mat(dim2DHA, phase2Pi, UnwrappedPhase);
               }
               else{
                   for(int cpt=0;cpt<NbPixU_Born;cpt++)
                        UnwrappedPhase[cpt]=phase2Pi[cpt];
               }
                //cout<<(double)t_init_deroul/CLOCKS_PER_SEC<<endl;
                //t_fin_deroul=clock();
               // t_total_deroul=(double)(t_fin_deroul-t_init_deroul)/CLOCKS_PER_SEC+t_total_deroul;
               // SAV(UnwrappedPhase,NbPixU_Born,"/home/mat/Projet_tomo/Tomo_Images/phase_deroul_avec_norm_spec.raw",FLOAT,"a+b");
                SAV2(UnwrappedPhase,NbPixU_Born,chemin_result+"/phase_deroul_avec_norm_spec.raw",FLOAT,"a+b");
                ///Correction aberration phase-------------------------------
                src=Mat(dim2DHA.x, dim2DHA.y,CV_64F, UnwrappedPhase);
                Mat Phase_corr(dim2DHA.x, dim2DHA.y, CV_64F);
                t_init_aber=clock();
                Phase_corr=aberCorr(src, mask,poly_aber,3,  NbPtOk);
                t_fin_aber=clock();
               t_total_aber=(double)(t_fin_aber-t_init_aber)/CLOCKS_PER_SEC+t_total_aber;
               // SAV(poly_aber,NbPixU_Born,"/home/mat/Projet_tomo/Tomo_Images/poly_aber_phase.raw",FLOAT,"a+b");
                SAV2(poly_aber,NbPixU_Born,chemin_result+"/poly_aber_phase.raw",FLOAT,"a+b");
                // (double*) Phasefinal = (double*)Phase_corr.data;//transtypage sur le pointeur par defaut uchar, mais contient des données en double)//SAV(poly_aber,NbPixU_Born,"/home/mat/Projet_tomo/Tomo_Images/poly_aber_phase.raw",FLOAT,"a+b");

                ///Correction amplitude----------------------------------------
                for(int cpt=0; cpt<(4*NXMAX*NYMAX); cpt++){
                    UBornAmp[cpt]=sqrt(pow(UBorn[cpt].Re,2)+pow(UBorn[cpt].Im,2));
                }
               // SAV(UBornAmp,NbPixU_Born,"/home/mat/Projet_tomo/Tomo_Images/UbornAmp.bin",FLOAT,"a+b");
                SAV2(UBornAmp,NbPixU_Born,chemin_result+"/UbornAmp.bin",FLOAT,"a+b");
                Mat src=Mat(dim2DHA.x, dim2DHA.y, CV_64F, UBornAmp);
                Mat UBornAmp_corr(dim2DHA.x, dim2DHA.y, CV_64F);
                UBornAmp_corr=ampliCorr(src, mask,poly_aber,3,  NbPtOk);
                ///Fin Correction amplitude----------------------------------------

                //reconstruire l'onde complexe/Recalculate the complex field
                for(size_t x=0;x<dim2DHA.x;x++){
                    for(size_t y=0;y<dim2DHA.y;y++){
                        size_t cpt=x+y*dim2DHA.x;
                        size_t cpt3D=x+y*dim2DHA.x+NbPixU_Born*(cpt_angle-1);
                        PhaseFinal[cpt]=Phase_corr.at<double>(y,x);//copie opencV->Tableau
                        //UBornAmp[cpt]=UBornAmp_corr.at<double>(y,x);

                       /* if (mask.at<uchar>(y,x) < 10)///remplacer les valeur aberrantes (diaphragme) par les valeur idéales.
                        {
                            PhaseFinal[cpt]=0;
                            UBornAmp[cpt]=1;
                        }*/
                        //phase_stack[cpt3D]=PhaseFinal[cpt];
                        //Amp_stack[cpt3D]=UBornAmp[cpt];
                        UBornAmpFinal[cpt]=UBornAmp[cpt];///ampli_ref[cpt];//-1;//division par ampli ref, puis soustraction Uinc (=1 car normalisé), déplacé ci dessous

                        if(b_Born==true){//BORN
                            UBornFinal[cpt].Re= (UBornAmpFinal[cpt]-1)*cos(PhaseFinal[cpt])*tukeyHA[cpt];
                            UBornFinal[cpt].Im=(UBornAmpFinal[cpt]-1)*sin(PhaseFinal[cpt])*tukeyHA[cpt];
                        }
                        else{//RYTOV
                            UBornFinal[cpt].Re= log(sqrt(UBornAmpFinal[cpt]*UBornAmpFinal[cpt]))*tukeyHA[cpt];
                            UBornFinal[cpt].Im=PhaseFinal[cpt]*tukeyHA[cpt];
                        }

                    }
                }
                //SAV(UBornAmp,NbPixU_Born,"/home/mat/Projet_tomo/Tomo_Images/UBornAmp_corr.raw",DOUBLE,"a+b");
               // SAV(UBornAmpFinal,NbPixU_Born,"/home/mat/Projet_tomo/Tomo_Images/UBornAmpFinal.raw",FLOAT,"a+b");
               // SAV(PhaseFinal,NbPixU_Born,"/home/mat/Projet_tomo/Tomo_Images/phase_final_apres_corr_et_norm_spec.raw",FLOAT,"a+b");
                SAV2(UBornAmpFinal,NbPixU_Born,chemin_result+"/UBornAmpFinal.raw",FLOAT,"a+b");
                SAV2(PhaseFinal,NbPixU_Born,chemin_result+"/phase_final_apres_corr_et_norm_spec.raw",FLOAT,"a+b");

                ///Recalculer la TF décalée pour le programme principal.
                Var2D recal={kxmi,kymi};
                decal2DCplxGen(UBornFinal,UBornFinalDecal, dim2DHA,NMAX);
                t_init_TF=clock();
                TF2Dcplx(UBornFinalDecal,TF_UBorn_norm,dim2DHA);
              //  SAV_Re(TF_UBorn_norm,NbPixU_Born, "/home/mat/Projet_tomo/Tomo_Images/TF2d/Tf_UBornfinalDecal_Re.raw", DOUBLE, "a+b");
                //t_fin_TF=clock();
                //t_total_TF=(double)(t_fin_TF-t_init_TF)/CLOCKS_PER_SEC+t_total_TF;

                ///redécalage de la TF pour projection.
                //decal2DCplxGen(TF_UBorn_norm,TF_UBorn_normC,dim2DHA,recal);
               // SAV_Re(UBornFinal,NbPixU_Born, "/home/mat/Projet_tomo/Tomo_Images/UBornfinal_Re.raw", DOUBLE, "a+b");
                SAV_Re2(UBornFinal,NbPixU_Born, chemin_result+"/UBornfinal_Re.raw", DOUBLE, "a+b");
               // SAV_Im(UBornFinal,NbPixU_Born, "/home/mat/Projet_tomo/Tomo_Images/UBornfinal_Im.raw", DOUBLE, "a+b");
                SAV_Im2(UBornFinal,NbPixU_Born,chemin_result+"/UBornfinal_Im.raw", DOUBLE, "a+b");

            }
            else{//sauvegarde onde avec aberration
                //cout<<"normalise spec"<<endl;
                for(int cpt=0; cpt<(4*NXMAX*NYMAX); cpt++)
                    {
                     TF_UBorn_norm[cpt].Re=(TF_UBorn[cpt].Re*max_part_reel+TF_UBorn[cpt].Im*max_part_imag)/max_module;
                     TF_UBorn_norm[cpt].Im=(TF_UBorn[cpt].Im*max_part_reel-TF_UBorn[cpt].Re*max_part_imag)/max_module;
                    }
               // SAV_Re(TF_UBorn,NbPixU_Born, "/home/mat/Projet_tomo/Tomo_Images/TF_UBorn_norm_Re.raw", DOUBLE, "a+b");
                calc_Uborn(TF_UBorn_norm,UBorn,dim2DHA,posSpec);///--/!\ recale le spectre dans support Uborn!

                //SAV_Re(UBorn,NbPixU_Born, "/home/mat/Projet_tomo/Tomo_Images/UBornfinal_Re.raw", DOUBLE, "a+b");
                SAV_Re2(UBorn,NbPixU_Born, chemin_result+"/UBornfinal_Re.raw", DOUBLE, "a+b");
               // SAV_Im(UBorn,NbPixU_Born, "/home/mat/Projet_tomo/Tomo_Images/UBornfinal_Im.raw", DOUBLE, "a+b");
                SAV_Im2(UBorn,NbPixU_Born, chemin_result+"/UBornfinal_Im.raw", DOUBLE, "a+b");
            }

        }//fin de boucle for sur tous les angles

        cout<<"déroulement, temps par angle="<<t_total_deroul/NbAngle<<endl;
        cout<<"calcul 2pi, temps par angle="<<t_total_2pi/NbAngle<<endl;
        cout<<"correction aberration phase, temps par angle="<<t_total_aber/NbAngle<<endl;
        cout<<"une TF, temps par angle="<<t_total_TF/NbAngle<<endl;
        cout<<"un decalage decal2DCplxGen, temps par angle="<<t_total_decal/NbAngle<<endl;
        ///SAUVER LES PARAMETRES UTILES À RECONSTRUCTION
        cout<<"NXMAX="<<NXMAX<<endl;
        double param[]={NXMAX,NbAngle,rayon,dimROI.x,tailleTheoPixelHolo};
        cout<<"param[1]="<<param[0]<<endl;;
        int NbParam=sizeof (param)/sizeof (double);
       // SAV(param,NbParam, "/home/mat/Projet_tomo/Tomo_Images/parametres.raw", DOUBLE, "wb");
       // SAV(tab_posSpec,NbAngle*2,"/home/mat/Projet_tomo/Tomo_Images/tab_posSpec.raw",DOUBLE,"wb");

        SAV2(param,NbParam, chemin_result+"/parametres.raw", DOUBLE, "wb");
        SAV2(tab_posSpec,NbAngle*2,chemin_result+"/tab_posSpec.raw",DOUBLE,"wb");
        // fic_posSpec.close();

        delete[] phase_stack; delete[] mediane_phase; delete[] mediane_Amp;delete[] col_Amp;delete[] col_phase;
        delete[] tukeyHA;
        delete[] TF_Holo_centre; delete[] TF_Holo;
        delete[] TF_UBorn_norm;
        delete[] TF_UBornTot;
        delete[] TF_UBorn_I;delete[] TF_UBorn_normC;delete[] TF_UBorn_normI;
        delete[] UBorn_I, delete[] UBorn,delete[] UBornFinal,delete[] UBornFinalDecal;delete[] UBornAmp,delete[] UBornAmpFinal;
        delete[] phase2Pi,delete[] UnwrappedPhase, delete[] PhaseFinal;
        delete[] TF_champMod;
        delete[] holo_shift;
        delete[] masque;
        delete[] cache_jumeau;
        delete[] holo1, delete[] holo;
        delete[] poly_aber;delete[] tab_posSpec;


        t_fin=clock();
        t_total=(double)(t_fin-t_init)/CLOCKS_PER_SEC;
        //printf("points_faux!! : %i,points_faux\n",points_faux);
        ///-------------exportation des centres---------------
       // SAV(centre, 4*NXMAX*NYMAX, "/home/mat/Projet_tomo/Tomo_Images/centres.bin", INT,"wb");
        SAV_Tiff2D(centre, "/home/mat/Projet_tomo/Tomo_Images/centres.tif", 2*NXMAX );
        delete[] centre;

        void fftw_cleanup_threads(void);//libération memoire allouée pour les threads
        // delete[] valeur_module_shift;

        return 0;

}

