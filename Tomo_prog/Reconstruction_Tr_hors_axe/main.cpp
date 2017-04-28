#include "projet.h"
#include "fonctions.h"
#include "deroulement.h"
#include "Correction_aberration.h"
//#include "CImg.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <octave-2.9.9/octave/oct.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <Magick++.h>
#include <fftw3.h>

#include <cstring>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <omp.h>

//using namespace cimg_library;
using namespace std;
using namespace std::chrono;
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
  /*     usage(argc, argv);

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
        }*/
        string home=getenv("HOME");
        string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
        string chemin_config_GUI=getenv("HOME")+fin_chemin_gui_tomo;
        //string repertoire_config=extract_string("CHEMIN_CONFIG",home+fin_chemin_gui_tomo);
        string chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);


        cout<<"##################"<<endl;
        cout<<"# RECONSTRUCTION #"<<endl;
         cout<<"##################"<<endl;
        ///-----Recharger les paramètres--------///
        const int nbParam=5;
        double parametre[nbParam];
        lire_bin(chemin_result+"parametres.raw",parametre,64,nbParam);

        const int NXMAX=parametre[0],NYMAX=parametre[0],NbAngle=parametre[1],rayon=parametre[2],tailleTheoPixelHolo=parametre[4]
       ;
        cout<<"NXMAX="<<parametre[0]<<" | NbAngle="<<NbAngle<<"| rayon="<<rayon<<" | TpHolo="<<tailleTheoPixelHolo<<endl;
        Var2D NMAX={NXMAX,NXMAX},dimROI={2*NXMAX,2*NXMAX};
        ///------Récupérer spéculaire---
        double *TabPosSpec=new double[2*NbAngle];
        cout<<"lecture posSpec"<<endl;

        lire_bin(chemin_result+"tab_posSpec.raw",TabPosSpec,64,2*NbAngle);


        float precision;//pour conversion de type en flottant (32 bits) lors de l'écriture finale

        ///--------------constantes multithread FFTW
        int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        int nthreads=3;
        printf("fftwThreadInit: %i\n",fftwThreadInit);

        ///---Variable de dimensions (taille pixel, volume, zoom)-------------------
        float tailleTheoPixelUborn=tailleTheoPixelHolo*(float)(1024/(2*NXMAX));
        cout<<"test: "<<(float)dimROI.x/(2*NXMAX)<<endl;
        cout<<"NMAX="<<NXMAX<<endl<<"Taille en pixel de UBorn="<<2*NXMAX<<endl;
        cout<<"taille théorique Uborn en nm, avant tomo : "<<tailleTheoPixelUborn<<"nm"<<endl;

        const  int dim_final=512;//peut être différent de 4*NXMAX, mais l'image final sera (dé)zoomée;
        cout<<"Dimension forcée à "<<dim_final<<endl;
        float K_tomo_final=(float)dim_final/(4*NXMAX)/2;//facteur de zoom final/tomo/facteur 2 pour synthèse ouverture.

        float tailleTheoPixelTomo=tailleTheoPixelUborn/2;
        cout<<"taille theorique pixel tomo="<<tailleTheoPixelTomo<<"nm"<<endl;
        double zoom=double(dim_final)/double(4*NXMAX);//double(dim_final)/double(4*n0*TpCam*dimROI.x/(Gt*lambda)*NA/n0);//dim_final/4NXMAX
        if(zoom!=1) {
                    cout<<"dim_final forcée à "<<dim_final<<" au lieu de "<<4*NXMAX<<"-->facteur de zoom numérique="<<zoom<<endl;
                    //cout<<"nouvelle taille pixel tomo imposée par le zoom numérique= "<<tailleTheoPixelTomo/zoom<<" nm"<<endl;
                    printf("nouvelle taille pixel tomo imposée par le zoom numérique= %f",tailleTheoPixelTomo/zoom);
                    cout<<endl<<"###########################################"<<endl;
        }



        ///---------------Calcul de quelques constantes, afin d'éviter de les calculer dans la boucle--------------
        Var3D  decal3DTF={(int)round(dim_final/2),(int)round(dim_final/2),(int)round(dim_final/2)},
        dimVol={dim_final,dim_final,dim_final};
        const int ///demi longueur pour recalculer les indices
        dimPlanFinal=round(dim_final*dim_final),
        N3D=dimVol.x*dimVol.x*dimVol.x;///nb pixel dans l'espace final tomographique 3D
        printf("dimx espace 3D, dimVolX: %i\n",dimVol.x);
        fflush(stdout);
        cout.flush();

        ///---------------Recharger le champ complexe depuis sav prétraitement-----------------------------------------------------
        // nbCplx *fft=new nbCplx[4*NXMAX*NYMAX];
        Var2D dim2DHA= {2*NXMAX,2*NYMAX},decal2D= {NXMAX,NYMAX};
        const int NbPixU_Born=dim2DHA.x*dim2DHA.y;
        double *UBornFinal3D_Re=new double[NbPixU_Born*NbAngle];
        double *UBornFinal3D_Im=new double[NbPixU_Born*NbAngle];
        double *mask_tukey2D=new double[NbPixU_Born];
        nbCplx *UBornFinal3D=new nbCplx[NbPixU_Born*NbAngle];
        cout<<"lecture uborn Re et Im"<<endl;

        lire_bin(chemin_result+"UBornfinal_Im.raw",UBornFinal3D_Im,64,NbPixU_Born*NbAngle);

        lire_bin(chemin_result+"UBornfinal_Re.raw",UBornFinal3D_Re,64,NbPixU_Born*NbAngle);

        for(int cpt=0;cpt<NbPixU_Born*NbAngle;cpt++){
            UBornFinal3D[cpt].Re=UBornFinal3D_Re[cpt];
            UBornFinal3D[cpt].Im=UBornFinal3D_Im[cpt];
        }

       // SAV_Re(UBornFinal3D,NbPixU_Born*NbAngle,"/home/aziz/Projet_tomo/Tomo_Images/test_recharge_UBornfinal_3D.raw",FLOAT,"wb");
        delete[] UBornFinal3D_Re;delete[] UBornFinal3D_Im;

        nbCplx *UBornFinal2D=new nbCplx[NbPixU_Born];//CHAMP2D
        nbCplx *UBorn2DFinalDecal=new nbCplx[NbPixU_Born];//CHAMP2D décalé pour TF
        nbCplx *TF_UBorn_normC=new nbCplx[NbPixU_Born], *TF_UBorn_normI=new nbCplx[NbPixU_Born];//FFT2D

        ///---------espace3D final --------------------------------------------------
        const int N_tab=dim_final*dim_final*dim_final;//Nombre de pixel dans l'esapce final (indice max du tableau 3D+1)
        printf("N_tab (indice espace 3D)= %i \n",N_tab);
        nbCplx *TF3D_PotObj=new nbCplx[N_tab];///FFT3D
        double *sup_redon=new double[N_tab];//pour moyennage fréquentiel

        ///-----Réservation chrono----------------

        clock_t
        temps_depart  = clock(), /*temps de départ de programme */
        temps_initial = clock (), /* temps initial en micro-secondes */
        temps_final=0,/* temps d'arrivée du programme */
        temps_arrivee;   /* temps final en micro-secondes */
        float temps_cpu=0,     /* temps total en secondes */
        temps_total=0;
        temps t1;


        ///Variable de masquage :  fenetre de Tukey,
            float alpha=0.1;//coeff pour le masque de tuckey
            mask_tukey2D=tukey2D(dim2DHA.x,dim2DHA.y,alpha);

        printf("*******************************************\n");
        printf("remplissage de l'espace réciproque\n");
        printf("*******************************************\n");
        printf("  \\,`//\n");
        printf(" _).. `_\n");
        printf("( __  -\\ \n");
        printf("    '`.\n");
        printf("   ( \\>\n");
        printf("   _||_ \n");
        const int premier_plan=0;
        Var2D posSpec={0,0};
        Var2D recal={0,0};
        string tmp=chemin_result+"/wisdom/test2D.wisdom";
        int bool_wisdom2D=fftw_import_wisdom_from_filename(tmp.c_str());//charger ou calculer le fichier wisdom
            if(bool_wisdom2D==0){
                prepare_wisdom2D(dimROI,tmp.c_str());
            }

        for(int cpt_angle=premier_plan; cpt_angle<NbAngle; cpt_angle++) //boucle sur tous les angles
        {
            //récupérer spéculaire puis champ cplx depuis sauvegarde prétraitement
            posSpec={(int)TabPosSpec[cpt_angle],(int)TabPosSpec[NbAngle+cpt_angle]};

                //cout<<"posSpec.x="<<posSpec.x<<" | posSpec.y="<<posSpec.y<<endl;
            for(int cpt=0;cpt<NbPixU_Born;cpt++){
                UBornFinal2D[cpt].Re=UBornFinal3D[cpt+cpt_angle*NbPixU_Born].Re*mask_tukey2D[cpt];
                UBornFinal2D[cpt].Im=UBornFinal3D[cpt+cpt_angle*NbPixU_Born].Im*mask_tukey2D[cpt];
            }

            decal2DCplxGen(UBornFinal2D,UBorn2DFinalDecal, dim2DHA,NMAX);
            //SAV_Im2(UBornFinal2D,NbPixU_Born,chemin_result+"Uborn_test.raw",FLOAT,"a+b");
            //SAV_Im2(UBorn2DFinalDecal,NbPixU_Born,chemin_result+"UBornFinalDecal.raw",FLOAT,"a+b");///sans doute supprimable si on supprime aussi le dernier circshift3D
            TF2Dcplx(UBorn2DFinalDecal,TF_UBorn_normI,dim2DHA);

            recal={posSpec.x,posSpec.y};
            ///recaler le spectre à la position du spéculaire (la variable est attendue ainsi par retroprpag)
            decal2DCplxGen(TF_UBorn_normI,TF_UBorn_normC,dim2DHA,recal);
            SAV_Im2(TF_UBorn_normC,NbPixU_Born,chemin_result+"/TF2d/TF_UBorn_normC.raw",FLOAT,"a+b");

            if((cpt_angle-100*(cpt_angle/100))==0)
                printf("cpt_angle=%i\n",cpt_angle);
            temps_initial = clock ();

            //cout<<"POsSPec.x="<<posSpec.x<<"|PosSpec.y="<<posSpec.y<<endl;
            //cout<<"TabPosSpec[cpt_angle]="<<TabPosSpec[cpt_angle]<<"|TabPosSpec[cpt_angle+NbAngle]"<<TabPosSpec[cpt_angle]<<endl;
            retroPropag_Born(TF3D_PotObj, TF_UBorn_normC, sup_redon, dim_final, posSpec, decal3DTF, NMAX, rayon);  ///--Mapping 3D=retropropagation
        }//fin de boucle for sur tous les angles


        delete[] TF_UBorn_normC;
        delete[] TabPosSpec;
        printf("temps_total proj : %f \n",temps_total);
        SAV2(sup_redon, N_tab, chemin_result+"/sup_redon_C.bin", FLOAT,"wb");

        ///---------------------mesure du temps de calcul--------------------------

        temps_final = clock ();
        temps_cpu = (temps_final - temps_initial) * 1e-6/NbAngle;
        printf("temps moyen pour 1 angle: %f\n",temps_cpu);
        temps_cpu = (temps_final - temps_initial) * 1e-6;
        printf("temps total pour %i angle(s): %f\n", NbAngle, temps_cpu);
        temps_initial=clock();//enclenchement du chronometre

        ///------------------moyennage par sup_redon------------------------------------------------
        for(int cpt=0; cpt<N_tab; cpt++) {
                if (sup_redon[cpt]==0) { ////////////////remplace les 0 de sup_redon par des 1---> evite la division par 0
                        //printf("sup_redon= %f \n" , sup_redon[cpt]);
                        sup_redon[cpt]=1;
                }
                TF3D_PotObj[cpt].Re = TF3D_PotObj[cpt].Re/sup_redon[cpt];//moyennage par sup_redon
                TF3D_PotObj[cpt].Im = TF3D_PotObj[cpt].Im/sup_redon[cpt];
        }

        temps_final = clock ();
        temps_cpu = (temps_final - temps_initial) * 1e-6;
        printf("temps apres normalisation : %lf\n",temps_cpu);
        delete[] sup_redon;
        //SAV(TF3D_PotObj_Re, N_tab, "/home/aziz/Projet_tomo/Tomo_Images/TF2d_apres_masquage/TF3D_PotObj_Re_norm_avant.bin", float,"wb");
        //interp3D(TF3D_PotObj_Re,  dimVolX,dimVolX, dimVolX);  //////////////interpolation dans Fourier
        //interp3D(TF3D_PotObj_Im,  dimVolX,dimVolX, dimVolX);

        //////////////////////////////papillon binarisé
        double *papillon_masque=new double[N_tab];
        for(int compteur=0; compteur<N_tab; compteur++) {
                if(TF3D_PotObj[compteur].Re==0)
                        papillon_masque[compteur]=0;
                else
                        papillon_masque[compteur]=1;
        }

        //SAV(papillon_masque, N_tab, "/home/aziz/Projet_tomo/Tomo_Images/TF2d_apres_masquage/papillon_masque.bin", FLOAT,"wb");
        delete[] papillon_masque;
        // SAV(TF3D_PotObj_Re, N_tab, "/home/aziz/Projet_tomo/Tomo_Images/TF2d_apres_masquage//home/mat/TF3D_PotObj_Re_norm.bin", FLOAT,"a+b");*/

        ///--------------------------------circshift avant TF3D
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        printf("*******************************************\n");
        printf("circshift avant TF3D \n");
        printf("*******************************************\n");

        //////////// tableaux qui recoivent le circshift
        nbCplx *TF3D_PotObj_shift=new nbCplx[N_tab];

        temps_initial = clock();//enclenchement chronometre
        Var3D    dimFinal= {dimVol.x,dimVol.y,dimVol.z};
        circshift3DCplx(TF3D_PotObj, TF3D_PotObj_shift, dimFinal, decal3DTF);
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

        delete[] TF3D_PotObj;
      /*  fftw_forget_wisdom();
            int bool_wisdom3D=fftw_import_wisdom_from_filename("/home/aziz/Projet_tomo/Tomo_Images/wisdom/test3D.wisdom");//charger ou calculer le fichier wisdom
                if(bool_wisdom3D==0){
                      cout<<"Calcul wisdom 3D (~8heures)"<<endl;
                      prepare_wisdom3D(dimVol,"/home/aziz/Projet_tomo/Tomo_Images/wisdom/test3D.wisdom");
                }*/

        ///------------------------ TF3D-------------------------------------
        nbCplx *PotObj_shift=new nbCplx[N_tab];

        high_resolution_clock::time_point t1v = high_resolution_clock::now();//temps vrai (pas CPU)
        TF3DCplx_INV(TF3D_PotObj_shift, PotObj_shift, dimFinal);

        high_resolution_clock::time_point t2v = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>( t2v - t1v ).count();

        cout << "temps TF3D="<< duration/1000<<" ms"<<endl;
        delete[] TF3D_PotObj_shift;
        ///--------------------------------circshift apres TF3D
        /////////////////////creation des tableaux qui recoivent le circshift
        nbCplx *PotObj3D=new nbCplx[N_tab];
        circshift3DCplx(PotObj_shift, PotObj3D, dimFinal, decal3DTF);

        delete[] PotObj_shift; //////////////////////////////suppression des tableaux non décalés

        ///--------------------------------chrono et écriture.
        temps_final = clock ();
        temps_cpu = (temps_final - temps_initial) * 1e-6;
        printf("temps apres circshift 3D final et calcul du module: %f\n",temps_cpu);
        temps_initial = clock();//enclenchement chronometre


        nbCplx *indice_cplx=new nbCplx[N_tab];
        float pot2ind=1;//-2*kv*kv/n0;
        float  pot2abs=1;//kv;

        for(size_t cpt=0;cpt<N_tab;cpt++){
            indice_cplx[cpt].Re=pot2ind*PotObj3D[cpt].Re; //réfraction
            indice_cplx[cpt].Im=pot2abs*PotObj3D[cpt].Im;
        }
        //Effacer précédent résultats.
        tmp=chemin_result+"/indice.tif";
        if( remove(tmp.c_str()) != 0 )
        perror( "Ancien fichier indice effacé" );
        tmp=chemin_result+"/absorption.tif";
        if( remove(tmp.c_str()) != 0 )
        perror( "Ancien fichier absorption effacé" );
        printf("ecriture des résultats\n");
        SAV_Tiff3D(indice_cplx,chemin_result+"/indice.tif",chemin_result+"/absorption.tif",dim_final);
        //SAV_Re2(indice_cplx, N_tab, chemin_result+"/indice.bin", FLOAT,"wb");
        //SAV_Im2(indice_cplx, N_tab, chemin_result+"/absorption.bin", FLOAT,"wb");

        delete[] PotObj3D; delete[] indice_cplx; delete[] mask_tukey2D;


        void fftw_cleanup_threads(void);//libération memoire allouée pour les threads
        printf("ecriture de PotObj3D.Re, PotObj3D.Im : OK \n");


        temps_arrivee = clock ();
        temps_cpu = (temps_arrivee-temps_depart )/CLOCKS_PER_SEC;
        printf("temps total: %f\n",temps_cpu);
        return 0;

}

