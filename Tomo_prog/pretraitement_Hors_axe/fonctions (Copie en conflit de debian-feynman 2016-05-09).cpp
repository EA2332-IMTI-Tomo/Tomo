#include "fonctions.h"
#define PI 3.14159

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */
 //extrait la valeur (entier) de la chaine "token" du fichier "chemin_fic"
float extract_val(string token,  string chemin_fic)
{
    ifstream fichier(chemin_fic.c_str(), ios::in);  // on ouvre en lecture
    string ligne,motcle,valeurMot,separ=" ";
    float valeur=0;
    vector<std::string> tokens;

    if(fichier)  // si l'ouverture a fonctionné
    {
        while(!fichier.eof()){
            getline(fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            if(ligne[0]!='#')//si pas ligne de commentaire
            tokens.push_back(ligne);
        }
    }
    else
        cerr << "Impossible d'ouvrir le fichier !" << endl;

    int nb_tok=tokens.size();
    for(int cpt=0;cpt<nb_tok;cpt++){
        ligne=tokens[cpt];
        if(ligne!=""){
            int pos_separ=ligne.find(separ);
            int long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,pos_fin)
            if(motcle==token){
            valeurMot=ligne.substr(pos_separ+long_separ,ligne.size()-(motcle.size()+long_separ));
            cout<<motcle<<"="<<valeurMot<<endl;
            valeur=atof(valeurMot.c_str());
            }
        }
    }
    fichier.close();
    return valeur;
}

 int coordSpec(nbCplx* TF_UBorn, double *TF_champMod,Var2D NMAX)
 {
    int cpt_max=0;
    TF_champMod[0]=pow(TF_UBorn[0].Re,2)+pow(TF_UBorn[0].Im,2);

    for(int cpt=1; cpt<(4*NMAX.x*NMAX.y); cpt++) {
        TF_champMod[cpt]=pow(TF_UBorn[cpt].Re,2)+pow(TF_UBorn[cpt].Im,2);
        if(TF_champMod[cpt]>TF_champMod[cpt_max]) {
        cpt_max=cpt;
        }
    }
   /* double  max_part_reel = TF_UBorn[cpt_max].Re,///sauvegarde de la valeur cplx des  spéculaires
    max_part_imag = TF_UBorn[cpt_max].Im,
    max_module = TF_champMod[cpt_max];*/

    //int kxmi=cpt_max%(2*NMAX.x), kymi=cpt_max/(2*NMAX.y);
    //Var2D posSpec={kxmi,kymi};///coord informatique speculaire
    //recalUBorn={-kxmi,-kymi};
    return cpt_max;
 }


double calc_mediane(double entree[], size_t nb_elem)
{
    int pos_med=round(nb_elem/2);
    //float a[] = {9, 8, 7, 6, 5, 0, 1, 2.5, 3, 4,11,125, -1};
    std::nth_element(entree, entree + pos_med, entree + nb_elem);
    double mediane=entree[pos_med];
   // cout<<"mediane : "<<mediane<<endl;
    return mediane;
}


void holo2TF_UBorn(unsigned char *holo1, nbCplx *TF_UBornTot,Var2D dimROI,Var2D decalROI, Var2D dim2DHA, Var2D coinHA, size_t NumAngle, double *masque)
{

        size_t NbPixROI2d=dimROI.x*dimROI.y;
        nbCplx *holo=new nbCplx[NbPixROI2d];
        nbCplx *holo_shift=new nbCplx[NbPixROI2d];
        nbCplx *TF_Holo=new nbCplx[NbPixROI2d];
        nbCplx *TF_Holo_centre=new nbCplx[NbPixROI2d];

        for(size_t pixel=0; pixel<NbPixROI2d; pixel++) {
                holo[pixel].Re=(double)holo1[pixel];//*masque[pixel];
        }
        //SAV_Re(holo, NbPixROI2d, "/home/mat/tomo_test/holo.bin", FLOAT,"a+b");
        ///--------Circshift et TF2D HOLOGRAMME------
       //
        circshift2DCplx(holo,holo_shift, dimROI,decalROI);

       //
       clock_t t_init_TF, t_fin_TF;
       double t_total_TF;
        //cout<<"circshift="<<t_total<<endl;

      // t_init_TF=clock();
                TF2Dcplx(holo_shift,TF_Holo,dimROI);
               // t_fin_TF=clock();
               // t_total_TF=(double)(t_fin_TF-t_init_TF)/CLOCKS_PER_SEC;


      //  cout<<"TF1024="<<t_total_TF<<endl;
        circshift2DCplx(TF_Holo, TF_Holo_centre, dimROI, decalROI);   //Décalage  sur fft_reel_tmp, pour recentrer le spectre avant découpe (pas obligatoire mais plus clair)

        ///--------Découpe hors axée------------------
        coupeCplx(TF_Holo_centre, TF_UBornTot, dimROI, dim2DHA, coinHA);///Découpe à [-Nxmax,+NXmax]
   delete[] holo, delete[] holo_shift, delete[] TF_Holo, delete[] TF_Holo_centre;
}

void calcPhase2pi(nbCplx* obj, Var2D taille,double* phaseMod2pi)///calcul phase -PI-PI
{ ///calcul phase de 0 à 2 PI
        size_t image_size=taille.x*taille.y;

     for(size_t pixel=0;pixel<image_size;pixel++)
     {
        double cos_phase=obj[pixel].Re;
        double sin_phase=obj[pixel].Im;

        if(sin_phase>0)
        {
            if(cos_phase>0)
            {
                double argTanPhase=sin_phase/cos_phase;
                phaseMod2pi[pixel]=atan(argTanPhase);
            }
            else if(cos_phase<0)
            {
                double argTanPhase=sin_phase/cos_phase;
                phaseMod2pi[pixel]=atan(argTanPhase)+PI;
            }
            else if(cos_phase==0)
            {
                phaseMod2pi[pixel]=PI/2;
            }
        }
        else if(sin_phase<0)
        {
            if(cos_phase<0)
            {
                double argTanPhase=sin_phase/cos_phase;
                phaseMod2pi[pixel]=atan(argTanPhase)+PI;
            }
            else if(cos_phase>0)
            {
                double argTanPhase=sin_phase/cos_phase;
                phaseMod2pi[pixel]=atan(argTanPhase)+2*PI;
            }
            else if(cos_phase==0)
            {
                phaseMod2pi[pixel]=3*PI/2;
            }
        }
        else if(sin_phase==0)
        {
            if(cos_phase>0)
            {
                phaseMod2pi[pixel]=0;
            }
            else if(cos_phase<0)
            {
                phaseMod2pi[pixel]=PI;
            }

        }
     }

}

void calc_Uborn(nbCplx *TF_UBorn,nbCplx *UBorn,Var2D dim2DHA,Var2D PosSpec)
{
    clock_t t_init,t_fin;
    double t_tot;
    Var2D recalUBorn={-PosSpec.x,-PosSpec.y},DecalU_Born={dim2DHA.x/2,dim2DHA.y/2};
    size_t NbPixU_Born=dim2DHA.x*dim2DHA.y;
    nbCplx TF_UBorn_I[NbPixU_Born];
    nbCplx UBorn_I[NbPixU_Born];
    decal2DCplxGen(TF_UBorn,TF_UBorn_I,dim2DHA,recalUBorn);
    TF2Dcplx_INV(TF_UBorn_I,UBorn_I,dim2DHA);
    t_init=clock();
    circshift2DCplx(UBorn_I,UBorn,dim2DHA,DecalU_Born);
    t_fin=clock();
    t_tot=(double)(t_fin-t_init)/CLOCKS_PER_SEC;
   // cout<<"temps cirshift 1 angle="<<t_tot<<endl;

}

void Chrono(temps *t, string message){
    t->fin = clock ();
    float temps_cpu = (t->fin - t->init) * 1e-6;
    t->total= t->total+temps_cpu;
   // cout<<"coucou total="<<t->fin<<endl;
    t->init = clock();
    cout<<message<<endl;

}
int chargeBin(float *objet, string chemin,  int NbPix)
{
    size_t precision=sizeof(objet[0]);//attention si on n'indique pas "0", sizeof donnera la taille du pointeur (64bits)!
    FILE *fichier_ID = NULL;
    fichier_ID= fopen(chemin.c_str(),"rb");
    if(fichier_ID==0)
    cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;
    cout<<"precision en nb d'octet="<<precision<<endl;
    fread(objet,precision,NbPix,fichier_ID);
    return 1;
}

void divCplx(nbCplx* imageA,nbCplx* imageB,nbCplx* resultat, int NbPix)
{
    for(int cpt=0;cpt<NbPix;cpt++){
        double modB_Q=pow(imageB[cpt].Re,2)+pow(imageB[cpt].Im,2);//Q pour quadratique=carré
        resultat[cpt].Re=double(imageA[cpt].Re*imageB[cpt].Re+imageA[cpt].Im*imageB[cpt].Im)/modB_Q;
        resultat[cpt].Im=double(imageA[cpt].Im*imageB[cpt].Re-imageA[cpt].Re*imageB[cpt].Im)/modB_Q;
    }
}
void multiplierCplx(nbCplx* image1,nbCplx* image2,nbCplx* resultat, int NbPix)
{
    for(int cpt=0;cpt<NbPix;cpt++)
    {
    resultat[cpt].Re=image1[cpt].Re*image2[cpt].Re-image1[cpt].Im*image2[cpt].Im;
    resultat[cpt].Im=image1[cpt].Re*image2[cpt].Im+image1[cpt].Im*image2[cpt].Re;
    }
}
///recalage par correlation croisée
void recale(nbCplx* obj,nbCplx* objDecal,nbCplx *objRecal, Var3D dimVol)
{
    unsigned int Npix3D=dimVol.x*dimVol.y*dimVol.z, cpt=0;
    nbCplx *A=new nbCplx[Npix3D];
    nbCplx *B=new nbCplx[Npix3D];
    nbCplx *C=new nbCplx[Npix3D];
    nbCplx *R=new nbCplx[Npix3D];
    TF3DCplx(obj, A, dimVol);
     TF3DCplx(objDecal, B, dimVol);

///Calcul de R, correlation de phase dans Fourier=A*conj(B)/(modA.modB)
   for(cpt=0;cpt<Npix3D;cpt++)
        {
            double modA=sqrt( pow(A[cpt].Re,2)+pow(A[cpt].Im,2));
            double modB=sqrt( pow(B[cpt].Re,2)+pow(B[cpt].Im,2));
            if(modA==0)
                modA=1;
            if(modB==0)
                modB=1;
            R[cpt].Re= (A[cpt].Re*B[cpt].Re + A[cpt].Im*B[cpt].Im)/(modA*modB);
            R[cpt].Im=(-B[cpt].Im*A[cpt].Re + A[cpt].Im*B[cpt].Re)/(modA*modB);
        }

    delete[] A;

    ///Recalage en multipliant B Par R
   for(cpt=0;cpt<Npix3D;cpt++)
        {
            C[cpt].Re=B[cpt].Re*R[cpt].Re  -  R[cpt].Im*B[cpt].Im;
            C[cpt].Im=B[cpt].Im*R[cpt].Re  +  B[cpt].Re*R[cpt].Im;
        }
    TF3DCplx_INV(C, objRecal, dimVol);

    delete[] B;
    delete[] C;
    delete[] R;

}

void mod_cplx(nbCplx *A, double *module, int NbPix3D)
{

       for(int cpt=0;cpt<NbPix3D;cpt++)
    {
        module[cpt]=sqrt( pow(A[cpt].Re,2)+pow(A[cpt].Im,2) );

    }

}



Var2D corr_crois(double* obj2D_A, double* obj2D_B, Var2D dim)
{
    int NbPts2D=dim.x*dim.y;

//copie pour passer en complexe...
    nbCplx *copie_obj2D_A=new nbCplx[NbPts2D];
    nbCplx *copie_obj2D_B=new nbCplx[NbPts2D];
        for(int cpt=0;cpt<NbPts2D;cpt++)//copies pour passer en complexe (la fonction Tf2D demandent des cplx)
    {
        copie_obj2D_A[cpt].Re=obj2D_A[cpt];
        copie_obj2D_A[cpt].Im=0;
        copie_obj2D_B[cpt].Re=obj2D_B[cpt];
        copie_obj2D_B[cpt].Im=0;
    }
    //stockage des spectres
    nbCplx *spect2D_A= new nbCplx[NbPts2D];
    nbCplx *spect2D_B= new nbCplx[NbPts2D];
    //stocker la correlation et son spectre
    nbCplx *spectCorr=new nbCplx[NbPts2D];
    nbCplx *Corr=new nbCplx[NbPts2D];


    //SAV_Re(copie_obj2D_B,NbPts2D,"/home/mat/tomo_test/copie_objetB.bin",FLOAT,"wb");
    TF2Dcplx(copie_obj2D_A,spect2D_A,dim);
    TF2Dcplx(copie_obj2D_B,spect2D_B,dim);

//Elimination du module pour isoler  la phase : A*conj(B)/(modA*modB)

for(int cpt=0;cpt<NbPts2D;cpt++){
    //double modA=sqrt(pow(spect2D_A[cpt].Re,2)+pow(spect2D_A[cpt].Im,2));
    //double modB=sqrt(pow(spect2D_B[cpt].Re,2)+pow(spect2D_B[cpt].Im,2));
    spectCorr[cpt].Re=double(spect2D_A[cpt].Re*spect2D_B[cpt].Re+spect2D_A[cpt].Im*spect2D_B[cpt].Im);
    spectCorr[cpt].Im=double(spect2D_A[cpt].Im*spect2D_B[cpt].Re-spect2D_A[cpt].Re*spect2D_B[cpt].Im);
    }
    //TF de l'expoentielle contenant le dépahsage->translation dans l'espace direct
    TF2Dcplx_INV(spectCorr,Corr,dim);
    int cptMax=0;
    double valMax=0;
    double valModCorr=0;
    for(int cpt=0;cpt<NbPts2D;cpt++)
    {
     valModCorr=pow(Corr[cpt].Re,2)+pow(Corr[cpt].Im,2);

    if(valModCorr>valMax)
        {
        valMax=valModCorr;
            cptMax=cpt;
        }
    }
    Var2D decal2D_I={0,0};
    decal2D_I.x=cptMax%dim.x;
    decal2D_I.y=cptMax/dim.y;
    if(decal2D_I.x>15 && decal2D_I.x<dim.x-15)
    decal2D_I.x=0;
     if( decal2D_I.y>15 &&  decal2D_I.y<dim.y-15)
    decal2D_I.y=0;
    cout<<"---------------------------"<<endl;
    cout<<"X="<<decal2D_I.x<<endl;
    cout<<"Y="<<decal2D_I.y<<endl;
   // SAV_Re(spectCorr,NbPts2D,"/home/mat/tomo_test/SpectCorrRE.bin",FLOAT,"wb");
    //SAV_Re(Corr,NbPts2D,"/home/mat/tomo_test/Corr.bin",FLOAT,"wb");
    return decal2D_I;
}
double bruit(int attenuation)
{
    double valBruit=0;
    valBruit=((double)rand()/(double)RAND_MAX);

       return valBruit;
}



void TF3DCplx(nbCplx *Objet3D_shift, nbCplx* Spect3D_shift,Var3D dimVol)
{       int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        int nthreads=4;
        fftw_plan_with_nthreads(nthreads);
        int NbPix3D=dimVol.x*dimVol.y*dimVol.z;
        //Déclaration des variables pour la FFT : entre,sortie et "fftplan"
        fftw_complex *in3D, *out3D;


        //Réservation memoire
        in3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPix3D);
        ///Récupération de l'image dans la partie reelle de l'entree
        for(int cpt=0; cpt<(NbPix3D); cpt++) {
                in3D[cpt][0]=Objet3D_shift[cpt].Re;
                in3D[cpt][1]=Objet3D_shift[cpt].Im;
        }
        //Réservation memoire
        out3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPix3D);
        fftw_plan p3D;
        p3D=fftw_plan_dft_3d(dimVol.x, dimVol.y, dimVol.z, in3D, out3D,FFTW_FORWARD, FFTW_ESTIMATE);
        //calcul du plan, parametre servant a calculer et optimiser le FFT
        fftw_execute(p3D); // repeat as needed
        fftw_destroy_plan(p3D);
        fftw_free(in3D);
        void fftw_cleanup_threads(void);
          for(int cpt=0; cpt<(NbPix3D); cpt++) {
                Spect3D_shift[cpt].Re=out3D[cpt][0]/(NbPix3D);
                Spect3D_shift[cpt].Im=out3D[cpt][1]/(NbPix3D);
        }
}

void TF3DCplx_INV(nbCplx *Spect3D_shift, nbCplx* Objet3D_shift,Var3D dimVol)
{       int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        int nthreads=4;
        fftw_plan_with_nthreads(nthreads);
        int NbPix3D=dimVol.x*dimVol.y*dimVol.z;
        //Déclaration des variables pour la FFT : entre,sortie et "fftplan"
        fftw_complex *in3D, *out3D;


        //Réservation memoire entrée
        in3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPix3D);
        ///Récupération de l'image dans la partie reelle de l'entree
        for(int cpt=0; cpt<(NbPix3D); cpt++) {
                in3D[cpt][0]=Spect3D_shift[cpt].Re;
                in3D[cpt][1]=Spect3D_shift[cpt].Im;
        }
        //Réservation memoire sortie
        out3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NbPix3D);
        fftw_plan p3D;
        p3D=fftw_plan_dft_3d(dimVol.x, dimVol.y, dimVol.z, in3D, out3D,FFTW_BACKWARD, FFTW_ESTIMATE);
        //calcul du plan, parametre servant a calculer et optimiser la FFT
        fftw_execute(p3D); // repeat as needed
        fftw_destroy_plan(p3D);
        fftw_free(in3D);
        void fftw_cleanup_threads(void);
          for(int cpt=0; cpt<(NbPix3D); cpt++) {
                Objet3D_shift[cpt].Re=out3D[cpt][0];
                Objet3D_shift[cpt].Im=out3D[cpt][1];
        }
}

void genere_rectang2D(double *objet,Var2D posI_Coin,Var2D dimRect,Var2D dim)
{
       for(int y=posI_Coin.y;y<posI_Coin.y+dimRect.y;y++)
        {
         int NbPixY=dim.x*y;
          for(int x=posI_Coin.x;x<posI_Coin.x+dimRect.x;x++)
            {
          int cpt3D=NbPixY+x;
          objet[cpt3D]=1.0;//+bruit(1);
            }
        }
}

void genere_rectang3D(nbCplx *objet,Var3D posI_Coin,Var3D dimRect,Var3D dim)
{
     for(int z=posI_Coin.z;z<posI_Coin.z+dimRect.z;z++)
    {
        int altitude=dim.x*dim.y*z;
       for(int y=posI_Coin.y;y<posI_Coin.y+dimRect.y;y++)
     {
         int NbPixY=dim.x*y;
          for(int x=posI_Coin.x;x<posI_Coin.x+dimRect.x;x++)
      {
          int cpt3D=altitude+NbPixY+x;
          objet[cpt3D].Re=1.0;//+bruit(1);
          objet[cpt3D].Im=7.0;//+bruit(1);
      }
     }
    }


}
void genere_OTF_T_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon)
{
  int k=0, dimPlanFinal=round(dim_final.x*dim_final.y), Kdx0=posSpec.x,Kdy0=posSpec.y, dimVolX=dim_final.x;
         ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
                                double r2=rayon*rayon, arg_z_arc=0,z_arc=0,Kdz0;
                                double Kdz0_carre = rayon*rayon-Kdx0*Kdx0-Kdy0*Kdy0;
                                if(round(Kdz0_carre)>-1) {
                                        Kdz0=sqrt(Kdz0_carre);
                                        int NXMAX_CARRE=NMAX.x*NMAX.x;
                                        for (int Kdy = -NMAX.y; Kdy < NMAX.y; Kdy++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
                                                int Kdy_carre=Kdy*Kdy;
                                                for (int Kdx = -NMAX.x; Kdx < NMAX.x; Kdx++) { //on balaye l'image 2D en y, centre au milieu
                                                        //int cpt=(Kdy+NMAX.y)*2*NMAX.x+Kdx+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D

                                                        if(Kdx*Kdx+Kdy_carre<NXMAX_CARRE)
                                                        { //ne pas depasser l'ouverture numérique pour 1 hologramme
                                                                double Kdz_carre=r2-Kdx*Kdx-Kdy_carre; //altitude au carré des données
                                                                double Kdz=round(sqrt(Kdz_carre)-Kdz0);///-Kdz->réflexion
                                                                double altitude=(Kdz+decal.z)*dimPlanFinal; //donne n'importequoi sans l'arrondi sur z!!
                                                                k=(-Kdx0+Kdx+decal.x)+(-Kdy0+Kdy+decal.y)*dimVolX+round(altitude);//indice du tableau 1D du volume 3D
                                                                OTFr[k].Re=1;//
                                                                OTFr[k].Im=1;//e
                                                        }
                                                } //fin for y
                                        }
                                }//fin if zm0>-1
}
void genere_OTF_RB_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon)
{
  int k=0, dimPlanFinal=round(dim_final.x*dim_final.y), Kdx0=posSpec.x,Kdy0=posSpec.y, dimVolX=dim_final.x;
         ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
                                double r2=rayon*rayon, arg_z_arc=0,z_arc=0,Kdz0;
                                double Kdz0_carre = rayon*rayon-Kdx0*Kdx0-Kdy0*Kdy0;
                                if(round(Kdz0_carre)>-1) {
                                        Kdz0=sqrt(Kdz0_carre);
                                        int NXMAX_CARRE=NMAX.x*NMAX.x;
                                        for (int Kdy = -NMAX.y; Kdy < NMAX.y; Kdy++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
                                                int Kdy_carre=Kdy*Kdy;
                                                for (int Kdx = -NMAX.x; Kdx < NMAX.x; Kdx++) { //on balaye l'image 2D en y, centre au milieu
                                                        int cpt=(Kdy+NMAX.y)*2*NMAX.x+Kdx+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D

                                                        if(Kdx*Kdx+Kdy_carre<NXMAX_CARRE)
                                                        { //ne pas depasser l'ouverture numérique pour 1 hologramme
                                                                double Kdz_carre=r2-Kdx*Kdx-Kdy_carre; //altitude au carré des données
                                                                double Kdz=round(sqrt(Kdz_carre)+Kdz0);///
                                                                double altitude=(Kdz+decal.z)*dimPlanFinal; //donne n'importequoi sans l'arrondi sur z!!
                                                                k=(+Kdx0-Kdx+decal.x)+(Kdy0-Kdy+decal.y)*dimVolX+round(altitude);//indice du tableau 1D du volume 3D
                                                                OTFr[k].Re=1;//
                                                                OTFr[k].Im=1;//e
                                                        }
                                                } //fin for y
                                        }
                                }//fin if zm0>-1
}
void genere_OTF_RH_Holo(nbCplx *OTFr, Var2D posSpec, Var3D dim_final, Var3D decal, Var2D NMAX, double rayon)
{
  int k=0, dimPlanFinal=round(dim_final.x*dim_final.y), Kdx0=posSpec.x,Kdy0=posSpec.y, dimVolX=dim_final.x;
         ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
                                double r2=rayon*rayon, arg_z_arc=0,z_arc=0,Kdz0;
                                double Kdz0_carre = rayon*rayon-Kdx0*Kdx0-Kdy0*Kdy0;
                                if(round(Kdz0_carre)>-1) {
                                        Kdz0=sqrt(Kdz0_carre);
                                        int NXMAX_CARRE=NMAX.x*NMAX.x;
                                        for (int Kdy = -NMAX.y; Kdy < NMAX.y; Kdy++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
                                                int Kdy_carre=Kdy*Kdy;
                                                for (int Kdx = -NMAX.x; Kdx < NMAX.x; Kdx++) { //on balaye l'image 2D en y, centre au milieu
                                                        int cpt=(Kdy+NMAX.y)*2*NMAX.x+Kdx+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D

                                                        if(Kdx*Kdx+Kdy_carre<NXMAX_CARRE)
                                                        { //ne pas depasser l'ouverture numérique pour 1 hologramme
                                                                double Kdz_carre=r2-Kdx*Kdx-Kdy_carre; //altitude au carré des données
                                                                double Kdz=round(-sqrt(Kdz_carre)-Kdz0);///-Kdz->réflexion
                                                                double altitude=(Kdz+decal.z)*dimPlanFinal; //donne n'importequoi sans l'arrondi sur z!!
                                                                k=(-Kdx0+Kdx+decal.x)+(-Kdy0+Kdy+decal.y)*dimVolX+round(altitude);//indice du tableau 1D du volume 3D
                                                                OTFr[k].Re=1;//
                                                                OTFr[k].Im=1;//e
                                                        }
                                                } //fin for y
                                        }
                                }//fin if zm0>-1
}

void multiplier_masqueCplx2(nbCplx *image, nbCplx *masque, int t_image, int t_mask, Var2D CentreI)
{
        int t_imageX=t_image;
        int t_imageY=t_image;
        int y=0;
        int x=0;
        int xMaskInf=round(t_mask/2);
        int yMaskInf=round(t_mask/2);
        //if((t_mask%2)!=0)
        //cout<<"fonction multiplier_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;
        if(CentreI.x<xMaskInf || CentreI.y<yMaskInf)
        cout<<"Débordement par le bord supérieur. ("<<CentreI.x<<","<<CentreI.y<<")"<<endl;
         if(t_image-CentreI.x>xMaskInf || t_image-CentreI.y<yMaskInf)
        cout<<"Débordement! Coordonnées supérieures du masque=("<<CentreI.x<<","<<CentreI.y<<")"<<endl;

        for(int xMask=-xMaskInf; xMask<xMaskInf; xMask++) {
                for(int yMask=-yMaskInf; yMask<yMaskInf; yMask++) {
                        //si le masque déborde de l'image (attention le masque repasse du coté opposé
                        ///vers coord info pour le masque
                        int yiMask=yMask+yMaskInf;
                        int xiMask=xMask+xMaskInf;
                        ///compteur 1D des coordonnées
                         y=yiMask+CentreI.y;
                         x=xiMask+CentreI.x;
                        int cptj=y*t_image+x;
                        int cptMask=yiMask*t_mask+xiMask;
                        cout<<cptMask<<endl;
                        //masque[cptMask].Re=0;
                        image[cptj].Re=image[cptj].Re*masque[cptMask].Re - image[cptj].Im*masque[cptMask].Im;
                        image[cptj].Im=image[cptj].Re*masque[cptMask].Im + image[cptj].Im*masque[cptMask].Re;
                }
        }
}

void Plan_ds_VolCplx(nbCplx *Vol3D, nbCplx *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di)
{
   int x3Di=0, y3Di=0;
   int taillePlan=dimPlan.x*dimPlan.y;
    int planVol3D=dimVol.x*dimVol.y;
    int ecartX=(dimVol.x-dimPlan.x)/2, ecartY=(dimVol.y-dimPlan.y)/2;
    int altitude_k=z3Di*planVol3D;
    int cpt2D=0, k=0;

    for (int yi = 0; yi < dimPlan.y; yi++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
              for (int xi = 0; xi < dimPlan.x; xi++) { //on balaye l'image 2D en y, centre au milieu
                    cpt2D=yi*dimPlan.x+xi;//calcul du cpt du tableau 1D de l'image 2D
                    x3Di=xi+ecartX;
                    y3Di=yi+ecartY;
                    k=altitude_k+y3Di*dimVol.x+x3Di;
                    Vol3D[k].Re=plan2D[cpt2D].Re;
                    Vol3D[k].Im=plan2D[cpt2D].Im;
              }
    }
}
void Plan_ds_Vol(double *Vol3D, double *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di)
{
   int x3Di=0, y3Di=0;
   int taillePlan=dimPlan.x*dimPlan.y;
    int planVol3D=dimVol.x*dimVol.y;
    int ecartX=(dimVol.x-dimPlan.x)/2, ecartY=(dimVol.y-dimPlan.y)/2;
    int altitude_k=z3Di*planVol3D;
    int cpt2D=0, k=0;

    for (int yi = 0; yi < dimPlan.y; yi++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
              for (int xi = 0; xi < dimPlan.x; xi++) { //on balaye l'image 2D en y, centre au milieu
                    cpt2D=yi*dimPlan.x+xi;//calcul du cpt du tableau 1D de l'image 2D
                    x3Di=xi+ecartX;
                    y3Di=yi+ecartY;
                    k=altitude_k+y3Di*dimVol.x+x3Di;
                    Vol3D[k]=plan2D[cpt2D];

              }
    }
}
void Vol_ds_Plan(double *Vol3D, double *plan2D, Var3D dimVol, Var2D dimPlan, int z3Di)
{
    int x3Di=0, y3Di=0;
    int taillePlan=dimPlan.x*dimPlan.y;
    int planVol3D=dimVol.x*dimVol.y;
    int ecartX=(dimVol.x-dimPlan.x)/2, ecartY=(dimVol.y-dimPlan.y)/2;
    int altitude_k=z3Di*planVol3D;
    int cpt2D=0, k=0;

    for (int yi = 0; yi < dimPlan.y; yi++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
              for (int xi = 0; xi < dimPlan.x; xi++) { //on balaye l'image 2D en y, centre au milieu

                    cpt2D=yi*dimPlan.x+xi;//calcul du cpt du tableau 1D de l'image 2D
                    x3Di=xi+ecartX;
                    y3Di=yi+ecartY;
                    k=altitude_k+y3Di*dimVol.x+x3Di;
                    plan2D[cpt2D]=Vol3D[k];
              }
    }
}


void  retroPropagSA(int deltaZ, nbCplx *fft_shift_norm, nbCplx * planObjet, Var3D decal, Var2D NMAX, double rayon)
 {
    double kz[2*NMAX.x*2*NMAX.y];
    nbCplx rephase[2*NMAX.x*2*NMAX.y];
    nbCplx spectrePropag[2*NMAX.x*2*NMAX.y];
    int dimPlan=round(2*NMAX.x*2*NMAX.y);

         ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
                                double r2=rayon*rayon;
                                double arg_z_arc=0,z_arc=0;
                                //printf("round(rayon*rayon-(xm0)^2-(ym0)^2: %i\n",round(rayon*rayon-(xm0)^2-(ym0)^2));
                                //double zm0_carre = rayon*rayon-xm0*xm0-ym0*ym0;

                                        int NXMAX_CARRE=NMAX.x*NMAX.x;
                                            ///CALCULER LES KZ-----------------------------------
                                        for (int y = -NMAX.y; y < NMAX.y; y++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
                                                int y_carre=y*y;
                                                for (int x = -NMAX.x; x < NMAX.x; x++) { //on balaye l'image 2D en y, centre au milieu
                                                        int cpt=(y+NMAX.y)*2*NMAX.x+x+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D
                                                        if(x*x+y_carre<NXMAX_CARRE)
                                                        { //ne pas depasser l'ouverture numérique pour 1 hologramme
                                                                double z_carre=r2-x*x-y_carre; //altitude au carré des données
                                                                int cpt=y*(NMAX.x)*2+x;
                                                                kz[cpt]=sqrt(z_carre);
                                                        }
                                                        else  {
                                                               int cpt=y*(NMAX.x)*2+x;
                                                                kz[cpt]=0;
                                                        }
                                                } //fin for y
                                        }
                                        ///---------------------CALCULER LE SPECTRE REPROPAGE-----------------------------------

                                                         for(int pix=0;pix<dimPlan;pix++)
                                                         {  int xc=0,yc=0;
                                                            rephase[pix].Re=cos(kz[pix]*deltaZ);
                                                            rephase[pix].Im=sin(kz[pix]*deltaZ);
                                                            spectrePropag[pix].Re=fft_shift_norm[pix].Re*rephase[pix].Re-fft_shift_norm[pix].Im*rephase[pix].Im;
                                                            spectrePropag[pix].Im=fft_shift_norm[pix].Re*rephase[pix].Im+fft_shift_norm[pix].Im*rephase[pix].Re;
                                                            //xc=pix%(2*NMAX.x)-NMAX.x, yc=pix/(2*NYMAX)-NMAX.y;
                                                         }
                                                          TF2Dcplx_INV(spectrePropag, planObjet, NMAX);

 }

void decalCoupeCplx(nbCplx *fft, nbCplx *fft_tmp, Var2D NMAX,Var2D dimCCD)
{
for (int xi=0; xi<NMAX.x; xi++) {
                                for (int yi=0; yi<NMAX.y; yi++) {
                                        int cpt1=yi*dimCCD.x+xi;
                                        int cpt2=yi*(NMAX.x)*2+xi;

                                        fft[cpt2].Re=fft_tmp[cpt1].Re;
                                        fft[cpt2].Im=fft_tmp[cpt1].Im;
                                }
                                for (int yi=dimCCD.y-NMAX.y; yi<dimCCD.y; yi++) {
                                        int cpt1=yi*dimCCD.x+xi;

                                        int cpt2 = xi+(yi-dimCCD.y+2*NMAX.y)*2*(NMAX.x);
                                        fft[cpt2].Re=fft_tmp[cpt1].Re;
                                        fft[cpt2].Im=fft_tmp[cpt1].Im;
                                }
                        }
                        ///---////////////////////////////////////////////deuxieme demi-espace
                        for(int xi=dimCCD.x-NMAX.x; xi<dimCCD.x; xi++) {
                                for (int yi=0; yi<NMAX.y; yi++) {
                                        int cpt1=yi*dimCCD.x+xi;
                                        int cpt2=yi*(NMAX.x)*2+(xi-dimCCD.x+2*NMAX.x);

                                        fft[cpt2].Re=fft_tmp[cpt1].Re;
                                        fft[cpt2].Im=fft_tmp[cpt1].Im;
                                }

                                for (int yi=dimCCD.y-NMAX.y; yi<dimCCD.y; yi++) {
                                        int cpt1=yi*dimCCD.x+xi;

                                        int cpt2 = (xi-dimCCD.x+2*NMAX.x)+(yi-dimCCD.y+2*NMAX.y)*2*(NMAX.x);
                                        fft[cpt2].Re=fft_tmp[cpt1].Re;
                                        fft[cpt2].Im=fft_tmp[cpt1].Im;
                                }
                        }
}

void decalCoupe(double *fft_reel, double *fft_imag, double *fft_reel_tmp, double *fft_imag_tmp, Var2D NMAX,Var2D dimCCD)
{
for (int xi=0; xi<NMAX.x; xi++) {
                                for (int yi=0; yi<NMAX.y; yi++) {
                                        int cpt1=yi*dimCCD.x+xi;
                                        int cpt2=yi*(NMAX.x)*2+xi;

                                        fft_reel[cpt2]=fft_reel_tmp[cpt1];
                                        fft_imag[cpt2]=fft_imag_tmp[cpt1];
                                }
                                for (int yi=dimCCD.y-NMAX.y; yi<dimCCD.y; yi++) {
                                        int cpt1=yi*dimCCD.x+xi;

                                        int cpt2 = xi+(yi-dimCCD.y+2*NMAX.y)*2*(NMAX.x);
                                        fft_reel[cpt2]=fft_reel_tmp[cpt1];
                                        fft_imag[cpt2]=fft_imag_tmp[cpt1];
                                }
                        }
                        ///---////////////////////////////////////////////deuxieme demi-espace
                        for(int xi=dimCCD.x-NMAX.x; xi<dimCCD.x; xi++) {
                                for (int yi=0; yi<NMAX.y; yi++) {
                                        int cpt1=yi*dimCCD.x+xi;
                                        int cpt2=yi*(NMAX.x)*2+(xi-dimCCD.x+2*NMAX.x);

                                        fft_reel[cpt2]=fft_reel_tmp[cpt1];
                                        fft_imag[cpt2]=fft_imag_tmp[cpt1];
                                }

                                for (int yi=dimCCD.y-NMAX.y; yi<dimCCD.y; yi++) {
                                        int cpt1=yi*dimCCD.x+xi;

                                        int cpt2 = (xi-dimCCD.x+2*NMAX.x)+(yi-dimCCD.y+2*NMAX.y)*2*(NMAX.x);
                                        fft_reel[cpt2]=fft_reel_tmp[cpt1];
                                        fft_imag[cpt2]=fft_imag_tmp[cpt1];
                                }
                        }
}
void InitTabCplx(nbCplx *z,int taille)//initailiser un tableau (taille totale="taille") de structure à zéro.
{
    for(int cpt=0;cpt<taille;cpt++)
    {
        z[cpt] = {0}; /* Tous les champs à zéro */
    }
}



int retroPropag_Born(nbCplx *TF3D_PotObj, nbCplx *TF_Uborn_norm, double * sup_redon, int dim_final, Var2D posSpec, Var3D decal3D, Var2D NMAX, double rayon)
{
                        int kxmi=posSpec.x, kymi=posSpec.y;
                        int kxm0=(kxmi-NMAX.x), kym0=(kymi-NMAX.x);//coordonnée dans l'image2D centrée (xm0,ym0)=(0,0)=au centre de l'image
                        int points_faux=0;
                        int dimVolX=round(dim_final), dimPlanFinal=round(dim_final*dim_final);
                        float n0=1.515, lambda=pow(633,10^(-9)),kv=2*3.1416/lambda;

                                //création de variable pou-9r éviter N calculs dans la boucle sur le volume 3D
                                int cptPot=0; //indice tableau 1d des données du potentiel3D
                                double r2=rayon*rayon, kzm0, kzm0_carre = rayon*rayon-kxm0*kxm0-kym0*kym0;
                                //printf("round(rayon*rayon-(xm0)^2-(ym0)^2: %i\n",round(rayon*rayon-(xm0)^2-(ym0)^2));

                                if(round(kzm0_carre)>-1) {

                                        kzm0=sqrt(kzm0_carre);


                                        int NMAX_CARRE=NMAX.x*NMAX.x;
                                        float ctePotUb=1;//2*kv*n0;
                                        for (int ky = -NMAX.y; ky < NMAX.y; ky++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
                                                int ky_carre=ky*ky;
                                                for (int kx = -NMAX.x; kx < NMAX.x; kx++) { //on balaye l'image 2D en y, centre au milieu
                                                        int cpt=(ky+NMAX.y)*2*NMAX.x+kx+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D
                                                        if(kx*kx+ky_carre<NMAX_CARRE) { //ne pas depasser l'ouverture numérique pour 1 hologramme
                                                                double kz_carre=r2-kx*kx-ky_carre; //altitude au carré des données
                                                                double kz=round(sqrt(kz_carre)-kzm0);
                                                                double altitude=(kz+decal3D.z)*dimPlanFinal; //donne n'importequoi sans l'arrondi sur z!!

                                                                cptPot=(-kxm0+kx+decal3D.x)+(-kym0+ky+decal3D.y)*dimVolX+round(altitude);//indice du tableau 1D du volume 3D
                                                                //cout<<"k"<<k<<endl;
                                                                TF3D_PotObj[cptPot].Re+=-ctePotUb*TF_Uborn_norm[cpt].Im;//Inversion Re->Im à cause du coefficient i entre potentiel et Uborn
                                                                TF3D_PotObj[cptPot].Im+=ctePotUb*TF_Uborn_norm[cpt].Re;//
                                                                sup_redon[cptPot]+=1;//pour calculer le support
                                                        } else
                                                                points_faux++;
                                                } //fin for y
                                        }
                                }//fin if zm0>-1
}



double max(double *entree, int tailleTab)
{
    double valMax=0;
    for (int cpt=0; cpt<tailleTab;cpt++)
    {
    if(entree[cpt]>valMax)
    valMax=entree[cpt];
    }
    return valMax;
}
int coordMaxMod2D(nbCplx *entree, int tailleTab)
{
    double moduleQ[tailleTab];
 for(int cpt=0;cpt<tailleTab;cpt++)
 {
     moduleQ[cpt]=pow(entree[cpt].Re,2)+pow(entree[cpt].Im,2);
 }
 int cptMax=0;
 double valMax=0;
  for (int cpt=0; cpt<tailleTab;cpt++)
    {
    if(moduleQ[cpt]>valMax)
        {
             cptMax=cpt;
             valMax=moduleQ[cpt];
        }
    }
    return cptMax;
}

///####################fonction################"
void SAV(double *var_sav, int NbPix2D, char *chemin, enum PRECISION precision, char options[])
{
        FILE *fichier_ID;
        fichier_ID= fopen(chemin, options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision) {
        case DOUBLE: //64 bit

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        double tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case FLOAT://32 bits float

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        float tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;

        case INT: //32 bit signé

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case UINT://32 bit non signé

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        unsigned int tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case CHAR: //8 bits

                for(unsigned int cpt=0; cpt<NbPix2D; cpt++) {
                        char tampon=var_sav[cpt];
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        default:
                break;
        }

        fclose(fichier_ID);
}

void SAV_Re(nbCplx *var_sav, int taille, char *chemin, enum PRECISION precision, char options[])
{
        FILE *fichier_ID;
        fichier_ID= fopen(chemin, options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision) {
        case DOUBLE: //64 bit

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        double tampon=var_sav[cpt].Re;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case FLOAT://32 bits float

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        float tampon=var_sav[cpt].Re;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;

        case INT: //32 bit signé

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        int tampon=var_sav[cpt].Re;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case UINT://32 bit non signé

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        unsigned int tampon=var_sav[cpt].Re;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case CHAR: //8 bits

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        char tampon=var_sav[cpt].Re;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        default:
                break;
        }

        fclose(fichier_ID);
}
void SAV_Im(nbCplx *var_sav, int taille, char *chemin, enum PRECISION precision, char options[])
{
        FILE *fichier_ID;
        fichier_ID= fopen(chemin, options);
        if(fichier_ID==0)
                cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

        switch(precision) {
        case DOUBLE: //64 bit

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        double tampon=var_sav[cpt].Im;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case FLOAT://32 bits float

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        float tampon=var_sav[cpt].Im;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;

        case INT: //32 bit signé

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        int tampon=var_sav[cpt].Im;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case UINT://32 bit non signé

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        unsigned int tampon=var_sav[cpt].Im;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        case CHAR: //8 bits

                for(unsigned int cpt=0; cpt<taille; cpt++) {
                        char tampon=var_sav[cpt].Im;
                        fwrite(&tampon,sizeof(tampon),1,fichier_ID);
                }
                break;
        default:
                break;
        }

        fclose(fichier_ID);
}
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
void changeDim2DCplx(nbCplx *tab, nbCplx* tabFinal, Var2D dimInit, Var2D dimFin)
{
        int diffX=round(dimFin.x-dimInit.x)/2;
        int diffY=round(dimFin.y-dimInit.y)/2;

        for(int x=0; x<dimInit.x; x++) {
                for( int y=0; y<dimInit.y; y++) {
                        int cpt_init=y*dimInit.x+x;
                        int cpt_final=(y+diffY)*dimFin.x+x+diffX;
                        tabFinal[cpt_final].Re=tab[cpt_init].Re;
                        tabFinal[cpt_final].Im=tab[cpt_init].Im;
                }
        }
}



/*void concatener(char chaine1[], char chaine2[], char resultat[])//fonction limitée mais intéret car copie locale (pas besoin de tampon)
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

}*/

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
        double *  tuk2D = new double[dimx*dimy];
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

void ecrire_rapport(int NXMAX,float rayon,float Rf,  int DIMX_CCD2,int coin_x, int coin_y,short int precision_exportation,string chemin,int nb_proj,float n1,float NA,float Tp, int G)
{
        time_t date;
        time(&date);
        char nom_rapport[12];
        string nom_rapport2=chemin+"/rapport.txt";
//        concatener(chemin,"/rapport.txt",nom_rapport);
        FILE *fichier_rapport ;
        cout<<"nom_rapport:" <<nom_rapport2<<endl;
        /*  ouverture pour ecriture (w) en mode texte (t) */
        fichier_rapport = fopen (nom_rapport2.c_str(), "wt") ;
        if (fichier_rapport == NULL)
                printf ("impossible de créer le fichier rapport_calcul.txt\n");
        //fprintf(fichier_rapport,"Date     : %s\nNXMAX    : %i\n,Rayon : %i\nRf       : %f\nK : %f\ndimx_ccd : %d\ncoin_x   : %d\ncoin_y   : %d\nPrecision: %i bits\nSession  : %s\nnb_proj  : %i\nindice n1: %f\nNA       : %f\nT_pixel  : %e\nG        : %i",ctime(&date),NXMAX,rayon,Rf,K,DIMX_CCD2,coin_x,coin_y,8*sizeof(precision),chemin,nb_proj,n1,NA,Tp,G);
        fprintf(fichier_rapport,"Date : %s\n NXMAX=%i\n Rf=%f\n Precision: %i bits\n Session  : %s\n Nombre de projections : %i",ctime(&date),NXMAX,Rf,8*precision_exportation,chemin.c_str(),nb_proj);
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


void charger_image2D(unsigned char* phasei, string imgFile, Var2D coin,Var2D taille)
{

        //rempli_tableau(phasei, imgFile, coin,taille);
        Image Monimage;
        //int i, currentImageWidth, currentImageHeight;
        Monimage.read(imgFile);////// chargement en memoire de l'image
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
        Monimage.writePixels(GrayQuantum, phasei);

}

/////////////////////////////////////////////////////////////////////////////////////
//Fonction transferant les pixels de l'image 2D de taille (nx,ny) dans un tableau 1D de taille (nx*ny)
// nécessite la bibliotheque imageMagick
/////////////////////////////////////////////////////////////////////////////////////

/*void rempli_tableau(unsigned char *finalArray, string path, Var2D coin, Var2D taille)
{

}*/

/////////////////////////////////////////////////////////////////////////////////////
//fonction de calcul TF2D
/////////////////////////////////////////////////////////////////////////////////////
void TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y)
{       int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        fftw_plan_with_nthreads(4);
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

                fft_reel[cpt]=out[cpt][0]/N; //division par N^2 pour normaliser la fftw qui n'est pas normalisée
                fft_imag[cpt]=out[cpt][1]/N;
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
}


void prepare_wisdom(Var2D dim, char *chemin)
{       int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        fftw_plan_with_nthreads(4);
        int N=dim.x*dim.y;

        fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
        fftw_plan p;
        //Réservation memoire
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_FORWARD, FFTW_EXHAUSTIVE);
        fftw_export_wisdom_to_filename(chemin);
        fftw_destroy_plan(p);
}


void TF2Dcplx(nbCplx *entree, nbCplx *fft, Var2D dim)
{
        int fftwThreadInit;
        fftwThreadInit=fftw_init_threads();
        //int nthreads=3;
        fftw_plan_with_nthreads(4);
        int N=dim.x*dim.y;
        fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
        fftw_plan p;
        //Réservation memoire
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        //Récupération de l'image dans la partie reelle de l'entree
        for(int cpt=0; cpt<N; cpt++) {
                in[cpt][0]=entree[cpt].Re;
                in[cpt][1]=entree[cpt].Im;
        }
        //calcul du plan, parametre servant a calculer et optimiser le FFT
       p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
           //p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_FORWARD,FFTW_WISDOM_ONLY);
        fftw_execute(p); /* repeat as needed */

        for(int cpt=0; cpt<(N); cpt++) {
                fft[cpt].Re=out[cpt][0]/N; //division par N (dim*dim) pour normaliser la fftw qui n'est pas normalisée
                fft[cpt].Im=out[cpt][1]/N;
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
}
void TF2Dcplx_INV(nbCplx *fft_entree, nbCplx *objet, Var2D dim)
{
        fftw_plan_with_nthreads(4);
        int N=dim.x*dim.y;
        fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
        fftw_plan p;
        //Réservation memoire
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        //Récupération de l'image dans la partie reelle de l'entree
        for(int cpt=0; cpt<N; cpt++) {
                in[cpt][0]=fft_entree[cpt].Re;
                in[cpt][1]=fft_entree[cpt].Im;
        }
        //calcul du plan, parametre servant a calculer et optimiser le FFT
        p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_BACKWARD, FFTW_ESTIMATE);
          // p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_BACKWARD,  FFTW_WISDOM_ONLY);
        fftw_execute(p); /* repeat as needed */

        for(int cpt=0; cpt<(N); cpt++) {
                objet[cpt].Re=out[cpt][0]; //pas de division par N^2 pour normaliser la fftw qui n'est pas normalisée
                objet[cpt].Im=out[cpt][1];
        }
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
}


/////////////////////////////////////////////////////////////////////////////////////
//fonction de circshift

/////////////////////////////////////////////////////////////////////////////////////

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
void decal2DGen(double* entree, double* result, Var2D dim,Var2D decal)
{
            //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
        decal.y=decal.y%dim.y;
        decal.x=decal.x%dim.x;
        if(decal.x<0)
        decal.x=dim.x+decal.x;
        if(decal.y<0)
        decal.y=dim.y+decal.y;
       for(int yi=0; yi<dim.y-decal.y; yi++) {
                for(int xi=0; xi<dim.x-decal.x; xi++)
                { //cout<<"xi,yi="<<xi<<","<<yi<<endl;
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
                        result[pixel_shift]=entree[pixel];
                }
                for(int xi=dim.x-decal.x; xi<dim.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
                        result[pixel_shift]=entree[pixel];
                }
        }
              for(int yi=dim.y-decal.y; yi<dim.y; yi++) {
                for(int xi=0; xi<dim.x-decal.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(-dim.y+yi+decal.y)*dim.x+xi+decal.x;
                        result[pixel_shift]=entree[pixel];
                }
                for(int xi=dim.x-decal.x; xi<dim.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(-dim.y+yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
                        result[pixel_shift]=entree[pixel];
                }
        }

}
void circshift2DCplx(nbCplx* entree, nbCplx* result, Var2D dim,Var2D decal)///___/!\ ne fonctionne qu'avec des demi espace??
{
        //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo

        decal.y=decal.y%dim.y;
        decal.x=decal.x%dim.x;
       // cout<<"decal.x="<<decal.x<<"; decal.y="<<decal.y<<endl;
       for(int yi=0; yi<decal.y; yi++) {
                for(int xi=0; xi<decal.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                      //  cout<<"pixel="<<pixel<<endl;
                       // cout<<"result[pixel]="<<result[pixel];
                        int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
                      //  cout<<"pixel_shift="<<pixel_shift<<endl;
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

void decal2DCplxGen(nbCplx* entree, nbCplx* result, Var2D dim,Var2D decal)
{
            //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
        decal.y=decal.y%dim.y;
        decal.x=decal.x%dim.x;
        if(decal.x<0)
        decal.x=dim.x+decal.x;
        if(decal.y<0)
        decal.y=dim.y+decal.y;
       for(int yi=0; yi<dim.y-decal.y; yi++) {
                for(int xi=0; xi<dim.x-decal.x; xi++)
                { //cout<<"xi,yi="<<xi<<","<<yi<<endl;
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
                        result[pixel_shift]=entree[pixel];
                }
                for(int xi=dim.x-decal.x; xi<dim.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
                        result[pixel_shift]=entree[pixel];
                }
        }
              for(int yi=dim.y-decal.y; yi<dim.y; yi++) {
                for(int xi=0; xi<dim.x-decal.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(-dim.y+yi+decal.y)*dim.x+xi+decal.x;
                        result[pixel_shift]=entree[pixel];
                }
                for(int xi=dim.x-decal.x; xi<dim.x; xi++)
                {
                        int pixel=yi*dim.x+xi;
                        int pixel_shift=(-dim.y+yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
                        result[pixel_shift]=entree[pixel];
                }
        }

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

                                volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2]=volume3D[nb_pixelz+nb_lignes+xi];
                        }
                }
        }
}
void circshift3DCplx(nbCplx *volume3D, nbCplx *volume3D_shift, Var3D dimFinal3D, Var3D decal3D)//(entrée, sortie_decalee, dim, decalage)
{
        decal3D.x=decal3D.x%dimFinal3D.x;//éliminer les "modulos"
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
                            volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2].Re=volume3D[nb_pixelz+nb_lignes+xi].Re;
                            volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2].Im=volume3D[nb_pixelz+nb_lignes+xi].Im;
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

void AXB_cplx(nbCplx *A, nbCplx *B, nbCplx *produit, int NbPix3D)
{

  for(int cpt=0;cpt<NbPix3D;cpt++)
        {
        produit[cpt].Re=A[cpt].Re*B[cpt].Re - A[cpt].Im*B[cpt].Im;
        produit[cpt].Im=B[cpt].Im*A[cpt].Re+A[cpt].Im*B[cpt].Re;
        }

}

void conj_cplx(nbCplx *A, nbCplx *conjA, int NbPix3D)
{

     for(int cpt=0;cpt<NbPix3D;cpt++)
        {
        conjA[cpt].Re=A[cpt].Re;
        conjA[cpt].Im=-A[cpt].Im;
        }

}
void recal_obj(nbCplx *a, nbCplx *b,nbCplx *objRecal, Var3D dimVol)
{
    int NPix3D=dimVol.x*dimVol.y*dimVol.z;
    nbCplx *A=new nbCplx[NPix3D];
    nbCplx *B=new nbCplx[NPix3D];



    TF3DCplx(a, A,dimVol);
    TF3DCplx(b, B,dimVol);


    nbCplx BConj[NPix3D];
    nbCplx prodAB[NPix3D];
    nbCplx decalPhi[NPix3D];
    conj_cplx(B, BConj, NPix3D);

    nbCplx produitABconj[NPix3D];

    AXB_cplx(A, BConj, prodAB, NPix3D);

for(int cpt=0;cpt<NPix3D;cpt++)
        {
            decalPhi[cpt].Re=prodAB[cpt].Re/sqrt(pow(A[cpt].Re,2))+sqrt(pow(BConj[cpt].Re,2));
            decalPhi[cpt].Im=prodAB[cpt].Im/sqrt(pow(A[cpt].Re,2))+sqrt(pow(BConj[cpt].Re,2));
        }

}

void multiplier_masque2Cplx(nbCplx *image, double masque[], int t_image, int t_mask, Var2D Centre)
{
        int t_imageX=t_image;
        int t_imageY=t_image;
        //if((t_mask%2)!=0)
        //cout<<"fonction multiplier_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;

        for(int x=Centre.x-t_mask/2; x<Centre.x+t_mask/2; x++) {
                for(int y=Centre.y-t_mask/2; y<Centre.y+t_mask/2; y++) {
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
                        //if((double(masque[t_mask*(y-Centre.y+t_mask/2)+(x-Centre.x+t_mask/2)])/255)!=1)
                        ///cout<<"masque:"<<(double(masque[t_mask*(y-Centre.y+t_mask/2)+(x-Centre.x+t_mask/2)])/255)<<endl;
                        image[cptj].Re=image[cptj].Re*masque[t_mask*(y-Centre.y+t_mask/2)+(x-Centre.x+t_mask/2)];
                        image[cptj].Im=image[cptj].Im*masque[t_mask*(y-Centre.y+t_mask/2)+(x-Centre.x+t_mask/2)];
                }
        }
}

///découpe une une fenetre de dimension dim_dest, coin haut gauche coin, dans src de taille dim_src.
void coupeCplx(nbCplx *src, nbCplx *dest, Var2D dim_src, Var2D dim_dest, Var2D coin)
{
    size_t nbPixSrc=dim_src.x*dim_src.y;
    size_t cpt_destX,cpt_destY, cpt_dest1D,
            cpt_srcX,cpt_srcY, cpt_src1D;

    for(cpt_destX=0; cpt_destX<dim_dest.x; cpt_destX++){
        for(cpt_destY=0; cpt_destY<dim_dest.y; cpt_destY++){

            cpt_dest1D=cpt_destX+cpt_destY*dim_dest.x;///coord 1D destination

            cpt_srcX=coin.x+cpt_destX;///coord X src
            cpt_srcY=coin.y+cpt_destY;///coord Y src
            cpt_src1D=cpt_srcX+cpt_srcY*dim_src.x;///coord 1D source

            dest[cpt_dest1D].Re=src[cpt_src1D].Re;
            dest[cpt_dest1D].Im=src[cpt_src1D].Im;

        }

    }

}


void methodeCarre(int NbPixROI2d, double *holo1,  double *holo2,  double *holo3,  double *holo4)
{
///Calcul du déphasage réel (méthode de Carré), et du taux de modulation

                        /// variable pour méthode de carré---------------------------
        double* holo1Re=new double[NbPixROI2d];
        double* holo2Re=new double[NbPixROI2d];
        double* holo3Re=new double[NbPixROI2d];
        double* holo4Re=new double[NbPixROI2d];
        double* txModulFrange=new double[NbPixROI2d];

        double* phaseMod2piIm=new double[NbPixROI2d];
        double* TfPhaseMod2piIm=new double[NbPixROI2d];
        double* TfPhaseMod2pi=new double[NbPixROI2d];
/*
        memset(TfPhaseMod2piIm,0,sizeof(TfPhaseMod2piIm));

        double txModulFrange=0;

        /*double* holo1Im=new double[NbPixROI2d];
        double* holo2Im=new double[NbPixROI2d];
        double* holo3Im=new double[NbPixROI2d];
        double* holo4Im=new double[NbPixROI2d];
        //partie imaginaire nulle en entrée
        memset(holo1Im,0,sizeof(holo1Im));
        memset(holo2Im,0,sizeof(holo2Im));
        memset(holo3Im,0,sizeof(holo3Im));
        memset(holo4Im,0,sizeof(holo4Im));*/
        double denominAlpha=0;
        double numeratAlpha=0;
        double argTanAlpha=0;//\delta=2*alpha, bref vaudrait mieux trouver 45° et tan alpha autour de 1

        double* dephasageCarre=new double[NbPixROI2d];

        double *phaseMod2pi=new double[NbPixROI2d];//calcul de la phase [2pi] de l'hologramme, méthode de carré
        memset(phaseMod2pi,0,sizeof(phaseMod2pi));
/// fin déclaration variable carré----------------------------
                                 double somDephasageCarre=0;
                                 int nb_tan=0, cptPhiAberrant=0;

                                 //calcul de la valeur du saut de phase et de la phase elle même
                                 for(int pixel=0;pixel<NbPixROI2d;pixel++)
                                 		{
                                 			holo1Re[pixel]=(double)holo1[pixel];
                                 			holo2Re[pixel]=(double)holo2[pixel];
                                 			holo3Re[pixel]=(double)holo3[pixel];
                                 			holo4Re[pixel]=(double)holo4[pixel];

                                            denominAlpha=((holo1Re[pixel]-holo4Re[pixel])+(holo2Re[pixel]-holo3Re[pixel]));
                                            if(denominAlpha!=0)
                                            {
                                                numeratAlpha=(3*(holo2Re[pixel]-holo3Re[pixel])-(holo1Re[pixel]-holo4Re[pixel]));
                                                argTanAlpha=numeratAlpha/denominAlpha;//formule de Carré
                                                if(argTanAlpha>0)
                                                {
                                                   // dephasageCarre[pixel]=2*atan(sqrt(argTanAlpha))*180/3.14;
                                                    dephasageCarre[pixel]=sqrt(argTanAlpha);
                                                }
                                                else
                                                {
                                                    cptPhiAberrant++;
                                                    dephasageCarre[pixel]=1;
                                                }
                                            }
                                            else
                                            {
                                                numeratAlpha=1;
                                                argTanAlpha=numeratAlpha/denominAlpha;//formule de Carré
                                                if(argTanAlpha>0)
                                                {
                                                  // dephasageCarre[pixel]=2*atan(sqrt(argTanAlpha))*180/3.14;
                                                   dephasageCarre[pixel]=sqrt(argTanAlpha);

                                                }
                                                else
                                                {
                                                    cptPhiAberrant++;
                                                    dephasageCarre[pixel]=1;
                                                }

                                            }//fin calcul Dphi
                                  			///calcul phase
                                  			double denominPhase=holo2Re[pixel]+holo3Re[pixel]-holo1Re[pixel]-holo4Re[pixel];
                                  			double argNumPhase=(3*(holo2Re[pixel]-holo3Re[pixel])-holo1Re[pixel]+holo4Re[pixel])*(holo1Re[pixel]-holo4Re[pixel]+holo2Re[pixel]-holo3Re[pixel]);
                                  			if(argNumPhase>0)
                                            {
                                                if (denominPhase==0)
                                                {
                                                    denominPhase=1;//solution la plus facile!
                                                    phaseMod2pi[pixel]=atan(sqrt(argNumPhase)/denominPhase);
                                                }
                                                else
                                                    {
                                                    phaseMod2pi[pixel]=phaseMod2pi[pixel]=atan(sqrt(argNumPhase)/denominPhase);
                                                     //  phaseMod2pi[pixel]=0;//sqrt((3*(holo2Re[pixel]-holo3Re[pixel])-holo1Re[pixel]+holo4Re[pixel])*(holo1Re[pixel]-holo4Re[pixel]+holo2Re[pixel]-holo3Re[pixel]))/denominPhase;
                                                   }
                                            }
                                            else
                                            {
                                                phaseMod2pi[pixel]=0;
                                            }
                                 			if(dephasageCarre[pixel]<1.4&&dephasageCarre[pixel]>0.6)//la valeur de la tangente doit valoir environ tan(0.5*pi/2)=1, on élimine les valeurs aberrantes>10
                                 			{
                                                somDephasageCarre=dephasageCarre[pixel]+somDephasageCarre;
                                                txModulFrange[pixel]=2*sqrt((pow(holo4Re[pixel]-holo2Re[pixel],2)+pow(holo1Re[pixel]-holo3Re[pixel],2)))/(holo1Re[pixel]+holo2Re[pixel]+holo3Re[pixel]+holo4Re[pixel]);
                                                nb_tan++;
                                 			}
                                 		}
                        //cout<<"Singularités de Delta phi : "<<cptPhiAberrant<<"%"<<endl;
                                        ///démodulation
                                       // phaseMod2pi=circshift(phaseMod2pi, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);
                                          /* holo1Re=circshift(holo1Re, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);
                                            holo2Re=circshift(holo2Re, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);
                                            holo3Re=circshift(holo3Re, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);
                                            holo4Re=circshift(holo4Re, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);*/
                        //TF2D(phaseMod2pi,phaseMod2piIm,TfPhaseMod2pi,TfPhaseMod2piIm,DIMX_CCD2,DIMY_CCD2);
                        //TfPhaseMod2pi=circshift(TfPhaseMod2pi, DIMX_CCD2,DIMY_CCD2,DIMX_CCD2/2,DIMY_CCD2/2);
///SAV tf phase
//SAV(TfPhaseMod2pi, NbPixROI2d, "/home/mat/tomo_test/IDP/TfPhaseMod2pi/TfPhaseMod2pi", CHAR,"a+b");
///SAV delta
//SAV(dephasageCarre, NbPixROI2d, "/home/mat/tomo_test/IDP/deltaCarre/deltaCarre", CHAR,"a+b");
///sav phase
//SAV(phaseMod2pi, NbPixROI2d, "/home/mat/tomo_test/IDP/phase/phase", CHAR,"a+b");
///sav txModulation
//SAV(txModulFrange, NbPixROI2d, "/home/mat/tomo_test/IDP/txModulation/txModulation", CHAR,"a+b");
///-------
                        // double moyenne_tan_alpha=somme_tan_alpha/N;
                        /*  for(int pixel=0;pixel<N;pixel++)
                          		{
                           			if(tan_alpha[pixel]<1.4&&tan_alpha[pixel]>0.6)//la valeur de la tangente doit valoir environ tan(0.5*pi/2)=1, on élimine les valeurs aberrantes>10
                          			{
                                     somme_ecart_carre=somme_ecart_carre+(tan_alpha[pixel]-moyenne_tan_alpha)*(tan_alpha[pixel]-moyenne_tan_alpha);//         calcul écart type de la valeur du saut de phase
                          			}
                          		}
                         double ecart_type=sqrt(somme_ecart_carre/N);*/
                        //affichage saut de phase
                        // printf("la moyenne du phase shifting sur cette image vaut: %f avec nb_tan :%f et ecart type : %f\n",moyenne_tan_alpha,sqrt(nb_tan),ecart_type);
                        //allocation résultat du phase shifting attention à ne pas désallouer hors de la boucle for
                        //cout<<"Décalage phase estimée : "<<atan(moyenne_tan_alpha)*180/3.14*2<<"°"<<endl;
                        //	cout<<"Taux moyen de modulation des Franges"<<txModulFrange<<endl;
                        //Libérer la mémoire des 4 images expérimentales+methode carré
                        ///Libération mémoire carré.
        delete[] holo1Re, holo2Re, holo3Re, holo4Re;
        delete[] holo1, holo2, holo3, holo4;
        delete[] phaseMod2piIm;
        delete[] TfPhaseMod2piIm;
        delete[] TfPhaseMod2pi;
        delete[] txModulFrange;
        /*delete[] holo1Im, holo2Im, holo3Im, holo4Im;
        delete[] TfHolo1Re, TfHolo2Re, TfHolo3Re, TfHolo4Re;
         delete[] TfHolo1Im, TfHolo2Im, TfHolo3Im, TfHolo4Im;*/

        delete[] dephasageCarre;
        delete[] phaseMod2pi;
///Fin méthode de carré
}
int CreerZoneFresnel(double *FresnelRe,double * FresnelIm, Var2D dim, Var2D centre,float d, float lambda)
{
        for(int cpty=0;cpty<dim.y;cpty++)
        {
            int hauteur=cpty*dim.x;
            for(int cptx=0;cptx<dim.x;cptx++)
            {
                int cpt=hauteur+cptx;
                FresnelRe[cpt]=cos(3.14159*(pow(cptx-centre.x,2)+pow(cpty-centre.y,2))/(lambda*d));
                FresnelIm[cpt]=sin(3.14159*(pow(cptx-centre.x,2)+pow(cpty-centre.y,2))/(lambda*d));
            }
}
    return 1;
}

///#########lecture propre d'un fichire binaire, connaissant sa taille et son type de données
    void lire_bin(string chemin, double resultat[], short int precision, Var3D dim_entree)
   {
       unsigned long lTaille;
       size_t nb_elmnt_lu;
       unsigned short int dimOctet=precision/8;//taille en octet d'un element.double=4 octet
       FILE* pFichier = NULL;
       pFichier = fopen(chemin.c_str(), "r");  //ouverture de ce fichier en écriture binaire

       if(pFichier==NULL)
       {
           fputs("Impossible d'ouvrir le fichier\n",stderr);
           exit (1);// obtenir la longueur du fichier, comparer avec donnée entrée.
           fseek(pFichier,0,SEEK_END);//trouver la fin de fichier
           lTaille = ftell (pFichier);//position de la fin du fichier
           printf("taille trouvé %li, taille estimée : %i\n",lTaille,dim_entree.x*dim_entree.y*dim_entree.z);
           rewind(pFichier);

           if(dimOctet*dim_entree.x*dim_entree.y*dim_entree.z!=lTaille)
            cout<<"Taille du fichier incompatible avec les dimensions\n"<<endl;

            nb_elmnt_lu = fread (resultat,1,lTaille,pFichier);//lecture

            if(nb_elmnt_lu!=lTaille)
                cout<<"Problème lors de la lecture du fichier"<<endl;

            fclose(pFichier);
        }

    }
#include <tiffio.h>

///Ã  ajouter dans fonctions.cpp

void SAV_Tiff2D(double *var_sav, char *chemin, int dim)
{
    TIFF *tif= TIFFOpen(chemin, "a");
    TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, dim);
    TIFFSetField (tif, TIFFTAG_IMAGELENGTH, dim);
    TIFFSetField (tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField (tif, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField (tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField (tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField (tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

    tsize_t strip_size = TIFFStripSize (tif);
    tstrip_t strips_num = TIFFNumberOfStrips (tif);

    float* strip_buf=(float*)_TIFFmalloc(strip_size);
    for (unsigned int s=0; s<strips_num; s++)
    {
        for (unsigned int col=0; col<dim; col++)
        {
            unsigned int cpt=col+dim*s;
            strip_buf[col]=(float)var_sav[cpt];
        }
        TIFFWriteEncodedStrip (tif, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif);
    TIFFClose(tif);
}

void SAV_Tiff_Re(nbCplx *var_sav, char *chemin, int dim)
{
    TIFF *tif= TIFFOpen(chemin, "a");
    TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, dim);
    TIFFSetField (tif, TIFFTAG_IMAGELENGTH, dim);
    TIFFSetField (tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField (tif, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField (tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField (tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField (tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

    tsize_t strip_size = TIFFStripSize (tif);
    tstrip_t strips_num = TIFFNumberOfStrips (tif);

    float* strip_buf=(float*)_TIFFmalloc(strip_size);
    for (unsigned int s=0; s<strips_num; s++)
    {
        for (unsigned int col=0; col<dim; col++)
        {
            unsigned int cpt=col+dim*s;
            strip_buf[col]=(float)var_sav[cpt].Re;
        }
        TIFFWriteEncodedStrip (tif, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif);
    TIFFClose(tif);
}

void SAV_Tiff_Im(nbCplx *var_sav, char *chemin, int dim)
{
    TIFF *tif= TIFFOpen(chemin, "a");
    TIFFSetField (tif, TIFFTAG_IMAGEWIDTH, dim);
    TIFFSetField (tif, TIFFTAG_IMAGELENGTH, dim);
    TIFFSetField (tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField (tif, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField (tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField (tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField (tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField (tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    TIFFSetField (tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

    tsize_t strip_size = TIFFStripSize (tif);
    tstrip_t strips_num = TIFFNumberOfStrips (tif);

    float* strip_buf=(float*)_TIFFmalloc(strip_size);
    for (unsigned int s=0; s<strips_num; s++)
    {
        for (unsigned int col=0; col<dim; col++)
        {
            unsigned int cpt=col+dim*s;
            strip_buf[col]=(float)var_sav[cpt].Im;
        }
        TIFFWriteEncodedStrip (tif, s, strip_buf, strip_size);
    }
    _TIFFfree(strip_buf);
    TIFFWriteDirectory(tif);
    TIFFClose(tif);
}

void SAV_Tiff3D(nbCplx *var_sav, char *chemin_ind, char *chemin_abs, int dim)
{
    double *bufferRe=new double[dim*dim];
    double *bufferIm=new double[dim*dim];
    for(int z=0; z < dim; z++)
    {
        for(int y=0; y < dim; y++)
        {
            for(int x=0; x < dim; x++)
            {
                int cpt1= x + y * dim + z * dim * dim;
                int cpt2= x + y * dim;
                bufferRe[cpt2]=var_sav[cpt1].Re;
                bufferIm[cpt2]=var_sav[cpt1].Im;
            }
        }
        SAV_Tiff2D(bufferRe, chemin_ind, dim);
        SAV_Tiff2D(bufferIm, chemin_abs, dim);
	    }
    delete[] bufferRe, bufferIm;
}
