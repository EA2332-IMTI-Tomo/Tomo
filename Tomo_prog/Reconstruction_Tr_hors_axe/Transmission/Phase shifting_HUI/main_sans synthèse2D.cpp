#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//#include <strstream>//obsolète!
#ifndef DEF_FONCTIONS// Si la constante n'a pas été définie` le fichier n'a jamais été inclus
#define DEF_FONCTIONS
#include "fonction.hpp"
#endif

int main(int argc, char *argv[])
{
int DIMX_CCD2=1016,  DIMY_CCD2=1016;
Var2D dimCCD={DIMX_CCD2,DIMY_CCD2};
//valeur controlant l'exclusion de certains centres
int xm0_limite=200; //centre maximum
int ym0_limite=200;
int rayon_inf=xm0_limite*xm0_limite+ym0_limite*ym0_limite;
//valeur du coin pour la découpe
//Var2D coin={440,220};//indentation
Var2D coin={150,0};
//if(coin.x+DIMX_CCD2>896 || coin.y+DIMY_CCD2>512)//largeur : 740, hauteur : 574.
if(coin.x+DIMX_CCD2>1280 || coin.y+DIMY_CCD2>1082)//largeur : 740, hauteur : 574.
{
    printf("Attention, la zone decoupee sort de l'image\n");
    return 0;
}

int jumeau_elimine=0;
int points_faux=0;
int fftwThreadInit;
fftwThreadInit=fftw_init_threads();
int nthreads=6;
printf("fftwThreadInit: %i\n",fftwThreadInit);


///###################################repertoire reflexion#########################################
//char Dossier_acquiz[]="/home/hui/maniptomo/reflexion_acquisition/2012031301/";
//char Dossier_acquiz[]="/home/hui/maniptomo/reflexion_acquisition/2012070502/IMAGE/";
char Dossier_acquiz[]="/home/hui/maniptomo/reflexion_acquisition/2012121004/";
char NumSession[]="i";
char NomHoloPart1[strlen(Dossier_acquiz)+strlen("i")]; //au max 9999 hologramme
concatener(Dossier_acquiz,NumSession,NomHoloPart1);
char Fichier_holo[strlen(NomHoloPart1)+strlen("XXXX_p0.pgm")+1];//Nom final

/*
//chemin pour référence d'amplitude
char Dossier_ref[]="/home/hui/maniptomo/reflexion_acquisition/2012120403/";
char CheminModRef[strlen(Dossier_ref)+strlen("mod_ref.bmp")+1];
concatener(Dossier_ref,"mod_ref.bmp",CheminModRef);
*/
/*
///############################repertoire transmission###########################################
char Dossier_acquiz[]="/home/hui/maniptomo/image_acquisition/session12060801/";
char NumSession[]="session12060801-";
char NomHoloPart1[strlen(Dossier_acquiz)+strlen("sessionXXXXXXXX")]; //au max 9999 hologramme
concatener(Dossier_acquiz,NumSession,NomHoloPart1);
char Fichier_holo[strlen(NomHoloPart1)+strlen("recordXXXX-004.bmp")+1];//Nom final
/*
//chemin pour référence d'amplitude
char Dossier_ref[]="/home/hui/maniptomo/image_reference/";
char CheminModRef[strlen(Dossier_ref)+strlen("mod_ref_XXXXXXXX.bmp")+1];
concatener(Dossier_ref,"mod_ref_18112011.bmp",CheminModRef);
*/

const int
premier_plan=1,
Num_Angle_final=479,
// int sans_signal=0;
NbAngle=Num_Angle_final-premier_plan+1,
SautAngle=1,
N=DIMX_CCD2*DIMY_CCD2;//nombre de pixel (pour le tableau 2D)
const float
n1=1.515,	//indice de l'huile
NA=1.40,	//ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
theta=asin(NA/n1),//theta_max
lambda=488*pow(10,-9), //longueur d'onde
G=100,	//grossissement telan+objectif
//Tp=11.2*pow(10,-6), //Taille des pixels
Tp=8*pow(10,-6),
Rf=1,//1.2;=1/facteur grossissement
K=lambda*G/(2*NA*Tp*Rf),// Facteur d'echelle
Tps=Tp,//*K; //Taille pixel apres mise à l'echelle important si RF<>1
tailleTheoPixelHolo=Tp/(G/Rf)*pow(10,9);//Pour info
float
rayon=n1*Rf*Tps*DIMX_CCD2/(G*lambda), //Rayon
delta_zmax=rayon*cos(theta);//float?
const int
NXMAX=round(n1*Rf*Tps*DIMX_CCD2/(G*lambda)*NA/n1);//NXMAX=round(rayon*NA/n1)=round(rayon*sin theta_max);
cout<<"NXMAX="<<NXMAX<<endl;
int
NYMAX=NXMAX,
N_tab=64*NXMAX*NXMAX*NXMAX,
centres_exclus=0,
nb_proj=0;
Var2D dimDecoup={NXMAX,NXMAX};

const int
dv0=round(4*NXMAX),//dimmension totale en x,y
dv0rf=round(4*NXMAX),
dv0s2rf=round(2*NXMAX),
dv1s2rf=round(2*NXMAX),
dv1xdv0rf=round(16*NXMAX*NXMAX),
dv2s2rf=round(2*NXMAX);

///reservation de 4 tableaux (image 2D)  pour le phase shifting
unsigned char* holo1=new unsigned char[N];
unsigned char* holo2=new unsigned char[N];
unsigned char* holo3=new unsigned char[N];
unsigned char* holo4=new unsigned char[N];

///variable pour image de référence
double* ampli_ref=new double[N]; //reservation reference et "noir" camera pour une "tentative" de correction de la ref
unsigned char* mod_ref=new unsigned char[N]; //reservation reference et "noir" camera pour une "tentative" de correction de la ref

///variable pour synthèse ouverture 2D
//on agrandit le support en doublant largeur et longueur (facteur 4)
Var2D dimInit={2*NXMAX,2*NYMAX},dimFin={2*NXMAX,2*NYMAX};

///variable pour TF
double *fft_reel_tmp=new double[N];
double *fft_imag_tmp=new double[N];
double *fft_reel=new double[4*NXMAX*NYMAX];
double *fft_imag=new double[4*NXMAX*NYMAX];
double *fft_reel_shift=new double[4*NXMAX*NYMAX];
double *fft_imag_shift=new double[4*NXMAX*NYMAX];
double *fft_module_shift=new double[4*NXMAX*NYMAX];
double *fft_reel_shift_norm=new double[4*NXMAX*NYMAX];
double *fft_imag_shift_norm=new double[4*NXMAX*NYMAX];
double  *spectre_reel=new double[dimFin.x*dimFin.y];
double  *spectre_imag=new double[dimFin.x*dimFin.y];

///variable pour plan
double *plan_reel=new double[N];
double *plan_imag=new double[N];
double *plan_reel_shift=new double[N];
double *plan_imag_shift=new double[N];
double *plan_reel_norm_shift=new double[dimFin.x*dimFin.y];
double *plan_imag_norm_shift=new double[dimFin.x*dimFin.y];
double *plan_reel_norm=new double[dimFin.x*dimFin.y];
double *plan_imag_norm=new double[dimFin.x*dimFin.y];
double *plan_module_norm=new double[dimFin.x*dimFin.y];

///variable pour phase
double *phaseMod1pi=new double[dimFin.x*dimFin.y];
double *centre=new double[dimFin.x*dimFin.y];
double *distCal = new double[dimFin.x*dimFin.y];

///variables qualité holo
double *qualiteHolo=new double[N];
double *sommeP_reel=new double[N];
double *sommeP_imag=new double[N];
double *plan_I42=new double[N];
double *plan_I13=new double[N];
double *supRedonQ=new double[N];

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
//rempli_tableau(mod_ref,CheminModRef,coin,dimCCD);

for(int cpt=0;cpt<N;cpt++)
{
    if (mod_ref[cpt]==0)
    {
      //cout<<"erreur: le minimum de l'intensité de image réf est à 0"<<endl;
      mod_ref[cpt]=mod_ref[cpt]+1;
    }
    else
    mod_ref[cpt]=mod_ref[cpt];

}

float alpha=0.1;//coeff pour le masque de tuckey
double* masque=new double[DIMX_CCD2*DIMY_CCD2];
tukey2D(masque,DIMX_CCD2,DIMY_CCD2,alpha);

///######################variable pour algo de déroulement de phase########################################
int image_height=2*NXMAX, image_width=2*NXMAX;
int image_size = image_height * image_width;
int No_of_Edges = (image_width)*(image_height-1) + (image_width-1)*(image_height);
double *phaseunwrap=new double[image_size];
PIXEL *Pixel = new PIXEL[image_size];
EDGE *Edge = new EDGE[No_of_Edges];

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

///########################################################debut boucle#############################################################
double coef_IDP=0;
for(int cpt_angle=premier_plan;cpt_angle<premier_plan+NbAngle;cpt_angle=cpt_angle+SautAngle)//boucle sur tous les angles
{
    //cout<<"cpt_angle : "<<cpt_angle<<endl;
	if((cpt_angle-100*(cpt_angle/100))==0)
	printf("cpt_angle=%i\n",cpt_angle);

///#################### reflexion############################
    char TamponChemin[sizeof(Fichier_holo)];
	char Num_angle[4];
	sprintf(Num_angle,"%04i",cpt_angle);
	//créer le chemin final
	concatener(NomHoloPart1,Num_angle,Fichier_holo);
    strcpy(TamponChemin,Fichier_holo);
	concatener(TamponChemin,"_p0.pgm",TamponChemin);
	test_existence = fopen(TamponChemin, "rb");

/*
    ///###################transmission#######################
    char TamponChemin[sizeof(Fichier_holo)];
	char Num_angle[4];
	sprintf(Num_angle,"record%i",cpt_angle);
	//créer le chemin final
	concatener(NomHoloPart1,Num_angle,Fichier_holo);
    strcpy(TamponChemin,Fichier_holo);
	concatener(TamponChemin,"-001.bmp",TamponChemin);
	test_existence = fopen(TamponChemin, "rb");
*/
    if(test_existence!=NULL)
    {
		fclose(test_existence);
		nb_proj++;

		holo1=charger_image2D(holo1,0,Fichier_holo, coin.x, coin.y,DIMX_CCD2,DIMY_CCD2);
		holo2=charger_image2D(holo2,1,Fichier_holo, coin.x, coin.y,DIMX_CCD2,DIMY_CCD2);
		holo3=charger_image2D(holo3,2,Fichier_holo, coin.x, coin.y,DIMX_CCD2,DIMY_CCD2);
		holo4=charger_image2D(holo4,3,Fichier_holo, coin.x, coin.y,DIMX_CCD2,DIMY_CCD2);
//
/////######################################################décalage de phase##########################################################
//
//		for(int pixel=0;pixel<N;pixel++)
//		{
//            plan_reel[pixel]=((double)holo1[pixel]-(double)holo3[pixel])*(masque[pixel]);
//            plan_imag[pixel]=((double)holo4[pixel]-(double)holo2[pixel])*(masque[pixel]);
//		}
//        SAV(plan_reel, N, "/home/hui/maniptomo/IDP/champ_2d/plan_reel.bin", FLOAT,"a+b");

///#########################################décalage de phase en supprimer les faibles SRB##########################################
        memset(plan_I42, 0, N*8);
        memset(plan_I13, 0, N*8);
		for(int pixel=0;pixel<N;pixel++)
		{
            plan_I13[pixel]=(double)holo1[pixel]-(double)holo3[pixel];
            plan_I42[pixel]=(double)holo4[pixel]-(double)holo2[pixel];
		}

        memset(qualiteHolo, 0, N*8);
        memset(sommeP_reel, 0, N*8);
        memset(sommeP_imag, 0, N*8);
        int pixelQ_max=0;
        double sommeQholo=0;
        double moyQholo=0;
		for(int pixel=0;pixel<N;pixel++)
		{
            ///qualité des hologrammes
            qualiteHolo[pixel]=2*sqrt((pow(holo4[pixel]-holo2[pixel],2)+pow(holo1[pixel]-holo3[pixel],2)))/(holo1[pixel]+holo2[pixel]+holo3[pixel]+holo4[pixel]);
            //cout<<"contraste pour chaque pixel: "<<qualiteHolo[pixel]<<endl;
            /*if(qualiteHolo[pixel]>qualiteHolo[pixelQ_max])
            {
                pixelQ_max=pixel;
            }*/
            if(qualiteHolo[pixel]>0)
            {
                sommeQholo=sommeQholo+qualiteHolo[pixel];
            }

            if(qualiteHolo[pixel]>0.2)
            {
                sommeP_reel[pixel]=sommeP_reel[pixel]+plan_I13[pixel];
                sommeP_imag[pixel]=sommeP_imag[pixel]+plan_I42[pixel];
                if(supRedonQ[pixel]!=0)
                supRedonQ[pixel]+=1;
                //cout<<"supRedonQ="<<supRedonQ[pixel]<<endl;
            }
        }
        //cout<<"pixelQ_max="<<pixelQ_max<<endl;
        //cout<<"maxQholo="<<qualiteHolo[pixelQ_max]<<endl;
        moyQholo=sommeQholo/N;
        cout<<"Moyen de qualite holo="<<moyQholo<<endl;

        SAV(qualiteHolo, N, "/home/hui/maniptomo/IDP/champ_2d/qualiteHolo.bin", FLOAT,"a+b");
        SAV(supRedonQ, N, "/home/hui/maniptomo/IDP/champ_2d/supRedonQ.bin", FLOAT,"a+b");

        //if(qualiteHolo[pixelQ_max]>0.8)
        if(moyQholo>0.32)
        {

        for(int pixel=0;pixel<N;pixel++)
        {
            if(supRedonQ[pixel]==0)
            supRedonQ[pixel]=1;
            plan_reel[pixel]=sommeP_reel[pixel]*(masque[pixel])/supRedonQ[pixel];
            plan_imag[pixel]=sommeP_imag[pixel]*(masque[pixel])/supRedonQ[pixel];

        }
        SAV(plan_reel, N, "/home/hui/maniptomo/IDP/champ_2d/plan_reel.bin", FLOAT,"a+b");

///######################################################Circshift avant TF2D#######################################################
        Var2D dimCCD={DIMX_CCD2,DIMY_CCD2},decalCCD={DIMX_CCD2/2,DIMY_CCD2/2};
		circshift2(plan_reel,plan_reel_shift, dimCCD,decalCCD);
		circshift2(plan_imag,plan_imag_shift, dimCCD,decalCCD);
        //SAV(plan_reel_shift, N, "/home/hui/maniptomo/IDP/champ_2d/plan_reel_shift.bin", FLOAT,"a+b");

///########################################TF 2D du front d'onde: passage du plan image au plan réciproque##########################
		TF2D(plan_reel_shift,plan_imag_shift,fft_reel_tmp,fft_imag_tmp,DIMX_CCD2,DIMY_CCD2);
        //SAV(fft_reel_tmp, N, "/home/hui/maniptomo/IDP/champ_2d/TF2D_reel.bin", FLOAT,"a+b");

///######################################################Découpage dans fourier à NXMAX#############################################
        decalCoupe(fft_reel,fft_imag,fft_reel_tmp,fft_imag_tmp,dimDecoup,dimCCD);
        //SAV(fft_reel, 4*NXMAX*NYMAX, "/home/hui/maniptomo/IDP/champ_2d/TF2D_reel_Cut.bin", FLOAT,"a+b");
		//printf("  NXMAX : %i, NYMAX : %i\n",NXMAX,NYMAX);

///######################################################circshift 2D apres TF######################################################
        Var2D dim2D={2*NXMAX,2*NYMAX},decal2D={NXMAX,NYMAX};
        circshift2(fft_reel,fft_reel_shift, dim2D,decal2D);
	    circshift2(fft_imag,fft_imag_shift, dim2D,decal2D);

///######################################################calcul maximum de spéculaire###############################################
		int cpt_max=0;
		fft_module_shift[0]=pow(fft_reel_shift[0],2)+pow(fft_imag_shift[0],2);

		//Recherche du MAX dansle module
		for(int cpt=1;cpt<(4*NXMAX*NYMAX);cpt++)
			{
				fft_module_shift[cpt]=pow(fft_reel_shift[cpt],2)+pow(fft_imag_shift[cpt],2);
				if(fft_module_shift[cpt]>fft_module_shift[cpt_max])
				{
					cpt_max=cpt;
				}
			}
		double max_part_reel = fft_reel_shift[cpt_max];
		double max_part_imag = fft_imag_shift[cpt_max];
		double max_module = fft_module_shift[cpt_max];
		//printf("max_part_reel : %f, max_part_imag: %f,  max_module %f \n ",max_part_reel, max_part_imag, max_module);

///#######Coordonnées de l'objet et jumeaux dans l'espace 2D à partir de l'indice 1D: (xc,yc)=(0,0)=en haut à gauche de l'image#####
		int xmi=cpt_max%(2*NXMAX);
		int ymi=cpt_max/(2*NYMAX);
		int xc=xmi-NXMAX;
		int yc=ymi-NYMAX;
		int xmij=2*NXMAX-xmi;//coordonnée objet jumeaux
		int ymij=2*NYMAX-ymi;
        //cout<<"Coordonnées jumeau angle "<<cpt_angle<<endl<<"xmij="<<xmij<<endl<<"ymij="<<ymij<<endl;
        int cptJumo_max=2*NYMAX*ymij+xmij;

///######################################################Virer l'objet jumeaux+ordre zero###########################################
		if((xc*xc+yc*yc)>9)//35*35=1225;objet jumeau pas trop près de l'objet
		{ jumeau_elimine++;
            int t_mask_var=round(sqrt(xc*xc+yc*yc)/7)*2+1;
           // cout<<"tmask="<<t_mask_var<<endl;
            double* cache_jumeau3=new double[t_mask_var*t_mask_var];
            antigaussienne(cache_jumeau3,t_mask_var,round(t_mask_var*2),1,0);
            multiplier_masque2(fft_imag_shift, cache_jumeau3, NXMAX, t_mask_var, xmij, ymij);//jumeau
            multiplier_masque2(fft_reel_shift, cache_jumeau3, NXMAX, t_mask_var, xmij, ymij);
            multiplier_masque2(fft_reel_shift, cache_jumeau3, NXMAX, t_mask_var, NXMAX, NXMAX);//ordre zéro
            multiplier_masque2(fft_imag_shift, cache_jumeau3, NXMAX, t_mask_var, NXMAX, NXMAX);
            delete[] cache_jumeau3;
		}//fin if objet jumeau pas trop pres

///######################################################normalisation par le pic central###########################################
		for(int cpt=0;cpt<(4*NXMAX*NYMAX);cpt++)
			{
                fft_reel_shift_norm[cpt]=(fft_reel_shift[cpt]*max_part_reel+fft_imag_shift[cpt]*max_part_imag)/max_module;
                fft_imag_shift_norm[cpt]=(fft_imag_shift[cpt]*max_part_reel-fft_reel_shift[cpt]*max_part_imag)/max_module;
                //printf("fft_reel_shift_norm: %f, fft_imag_shift_norm: %f \n ",fft_reel_shift_norm[cpt], fft_imag_shift_norm[cpt]);
			}
        SAV(fft_reel_shift_norm, 4*NXMAX*NYMAX, "/home/hui/maniptomo/IDP/champ_2d/fft_reel_shift_norm.bin", FLOAT,"a+b");

///on redécale le spectre pour éliminer les franges d inclinaisons (-xc), on en profite pour revenir au zero informatique (+2*NXMAX)
        Var2D decal{-xc+NXMAX,-yc+NXMAX};
        circshift2(fft_reel_shift_norm,spectre_reel,dimFin,decal);
        circshift2(fft_imag_shift_norm,spectre_imag,dimFin,decal);

///#####################################exportation spectre_2D######################################################################
        SAV(spectre_reel, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/spectre2D.bin", FLOAT,"a+b");
        TF2D_INV(spectre_reel,spectre_imag,plan_reel_norm_shift,plan_imag_norm_shift,dimFin.x,dimFin.y);
        Var2D fftDecal{NXMAX,NXMAX};
        circshift2(plan_reel_norm_shift,plan_reel_norm,dimFin,fftDecal);
        circshift2(plan_imag_norm_shift,plan_imag_norm,dimFin,fftDecal);

///#####################################calcul module de champ2D####################################################################
        for(int cpt=0;cpt<(dimFin.x*dimFin.y);cpt++)
        {
            plan_module_norm[cpt]=pow(plan_reel_norm[cpt],2)+pow(plan_imag_norm[cpt],2);
        }

///#####################################déroulement de phase sur tous les angles####################################################
        phaseUnwrapping(plan_reel_norm, plan_imag_norm, phaseMod1pi, phaseunwrap, Pixel, Edge, dimFin.x);
        SAV(phaseMod1pi, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/phaseMod1pi.bin", FLOAT,"a+b");
        write_data("/home/hui/maniptomo/IDP/champ_2d/phaseUnwrapped.bin",phaseunwrap,image_size);

        for(int cpt=0; cpt<dimFin.x*dimFin.y; cpt++)
        {
            distCal[cpt] = phaseunwrap[cpt]*lambda*pow(10,9)/(4*PI*n1);
            //cout<<"d="<<d<<endl;
        }
    ///exportation distance
        SAV(distCal, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/distanceTAngles.bin", FLOAT,"a+b");

    ///exportation champ_2D
        SAV(plan_reel_norm, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/champ2D.bin", FLOAT,"a+b");

    temps_final = clock ();
	temps_cpu = (temps_final - temps_initial) * 1e-6;
	//printf("temps apres lecture : %f\n",temps_cpu);
    temps_total=temps_total+temps_cpu;
    temps_initial = clock();

} //fin de test de validite du nom de fichier
    else{printf("qualité d'hologramme %i insuffisante\n",cpt_angle);
    }
}
else{
	printf("fichier %i inexistant\n",cpt_angle);
	//fclose(test_existence);
	}
}//fin de boucle for sur tous les angles  on peut désallouer les variables définit hors de la boucle

delete[] centre, ampli_ref, mod_ref, supRedonQ;
delete[] fft_reel, fft_imag, fft_reel_tmp, fft_imag_tmp;
delete[] sommeP_imag, sommeP_reel;
delete[] qualiteHolo;
delete[] plan_I13, plan_I42;
delete[] masque;
delete[] holo1, holo2, holo3, holo4;
delete[] fft_reel_shift_norm, fft_imag_shift_norm;
delete[] plan_reel, plan_imag, plan_reel_norm, plan_imag_norm, plan_module_norm, plan_reel_norm_shift, plan_imag_norm_shift;
delete[] fft_reel_shift, fft_imag_shift, fft_module_shift, plan_reel_shift, plan_imag_shift, spectre_reel, spectre_imag;
delete[] distCal, phaseMod1pi, phaseunwrap, Edge, Pixel;
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

//libération memoire allouée pour les threads
void fftw_cleanup_threads(void);
return 0;
}


