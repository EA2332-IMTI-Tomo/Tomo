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
int DIMX_CCD2=356,DIMY_CCD2=356;
Var2D dimCCD={DIMX_CCD2,DIMY_CCD2};
const int N=DIMX_CCD2*DIMY_CCD2;//nombre de pixel (pour le tableau 2D)
//valeur controlant l'exclusion de certains centres
int xm0_limite=200; //centre maximum
int ym0_limite=200;
int rayon_inf=xm0_limite*xm0_limite+ym0_limite*ym0_limite;
//valeur du coin pour la découpe
Var2D coin={460,460};
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
char Dossier_acquiz[]="/home/hui/maniptomo/reflexion_acquisition/2013060401_VLSI_def/";
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
char Dossier_acquiz[]="/home/hui/maniptomo/image_acquisition/session13022601/syn2d/";
char NumSession[]="session13022601-";
char NomHoloPart1[strlen(Dossier_acquiz)+strlen("sessionXXXXXXXX")]; //au max 9999 hologramme
concatener(Dossier_acquiz,NumSession,NomHoloPart1);
char Fichier_holo[strlen(NomHoloPart1)+strlen("recordXXXX-004.bmp")+1];//Nom final

//chemin pour référence d'amplitude
char Dossier_ref[]="/home/hui/maniptomo/image_acquisition/session13020601/";
char CheminModRef[strlen(Dossier_ref)+strlen("mod_ref.bmp")+1];
concatener(Dossier_ref,"mod_ref.bmp",CheminModRef);
*/
///variable pour image de référence
double* ampli_ref=new double[N]; //reservation reference et "noir" camera pour une "tentative" de correction de la ref
unsigned char* mod_ref=new unsigned char[N]; //reservation reference et "noir" camera pour une "tentative" de correction de la ref

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


const int
premier_plan=0,
Num_Angle_final=99,
// int sans_signal=0;
NbAngle=Num_Angle_final-premier_plan+1,
SautAngle=1;

const float
n1=1.33,	//indice de l'huile
NA=1.2,	//ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
theta=asin(NA/n1),//theta_max
lambda=475*pow(10,-9), //longueur d'onde
G=60,	//grossissement telan+objectif
//Tp=11.2*pow(10,-6), //Taille des pixels
Tp=8*pow(10,-6),
Rf=0.5,//1.2;=1/facteur grossissement
K=lambda*G/(2*NA*Tp*Rf),// Facteur d'echelle
Tps=Tp,//*K; //Taille pixel apres mise à l'echelle important si RF<>1
tailleTheoPixelHolo=Tp/(G/Rf)*pow(10,9);//Pour info
float
rayon=n1*Rf*Tps*DIMX_CCD2/(G*lambda), //Rayon
delta_zmax=rayon*cos(theta);//float?
cout<<"rayon="<<rayon<<endl;
cout<<"taillePTheo="<<tailleTheoPixelHolo<<endl;

const int
NXMAX=round(n1*Rf*Tps*DIMX_CCD2/(G*lambda)*NA/n1);//NXMAX=round(rayon*NA/n1)=round(rayon*sin theta_max);
cout<<"NXMAX="<<NXMAX<<endl;

float tailleTheoPixelTomo=tailleTheoPixelHolo*DIMX_CCD2/(4*NXMAX);
printf("Taille théorique des pixels Tomo: %f\n",tailleTheoPixelTomo);

int
NYMAX=NXMAX,
N_tab=64*NXMAX*NXMAX*NXMAX,
centres_exclus=0,
nb_proj=0;
Var2D dimDecoup={NXMAX,NXMAX};

float obj,
coef,
coefSomme=0;

int nomfinal=0;
Var2D posDecalAmpli={0,0};
Var2D posDecalPhase={0,0};

///reservation de 4 tableaux (image 2D)  pour le phase shifting
unsigned char* holo1=new unsigned char[N];
unsigned char* holo2=new unsigned char[N];
unsigned char* holo3=new unsigned char[N];
unsigned char* holo4=new unsigned char[N];

///variable pour synthèse ouverture 2D
//on agrandit le support en doublant largeur et longueur (facteur 4)
Var2D dimInit={2*NXMAX,2*NYMAX},dimFin={4*NXMAX,4*NYMAX};

///variable pour phase

double *centre=new double[4*NXMAX*NYMAX];//pour mettre la position des centres translatés, on crée une variable 2D de la taille d'un plan apres tomo

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
double *fft_reel_shift_refocal=new double[4*NXMAX*NYMAX];
double *fft_imag_shift_refocal=new double[4*NXMAX*NYMAX];

///variable pour plan
double *plan_reel=new double[N];
double *plan_imag=new double[N];
double *plan_reel_shift=new double[N];
double *plan_imag_shift=new double[N];
double *plan_reel_norm=new double[dimFin.x*dimFin.y];
double *plan_imag_norm=new double[dimFin.x*dimFin.y];
double *plan_module_norm=new double[dimFin.x*dimFin.y];
double *plan_reel_norm_shift=new double[dimFin.x*dimFin.y];
double *plan_imag_norm_shift=new double[dimFin.x*dimFin.y];
//spectre repère centré
double *spectreC_reel=new double[dimFin.x*dimFin.y];
double *spectreC_imag=new double[dimFin.x*dimFin.y];
double *supRedon2D=new double[dimFin.x*dimFin.y];
//memset(supRedon2D, 0, dimFin.x*dimFin.y*8);

//spectre repère informatique
double  *spectreI_reel=new double[dimFin.x*dimFin.y];
double  *spectreI_imag=new double[dimFin.x*dimFin.y];
double  *spectreI_reel_fit=new double[dimFin.x*dimFin.y];
double  *spectreI_imag_fit=new double[dimFin.x*dimFin.y];

//spectre repère informatique total
double  *spectreI_reelTot=new double[dimFin.x*dimFin.y];
double  *spectreI_imagTot=new double[dimFin.x*dimFin.y];
//memset(spectreI_reelTot, 0, dimFin.x*dimFin.y*8);
//memset(spectreI_imagTot, 0, dimFin.x*dimFin.y*8);

double *distCal_fit= new double[dimFin.x*dimFin.y];
double *distCal_resharp= new double[dimFin.x*dimFin.y];
double *distRef=new double[dimFin.x*dimFin.y];
double *distCal_recal=new double[dimFin.x*dimFin.y];

double *plan_imag_normRef=new double[dimFin.x*dimFin.y];

double *phaseunwrap_recal=new double[dimFin.x*dimFin.y];
double *phaseMod1piRecal=new double[dimFin.x*dimFin.y];
double *phaseunwrapRecal=new double[dimFin.x*dimFin.y];

double *spectreDist=new double[dimFin.x*dimFin.y];
double *spectreDistTot=new double[dimFin.x*dimFin.y];
double *supRedon2D_dist=new double[dimFin.x*dimFin.y];

double *spectreI_reel_shift=new double[dimFin.x*dimFin.y];
double *spectreI_imag_shift=new double[dimFin.x*dimFin.y];

double *ampli_norm= new double[dimFin.x*dimFin.y];

///Chargement de la réference en module et calcul amplitude
//rempli_tableau(mod_ref,CheminModRef,coin,dimCCD);

float alpha=0.1;//coeff pour le masque de tuckey
double* masque=new double[DIMX_CCD2*DIMY_CCD2];
tukey2D(masque,DIMX_CCD2,DIMY_CCD2,alpha);

///######################variable pour algo de déroulement de phase########################################
int image_height=dimFin.x, image_width=dimFin.y;
int image_size = image_height * image_width;
int No_of_Edges = (image_width)*(image_height-1) + (image_width-1)*(image_height);
double *distCal = new double[image_size];
double *phaseMod1pi=new double[image_size];
double *phaseunwrap=new double[image_size];
PIXEL *Pixel = new PIXEL[image_size];
EDGE *Edge = new EDGE[No_of_Edges];
Mat masqueC=imread("vlsi.jpg",0);

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

double *phaseRef=new double[dimFin.x*dimFin.y];
double *ampliRef=new double[dimFin.x*dimFin.y];

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
    //cout<<TamponChemin<<endl;
    //cout<<test_existence<<endl;
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

		//calculPlan(holo1, holo2, holo3, holo4, plan_reel, plan_imag, DIMY_CCD2, coin, Fichier_holo, masque, mod_ref);
		///Calcul plan_reel et plan_imag avec critère local et global sur signal sur bruit
		calculPlanCritLG(holo1, holo2, holo3, holo4, plan_reel, plan_imag, DIMY_CCD2, coin, Fichier_holo, masque, mod_ref, cpt_angle);

///######################################################Circshift avant TF2D#######################################################
        Var2D dimCCD={DIMX_CCD2,DIMY_CCD2},decalCCD={DIMX_CCD2/2,DIMY_CCD2/2};
		circshift2(plan_reel,plan_reel_shift, dimCCD,decalCCD);
		circshift2(plan_imag,plan_imag_shift, dimCCD,decalCCD);

///########################################TF 2D du front d'onde: passage du plan image au plan réciproque##########################
		TF2D(plan_reel_shift,plan_imag_shift,fft_reel_tmp,fft_imag_tmp,DIMX_CCD2,DIMY_CCD2);

///######################################################Découpage dans fourier à NXMAX#############################################
        decalCoupe(fft_reel,fft_imag,fft_reel_tmp,fft_imag_tmp,dimDecoup,dimCCD);
        //SAV(fft_reel, 4*NXMAX*NYMAX, "/home/hui/maniptomo/IDP/TF2d/TF_reel_apres.bin", FLOAT,"a+b");

///######################################################circshift 2D apres TF######################################################
        Var2D dim2D={2*NXMAX,2*NYMAX},decal2D={NXMAX,NYMAX};
        circshift2(fft_reel,fft_reel_shift,dim2D,decal2D);
	    circshift2(fft_imag,fft_imag_shift,dim2D,decal2D);
        SAV(fft_reel_shift, 4*NXMAX*NYMAX, "/home/hui/maniptomo/IDP/champ_2d/TF_reel.bin", FLOAT,"a+b");
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

        //cout<<"Coordonnées speculaire"<<cpt_angle<<endl<<"xc="<<xc<<endl<<"yc="<<yc<<endl;
		float zc=round(sqrt(rayon*rayon-xc*xc-yc*yc));
		//cout<<"Coordonnées speculaire"<<"zc="<<zc<<endl;

		int xmij=2*NXMAX-xmi;//coordonnée objet jumeaux
		int ymij=2*NYMAX-ymi;
        //cout<<"Coordonnées jumeau angle "<<cpt_angle<<endl<<"xmij="<<xmij<<endl<<"ymij="<<ymij<<endl;
        //int cptJumo_max=2*NYMAX*ymij+xmij;

        centre[ymi*2*NXMAX+xmi]=cpt_angle;
        //cout<<"Coordonnées centre"<<cpt_angle<<endl<<"xmi="<<xmi<<endl<<"ymi="<<ymi<<endl;
///######################################################Virer l'objet jumeaux+ordre zero###########################################
		int t_mask_varZ=round(sqrt(xc*xc+yc*yc)/10)*2+1;
		double* cache_jumeauZ=new double[t_mask_varZ*t_mask_varZ];
		antigaussienne(cache_jumeauZ,t_mask_varZ,round(t_mask_varZ)*2,1,0);
        multiplier_masque2(fft_reel_shift, cache_jumeauZ, NXMAX, t_mask_varZ, NXMAX, NXMAX);//ordre zéro
        multiplier_masque2(fft_imag_shift, cache_jumeauZ, NXMAX, t_mask_varZ, NXMAX, NXMAX);

		if((xc*xc+yc*yc)>10)//35*35=1225;objet jumeau pas trop près de l'objet
		{ jumeau_elimine++;
            int t_mask_varJ=round(sqrt(xc*xc+yc*yc))*2+1;
           // cout<<"tmask="<<t_mask_var<<endl;
            double* cache_jumeauJ=new double[t_mask_varJ*t_mask_varJ];
            antigaussienne(cache_jumeauJ,t_mask_varJ,round(t_mask_varJ*2),1,0);//taille de mask
            multiplier_masque2(fft_imag_shift, cache_jumeauJ, NXMAX, t_mask_varJ, xmij, ymij);//jumeau
            multiplier_masque2(fft_reel_shift, cache_jumeauJ, NXMAX, t_mask_varJ, xmij, ymij);
            delete[] cache_jumeauJ;
		}//fin if objet jumeau pas trop pres
        delete[] cache_jumeauZ;

/*
///######################################################Virer l'objet jumeaux+ordre zero###########################################

        int t_mask_varZ=round(sqrt(xc*xc+yc*yc)/10)*2+1;
		double* cache_jumeauZ=new double[t_mask_varZ*t_mask_varZ];
		antigaussienne(cache_jumeauZ,t_mask_varZ,round(t_mask_varZ)*2,1,0);
        multiplier_masque2(fft_reel_shift, cache_jumeauZ, NXMAX, t_mask_varZ, NXMAX, NXMAX);//ordre zéro
        multiplier_masque2(fft_imag_shift, cache_jumeauZ, NXMAX, t_mask_varZ, NXMAX, NXMAX);

		int t_mask_var=2*NXMAX+1;
        double* cache_objet=new double[t_mask_var*t_mask_var];
        gaussienne(cache_objet,t_mask_var,round(t_mask_var*10),1,0);//taille de mask
        multiplier_masque2(fft_imag_shift, cache_objet, NXMAX, t_mask_var, xmi, ymi);//objet
        multiplier_masque2(fft_reel_shift, cache_objet, NXMAX, t_mask_var, xmi, ymi);
        delete[] cache_objet;
*/
///######################################################normalisation par le pic central###########################################
		for(int cpt=0;cpt<(4*NXMAX*NYMAX);cpt++)
			{
                fft_reel_shift_norm[cpt]=(fft_reel_shift[cpt]*max_part_reel+fft_imag_shift[cpt]*max_part_imag)/max_module;
                fft_imag_shift_norm[cpt]=(fft_imag_shift[cpt]*max_part_reel-fft_reel_shift[cpt]*max_part_imag)/max_module;
                //printf("fft_reel_shift_norm: %f, fft_imag_shift_norm: %f \n ",fft_reel_shift_norm[cpt], fft_imag_shift_norm[cpt]);
			}
        //SAV(fft_reel_shift_norm, 4*NXMAX*NYMAX, "/home/hui/maniptomo/IDP/champ_2d/fft_reel_shift_norm.bin", FLOAT,"a+b");

///########################################Utiliser spectre amplitude pour retrouver spectre refocalisée############################

        Var2D CoordSpecI={xmi,ymi};
        Autofocus(fft_reel_shift_norm, fft_imag_shift_norm, fft_reel_shift_refocal, fft_imag_shift_refocal, dimInit, rayon, CoordSpecI);
        changeDim2D(fft_reel_shift_refocal, spectreC_reel, dimInit, dimFin); ///utiliser les plans refocalisés
        changeDim2D(fft_imag_shift_refocal, spectreC_imag, dimInit, dimFin);

///######################################Utiliser spectre normal (sans refocaliser)#################################################

        //changeDim2D(fft_reel_shift_norm, spectreC_reel, dimInit, dimFin);
        //changeDim2D(fft_imag_shift_norm, spectreC_imag, dimInit, dimFin);
        //SAV(spectreC_reel, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/spectreC_reel.bin", FLOAT,"a+b");

///on redécale le spectre pour éliminer les franges d inclinaisons (-xc), on en profite pour revenir au zero informatique (+2*NXMAX)

        Var2D decal={-xc+dimFin.x/2,-yc+dimFin.x/2};
        circshift2(spectreC_reel,spectreI_reel,dimFin,decal);
        circshift2(spectreC_imag,spectreI_imag,dimFin,decal);

        Var2D fftDecal={dimFin.x/2,dimFin.x/2};
        circshift2(spectreI_reel,spectreI_reel_shift,dimFin,fftDecal);
        circshift2(spectreI_imag,spectreI_imag_shift,dimFin,fftDecal);
        SAV(spectreI_reel_shift, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/spectre2D.bin", FLOAT,"a+b");

        TF2D_INV(spectreI_reel,spectreI_imag,plan_reel_norm,plan_imag_norm,dimFin.x,dimFin.y);

        circshift2(plan_reel_norm,plan_reel_norm_shift,dimFin,fftDecal);
        circshift2(plan_imag_norm,plan_imag_norm_shift,dimFin,fftDecal);

///#####################################calcul module de champ2D####################################################################

        for(int cpt=0;cpt<(dimFin.x*dimFin.y);cpt++)
        {
            //plan_module_norm[cpt]=pow(plan_reel_norm_shift[cpt],2)+pow(plan_imag_norm_shift[cpt],2);  //Utiliser pour les plans non refocalisés
            plan_module_norm[cpt]=pow(plan_reel_norm[cpt],2)+pow(plan_imag_norm[cpt],2);  //Utiliser pour les plans refocalisés
            ampli_norm[cpt]=sqrt(plan_module_norm[cpt]);
        }
        SAV(ampli_norm, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/amplitude2D.bin", FLOAT,"a+b");

///################################################sans corrélation croisée#########################################################

        //phaseUnwrapping(plan_reel_norm_shift, plan_imag_norm_shift, phaseMod1pi, phaseunwrap, Pixel, Edge, dimFin.x); //Utiliser pour les plans non refocalisés
        phaseUnwrapping(plan_reel_norm, plan_imag_norm, phaseMod1pi, phaseunwrap, Pixel, Edge, dimFin.x); //Utiliser pour les plans refocalisés
        SAV(phaseMod1pi, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/phaseMod1pi.bin", FLOAT,"a+b");
        write_data("/home/hui/maniptomo/IDP/champ_2d/phaseUnwrapped.bin",phaseunwrap,image_size);

        for(int cpt=0;cpt<dimFin.x*dimFin.y;cpt++)
        {
            spectreI_reelTot[cpt]=spectreI_reelTot[cpt]+spectreI_reel[cpt];
            spectreI_imagTot[cpt]=spectreI_imagTot[cpt]+spectreI_imag[cpt];
            if(spectreI_reel[cpt]!=0 || spectreI_imag[cpt]!=0)
            supRedon2D[cpt]+=1;
        }

///#####################################corrélation croisée, retrouver plan reel et imaginaire######################################

//        double *plan_reel_norm_recal= new double[dimFin.x*dimFin.y];
//        double *plan_imag_norm_recal= new double[dimFin.x*dimFin.y];
//        recalAmpPha(ampli_norm, phaseMod1pi, ampliRef, phaseRef, dimFin, plan_reel_norm_recal, plan_imag_norm_recal, cpt_angle, posDecalAmpli, posDecalPhase, NbAngle, Gradient);
//        phaseUnwrapping(plan_reel_norm_recal, plan_imag_norm_recal, phaseMod1piRecal, phaseunwrapRecal, Pixel, Edge, dimFin.x);
//
//        double *spectreI_reel_recal= new double[dimFin.x*dimFin.y];
//        double *spectreI_imag_recal= new double[dimFin.x*dimFin.y];
//        TF2D(plan_reel_norm_recal,plan_imag_norm_recal,spectreI_reel_recal,spectreI_imag_recal,dimFin.x,dimFin.y);   ///utiliser plans recal avec correlation croisée amplitude et phase
//        delete[] plan_reel_norm_recal, plan_imag_norm_recal;
//
//        for(int cpt=0;cpt<dimFin.x*dimFin.y;cpt++)
//        {
//            spectreI_reelTot[cpt]=spectreI_reelTot[cpt]+spectreI_reel_recal[cpt];
//            spectreI_imagTot[cpt]=spectreI_imagTot[cpt]+spectreI_imag_recal[cpt];
//            if(spectreI_reel_recal[cpt]!=0 || spectreI_imag_recal[cpt]!=0)
//            supRedon2D[cpt]+=1;
//        }

 ///#####################################calcul la distance à partir de phase#######################################################

        //obj=sqrt(xc*xc+yc*yc);
        //coef=cos(asin(obj/NXMAX*NA/n1));
        coef=sqrt(rayon*rayon-xc*xc-yc*yc)/rayon;
        cout<<"coef="<<coef<<endl;

        for(int cpt=0; cpt<dimFin.x*dimFin.y; cpt++)
        {
           distCal[cpt] = phaseunwrap[cpt]*lambda*pow(10,9)/(4*PI*n1)/coef;
        }
        corrPhase(distCal,distCal_fit,dimFin.x, masqueC);   /// compensation aberration

        ///exportation distance
        SAV(distCal, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/distanceTAngles.bin", FLOAT,"a+b");
        SAV(distCal_fit, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/distanceTAngles_fit.bin", FLOAT,"a+b");

/*
///______________________________________________________________________________///

        for(int cpt=0;cpt<dimFin.x*dimFin.y;cpt++)
        {
            distCal_fit[cpt]=(double)sqrt(pow((distCal_fit[cpt]),2));
        }

        if(cpt_angle==1)
        {
            for(int cpt=0;cpt<dimFin.x*dimFin.y;cpt++)
            {
                distRef[cpt]=distCal_fit[cpt];
            }
            //SAV(distRef, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/distanceREF.bin", FLOAT,"a+b");
        }

        posDecal=corr_crois(distRef, distCal_fit, dimFin);
        //Tlérance sur le décalage. Si décalage supérieure à la tolérance, remttre à zero le décalage trouvé.
        decal2DGen(distCal_fit, distCal_recal,dimFin, posDecal);
        SAV(distCal_recal, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/distanceTAngles_recal.bin", FLOAT,"a+b");
///______________________________________________________________________________///
*/

SAV(supRedon2D, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/supRedon2D.bin", FLOAT,"a+b");
coefSomme=coef+coefSomme;
nomfinal++;}
else{
	printf("fichier %i inexistant\n",cpt_angle);
	//fclose(test_existence);
	}
}//fin de boucle for sur tous les angles  on peut désallouer les variables définit hors de la boucle

delete[] masque;
delete[] holo1, holo2, holo3, holo4;
delete[] fft_reel_shift_norm, fft_imag_shift_norm;
delete[] plan_reel, plan_imag, plan_reel_norm, plan_imag_norm, plan_module_norm, plan_reel_norm_shift, plan_imag_norm_shift;
delete[] spectreI_imag_shift, spectreI_reel_shift, ampli_norm;
delete[] distCal;//, distCal_fit, distCal_recal;
//delete[] fft_reel_shift, fft_imag_shift, fft_module_shift;
delete[] plan_reel_shift, plan_imag_shift;
delete[] spectreC_reel, spectreC_imag, spectreI_reel, spectreI_imag, spectreI_reel_fit, spectreI_imag_fit;

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

///###################################################################Normaliser par supredon###############################################################

///exportation du plan corrigé de supredon
for(int cpt=0;cpt<dimFin.x*dimFin.y;cpt++)
{   if(supRedon2D[cpt]!=0)
    {
        spectreI_reelTot[cpt]=spectreI_reelTot[cpt]/supRedon2D[cpt];
        spectreI_imagTot[cpt]=spectreI_imagTot[cpt]/supRedon2D[cpt];
    }
}

/////____________________________________________________________________________________
//int maxMod=0;
//double *spectreI_mod= new double[dimFin.x*dimFin.y];
//for(int cpt=0; cpt<dimFin.x*dimFin.y;cpt++)
//{
//    spectreI_mod[cpt]=pow(spectreI_reelTot[cpt],2)+pow(spectreI_imagTot[cpt],2);
//    if(spectreI_mod[cpt]>spectreI_mod[maxMod])
//    {
//        maxMod=cpt;
//    }
//}
//int xMod=maxMod%(dimFin.x/2);
//int yMod=maxMod/(dimFin.x/2);
//int xcMod=xMod-dimFin.x/4;
//int ycMod=yMod-dimFin.x/4;
//int t_mask_var=round(sqrt(xcMod*xcMod+ycMod*ycMod)/10)*2+1;
//double* cache_jumeau=new double[t_mask_var*t_mask_var];
//antigaussienne(cache_jumeau,t_mask_var,round(t_mask_var)*2,1,0);
//multiplier_masque2(spectreI_reelTot, cache_jumeau, dimFin.x/2, t_mask_var, xMod, yMod);//jumeau
//multiplier_masque2(spectreI_imagTot, cache_jumeau, dimFin.x/2, t_mask_var, xMod, yMod);//jumeau
/////____________________________________________________________________________________


///____________________________________________________________________________________
double  *spectreI_imagTot_shift=new double[dimFin.x*dimFin.y];
double  *spectreI_reelTot_shift=new double[dimFin.x*dimFin.y];
Var2D fftDecal={dimFin.x/2,dimFin.x/2};
circshift2(spectreI_reelTot,spectreI_reelTot_shift,dimFin,fftDecal);
circshift2(spectreI_imagTot,spectreI_imagTot_shift,dimFin,fftDecal);

double  *spectreI_moduleTot_shift=new double[dimFin.x*dimFin.y];
for(int cpt=0; cpt<dimFin.x*dimFin.y;cpt++)
{
    spectreI_moduleTot_shift[cpt]=pow(spectreI_imagTot_shift[cpt],2)+pow(spectreI_reelTot_shift[cpt],2);
}
SAV(spectreI_moduleTot_shift, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/spectresynth2D.bin", FLOAT,"a+b");

//float alpha1=0.9;//coeff pour le masque de tuckey
//double* masque1=new double[dimFin.x*dimFin.x];
//tukey2D(masque1,dimFin.x,dimFin.x,alpha1);
//for(int pixel=0;pixel<dimFin.x*dimFin.x;pixel++)
//{
//spectreI_reelTot_shift[pixel]=spectreI_reelTot_shift[pixel]*masque1[pixel];
//spectreI_imagTot_shift[pixel]=spectreI_imagTot_shift[pixel]*masque1[pixel];
//}

//SAV(spectreI_imagTot_shift, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/spectre_2d_apress.bin", FLOAT,"a+b");
//delete[] spectreI_imagTot_shift, spectreI_reelTot_shift;
///____________________________________________________________________________________


double  *plan_reel_normMoy=new double[dimFin.x*dimFin.y];
double  *plan_imag_normMoy=new double[dimFin.x*dimFin.y];

TF2D_INV(spectreI_reelTot,spectreI_imagTot,plan_reel_normMoy,plan_imag_normMoy,dimFin.x,dimFin.y);

double  *plan_reel_MoyFin=new double[dimFin.x*dimFin.y];
double  *plan_imag_MoyFin=new double[dimFin.x*dimFin.y];
circshift2(plan_imag_normMoy,plan_imag_MoyFin,dimFin,fftDecal);
circshift2(plan_reel_normMoy,plan_reel_MoyFin,dimFin,fftDecal);

///#####################################exportation Moy_champ_2D##################################################################
SAV(plan_reel_normMoy, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/Moy_champ2D.bin", FLOAT,"a+b");

///#####################################déroulement de phase après synthèse2D#######################################################
memset(phaseMod1pi, 0, dimFin.x*dimFin.y*8);
memset(phaseunwrap, 0, dimFin.x*dimFin.y*8);
phaseUnwrapping(plan_reel_normMoy, plan_imag_normMoy, phaseMod1pi, phaseunwrap, Pixel, Edge, dimFin.x);
//phaseUnwrapping(plan_reel_MoyFin, plan_imag_MoyFin, phaseMod1pi, phaseunwrap, Pixel, Edge, dimFin.x);
write_data("/home/hui/maniptomo/IDP/champ_2d/phaseSyn2D.bin",phaseMod1pi,image_size);
write_data("/home/hui/maniptomo/IDP/champ_2d/phaseUnwrappedSyn2D.bin",phaseunwrap,image_size);

//double  *fft_phase=new double[dimFin.x*dimFin.y];
//double  *fft_phase_shift=new double[dimFin.x*dimFin.y];
//TF2D_1P(phaseMod1pi,fft_phase,dimFin.x,dimFin.y);
//circshift2(fft_phase,fft_phase_shift,dimFin,fftDecal);
////write_data("/home/hui/maniptomo/IDP/champ_2d/fft_phase_shift_avant.bin",fft_phase_shift,image_size);
//
//int maskP=4*NXMAX;
//double* cache_maskP=new double[maskP*maskP];
//gaussienne(cache_maskP,maskP,round(maskP)*10,1,0);
//multiplier_masque2(fft_phase_shift, cache_maskP, 2*NXMAX, maskP, 2*NXMAX, 2*NXMAX);
//
//double  *phaseC=new double[dimFin.x*dimFin.y];
//double  *phaseC_shift=new double[dimFin.x*dimFin.y];
//TF2D_INV_1P(fft_phase,phaseC_shift,dimFin.x,dimFin.y);
//circshift2(phaseC_shift,phaseC,dimFin,fftDecal);
//
////write_data("/home/hui/maniptomo/IDP/champ_2d/fft_phase_shift_apres.bin",fft_phase_shift,image_size);
////write_data("/home/hui/maniptomo/IDP/champ_2d/phaseC.bin",phaseC_shift,image_size);
//delete[] fft_phase,fft_phase_shift;
//delete[] phaseC, phaseC_shift;

float coef_syn2D=coefSomme/nomfinal;
double *distCalSyn2D = new double[dimFin.x*dimFin.y];
for(int cpt=0; cpt<dimFin.x*dimFin.y; cpt++)
{
    distCalSyn2D[cpt] = phaseunwrap[cpt]*lambda*pow(10,9)/(4*PI*n1)/coef_syn2D;
}

double *distCalSyn2D_fit= new double[dimFin.x*dimFin.y];
corrPhase(distCalSyn2D,distCalSyn2D_fit,dimFin.x, masqueC);

///exportation distance
SAV(distCalSyn2D_fit, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/distanceSyn2D.bin", FLOAT,"a+b");

delete[] distCalSyn2D, distCalSyn2D_fit;
delete[] plan_reel_normMoy, plan_imag_normMoy;

///fin exporation plan corrigé supredon
delete[] spectreI_reelTot, spectreI_imagTot, supRedon2D;

///######################################################exportation de centre######################################################
SAV(centre, 4*NXMAX*NYMAX, "/home/hui/maniptomo/IDP/champ_2d/centre.bin", UINT,"wb");
printf("ecriture de centre en 32 bits\nDimension centre : %d x%d\n",2*NXMAX,2*NYMAX);

delete[] centre, ampli_ref, mod_ref;
delete[] fft_reel, fft_imag, fft_reel_tmp, fft_imag_tmp;
delete[] distCal, phaseMod1pi, phaseunwrap, Edge, Pixel;
//getchar();
//libération memoire allouée pour les threads
void fftw_cleanup_threads(void);
return 0;
}

