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
int DIMX_CCD2=508,  DIMY_CCD2=508;
Var2D dimCCD={DIMX_CCD2,DIMY_CCD2};
const int N=DIMX_CCD2*DIMY_CCD2;//nombre de pixel (pour le tableau 2D)
//valeur controlant l'exclusion de certains centres
int xm0_limite=200; //centre maximum
int ym0_limite=200;
int rayon_inf=xm0_limite*xm0_limite+ym0_limite*ym0_limite;
//valeur du coin pour la découpe
Var2D coin={440,220};
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

/*
///###################################repertoire reflexion#########################################
//char Dossier_acquiz[]="/home/hui/maniptomo/reflexion_acquisition/2012031301/";
//char Dossier_acquiz[]="/home/hui/maniptomo/reflexion_acquisition/2012070502/IMAGE/";
char Dossier_acquiz[]="/home/hui/maniptomo/reflexion_acquisition/2012120602/";
char NumSession[]="i";
char NomHoloPart1[strlen(Dossier_acquiz)+strlen("i")]; //au max 9999 hologramme
concatener(Dossier_acquiz,NumSession,NomHoloPart1);
char Fichier_holo[strlen(NomHoloPart1)+strlen("XXXX_p0.pgm")+1];//Nom final
*/
/*
//chemin pour référence d'amplitude
char Dossier_ref[]="/home/hui/maniptomo/reflexion_acquisition/2012120403/";
char CheminModRef[strlen(Dossier_ref)+strlen("mod_ref.bmp")+1];
concatener(Dossier_ref,"mod_ref.bmp",CheminModRef);
*/

///############################repertoire transmission###########################################
char Dossier_acquiz[]="/home/hui/maniptomo/image_acquisition/session13020602/";
char NumSession[]="session13020602-";
char NomHoloPart1[strlen(Dossier_acquiz)+strlen("sessionXXXXXXXX")]; //au max 9999 hologramme
concatener(Dossier_acquiz,NumSession,NomHoloPart1);
char Fichier_holo[strlen(NomHoloPart1)+strlen("recordXXXX-004.bmp")+1];//Nom final

//chemin pour référence d'amplitude
char Dossier_ref[]="/home/hui/maniptomo/image_acquisition/session13020601/";
char CheminModRef[strlen(Dossier_ref)+strlen("mod_ref.bmp")+1];
concatener(Dossier_ref,"mod_ref.bmp",CheminModRef);
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
premier_plan=8,
Num_Angle_final=57,
// int sans_signal=0;
NbAngle=Num_Angle_final-premier_plan+1,
SautAngle=1;

const float
n1=1.515,	//indice de l'huile
NA=1.40,	//ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon)
theta=asin(NA/n1),//theta_max
lambda=475*pow(10,-9), //longueur d'onde
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
double *plan_reel_norm_fit=new double[dimFin.x*dimFin.y];
double *plan_imag_norm_fit=new double[dimFin.x*dimFin.y];

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

///variables qualité holo
double *qualiteHolo=new double[N];
double *sommeP_reel=new double[N];
double *sommeP_imag=new double[N];
double *plan_I42=new double[N];
double *plan_I13=new double[N];
double *supRedonQ=new double[N];

double *reel_arc=new double[N_tab];
double *imag_arc=new double[N_tab];
double *sup_redon=new double[N_tab];

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
Mat masqueC=imread("Contours3.jpg",0);

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
/*
///#################### reflexion############################
    char TamponChemin[sizeof(Fichier_holo)];
	char Num_angle[4];
	sprintf(Num_angle,"%04i",cpt_angle);
	//créer le chemin final
	concatener(NomHoloPart1,Num_angle,Fichier_holo);
    strcpy(TamponChemin,Fichier_holo);
	concatener(TamponChemin,"_p0.pgm",TamponChemin);
	test_existence = fopen(TamponChemin, "rb");
*/

    ///###################transmission#######################
    char TamponChemin[sizeof(Fichier_holo)];
	char Num_angle[4];
	sprintf(Num_angle,"record%i",cpt_angle);
	//créer le chemin final
	concatener(NomHoloPart1,Num_angle,Fichier_holo);
    strcpy(TamponChemin,Fichier_holo);
	concatener(TamponChemin,"-001.bmp",TamponChemin);
	test_existence = fopen(TamponChemin, "rb");

    if(test_existence!=NULL)
    {
		fclose(test_existence);
		nb_proj++;

		holo1=charger_image2D(holo1,1,Fichier_holo, coin.x, coin.y,DIMX_CCD2,DIMY_CCD2);
		holo2=charger_image2D(holo2,2,Fichier_holo, coin.x, coin.y,DIMX_CCD2,DIMY_CCD2);
		holo3=charger_image2D(holo3,3,Fichier_holo, coin.x, coin.y,DIMX_CCD2,DIMY_CCD2);
		holo4=charger_image2D(holo4,4,Fichier_holo, coin.x, coin.y,DIMX_CCD2,DIMY_CCD2);
/*
///######################################################décalage de phase##########################################################
        memset(plan_reel, 0, N*8);
        memset(plan_imag, 0, N*8);
		for(int pixel=0;pixel<N;pixel++)
		{
            plan_reel[pixel]=((double)holo1[pixel]-(double)holo3[pixel])*(masque[pixel]);
            plan_imag[pixel]=((double)holo4[pixel]-(double)holo2[pixel])*(masque[pixel]);
		}
        SAV(plan_reel, N, "/home/hui/maniptomo/IDP/champ_2d/plan_reel.bin", FLOAT,"a+b");
*/

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
        //int pixelQ_max=0;
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

            if(qualiteHolo[pixel]>0)
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

        memset(plan_reel, 0, N*8);
        memset(plan_imag, 0, N*8);

        //if(qualiteHolo[pixelQ_max]>0.8)
        if(moyQholo>0)
        {

        for(int pixel=0;pixel<N;pixel++)
        {
            if(supRedonQ[pixel]==0)
            supRedonQ[pixel]=1;
            //plan_reel[pixel]=sommeP_reel[pixel]*(masque[pixel])/supRedonQ[pixel];
            //plan_imag[pixel]=sommeP_imag[pixel]*(masque[pixel])/supRedonQ[pixel];

            ampli_ref[pixel]= sqrt((double)mod_ref[pixel]);
            plan_reel[pixel]=sommeP_reel[pixel]*(masque[pixel])/supRedonQ[pixel]/ampli_ref[pixel];
            plan_imag[pixel]=sommeP_imag[pixel]*(masque[pixel])/supRedonQ[pixel]/ampli_ref[pixel];

        }
        SAV(plan_reel, N, "/home/hui/maniptomo/IDP/champ_2d/plan_reel.bin", FLOAT,"a+b");

///######################################################Circshift avant TF2D#######################################################
        Var2D dimCCD={DIMX_CCD2,DIMY_CCD2},decalCCD={DIMX_CCD2/2,DIMY_CCD2/2};
		circshift2(plan_reel,plan_reel_shift, dimCCD,decalCCD);
		circshift2(plan_imag,plan_imag_shift, dimCCD,decalCCD);

///########################################TF 2D du front d'onde: passage du plan image au plan réciproque##########################
		TF2D(plan_reel_shift,plan_imag_shift,fft_reel_tmp,fft_imag_tmp,DIMX_CCD2,DIMY_CCD2);

///######################################################Découpage dans fourier à NXMAX#############################################
        SAV(fft_reel_tmp, N, "/home/hui/maniptomo/IDP/champ_2d/TF_reel.bin", FLOAT,"a+b");
        decalCoupe(fft_reel,fft_imag,fft_reel_tmp,fft_imag_tmp,dimDecoup,dimCCD);
        //SAV(fft_reel, 4*NXMAX*NYMAX, "/home/hui/maniptomo/IDP/TF2d/TF_reel_apres.bin", FLOAT,"a+b");

///######################################################circshift 2D apres TF######################################################
        Var2D dim2D={2*NXMAX,2*NYMAX},decal2D={NXMAX,NYMAX};
        circshift2(fft_reel,fft_reel_shift,dim2D,decal2D);
	    circshift2(fft_imag,fft_imag_shift,dim2D,decal2D);

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
        //int cptJumo_max=2*NYMAX*ymij+xmij;
        int Obj= xc*xc+yc*yc;
        int Jum= xmij*xmij+ymij*ymij;
        //cout<<"objet= \n"<<Obj<<"Jumeau="<<Jum<<endl;

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
            antigaussienne(cache_jumeauJ,t_mask_varJ,round(t_mask_varJ*10),1,0);//taille de mask
            multiplier_masque2(fft_imag_shift, cache_jumeauJ, NXMAX, t_mask_varJ, xmij, ymij);//jumeau
            multiplier_masque2(fft_reel_shift, cache_jumeauJ, NXMAX, t_mask_varJ, xmij, ymij);
            delete[] cache_jumeauJ;
		}//fin if objet jumeau pas trop pres
        delete[] cache_jumeauZ;

///######################################################normalisation par le pic central###########################################
		for(int cpt=0;cpt<(4*NXMAX*NYMAX);cpt++)
			{
                fft_reel_shift_norm[cpt]=(fft_reel_shift[cpt]*max_part_reel+fft_imag_shift[cpt]*max_part_imag)/max_module;
                fft_imag_shift_norm[cpt]=(fft_imag_shift[cpt]*max_part_reel-fft_reel_shift[cpt]*max_part_imag)/max_module;
                //printf("fft_reel_shift_norm: %f, fft_imag_shift_norm: %f \n ",fft_reel_shift_norm[cpt], fft_imag_shift_norm[cpt]);
			}
        SAV(fft_reel_shift_norm, 4*NXMAX*NYMAX, "/home/hui/maniptomo/IDP/champ_2d/fft_reel_shift_norm.bin", FLOAT,"a+b");
///########################Exportation du champ complexe pour vérification##########################################################
        changeDim2D(fft_reel_shift_norm, spectreC_reel, dimInit, dimFin);
        changeDim2D(fft_imag_shift_norm, spectreC_imag, dimInit, dimFin);

///on redécale le spectre pour éliminer les franges d inclinaisons (-xc), on en profite pour revenir au zero informatique (+2*NXMAX)
        Var2D decal={-xc+2*NXMAX,-yc+2*NXMAX};
        circshift2(spectreC_reel,spectreI_reel,dimFin,decal);
        circshift2(spectreC_imag,spectreI_imag,dimFin,decal);

        TF2D_INV(spectreI_reel,spectreI_imag,plan_reel_norm_shift,plan_imag_norm_shift,dimFin.x,dimFin.y);

        Var2D fftDecal={2*NXMAX,2*NXMAX};
        circshift2(plan_reel_norm_shift,plan_reel_norm,dimFin,fftDecal);
        circshift2(plan_imag_norm_shift,plan_imag_norm,dimFin,fftDecal);

///####################################### Correction de fond avec ajustement polynomiale ##########################################
        corrplan_I(plan_reel_norm, plan_imag_norm, plan_reel_norm_fit, plan_imag_norm_fit, dimFin.x, masqueC);

///#####################################calcul module de champ2D####################################################################
        for(int cpt=0;cpt<(dimFin.x*dimFin.y);cpt++)
        {
            plan_module_norm[cpt]=pow(plan_reel_norm[cpt],2)+pow(plan_imag_norm[cpt],2);
        }

///#####################################déroulement de phase sur tous les angles####################################################
        phaseUnwrapping(plan_reel_norm_fit, plan_imag_norm_fit, phaseMod1pi, phaseunwrap, Pixel, Edge, dimFin.x);
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
        SAV(plan_reel_norm_fit, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/champ2D.bin", FLOAT,"a+b");

///#####################################exportation spectre_2D######################################################################
        TF2D(plan_reel_norm_fit,plan_imag_norm_fit,spectreI_reel_fit,spectreI_imag_fit,dimFin.x,dimFin.y);
        //SAV(spectreC_reel, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/spectre2D.bin", FLOAT,"a+b");
        for(int cpt=0;cpt<dimFin.x*dimFin.y;cpt++)
        {
            spectreI_reelTot[cpt]=spectreI_reelTot[cpt]+spectreI_reel_fit[cpt];
            spectreI_imagTot[cpt]=spectreI_imagTot[cpt]+spectreI_imag_fit[cpt];
            if(spectreI_reel_fit[cpt]!=0)
            supRedon2D[cpt]=supRedon2D[cpt]++;
        }

    } //fin de test de validite du nom de fichier
    else{
        printf("qualité d'hologramme %i insuffisante\n",cpt_angle);
        }
}
else{
	printf("fichier %i inexistant\n",cpt_angle);
	//fclose(test_existence);
	}
}//fin de boucle for sur tous les angles  on peut désallouer les variables définit hors de la boucle
delete[] masque;
delete[] supRedonQ, sommeP_imag, sommeP_reel, qualiteHolo, plan_I13, plan_I42;
delete[] holo1, holo2, holo3, holo4;
delete[] fft_reel_shift_norm, fft_imag_shift_norm, plan_reel_norm_fit, plan_imag_norm_fit;
delete[] plan_reel, plan_imag, plan_reel_norm, plan_imag_norm, plan_module_norm, plan_reel_norm_shift, plan_imag_norm_shift;
delete[] distCal;
delete[] fft_reel_shift, fft_imag_shift, fft_module_shift;
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

///exportation du plan corrigé de supredon
for(int cpt=0;cpt<dimFin.x*dimFin.y;cpt++)
{   if(supRedon2D[cpt]!=0)
    {
        spectreI_reelTot[cpt]=spectreI_reelTot[cpt]/supRedon2D[cpt];
        spectreI_imagTot[cpt]=spectreI_imagTot[cpt]/supRedon2D[cpt];
    }
}
double  *plan_reel_normMoy=new double[dimFin.x*dimFin.y];
double  *plan_imag_normMoy=new double[dimFin.x*dimFin.y];

//SAV(spectreI_imagTot, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/spectre_2d_apress.bin", FLOAT,"a+b");

TF2D_INV(spectreI_reelTot,spectreI_imagTot,plan_reel_normMoy,plan_imag_normMoy,dimFin.x,dimFin.y);

double  *plan_reel_MoyFin=new double[dimFin.x*dimFin.y];
double  *plan_imag_MoyFin=new double[dimFin.x*dimFin.y];

Var2D fftDecal={2*NXMAX,2*NXMAX};

circshift2(plan_imag_normMoy,plan_imag_MoyFin,dimFin,fftDecal);
circshift2(plan_reel_normMoy,plan_reel_MoyFin,dimFin,fftDecal);

///#####################################exportation Moy_champ_2D##################################################################
SAV(plan_imag_MoyFin, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/Moy_champ2D.bin", FLOAT,"a+b");

///#####################################déroulement de phase après synthèse2D#######################################################
memset(phaseMod1pi, 0, dimFin.x*dimFin.y*8);
memset(phaseunwrap, 0, dimFin.x*dimFin.y*8);
phaseUnwrapping(plan_reel_MoyFin, plan_imag_MoyFin, phaseMod1pi, phaseunwrap, Pixel, Edge, dimFin.x);
write_data("/home/hui/maniptomo/IDP/champ_2d/phaseUnwrappedSyn2D.bin",phaseunwrap,image_size);
double *distCalSyn2D = new double[dimFin.x*dimFin.y];
for(int cpt=0; cpt<dimFin.x*dimFin.y; cpt++)
{
    distCalSyn2D[cpt] = phaseunwrap[cpt]*lambda*pow(10,9)/(4*PI*n1);
    //cout<<"d2="<<d<<endl;
}
///exportation distance
SAV(distCalSyn2D, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/distanceSyn2D.bin", FLOAT,"a+b");

delete[] distCalSyn2D, spectreI_reelTot, spectreI_imagTot, supRedon2D;
delete[] plan_reel_MoyFin, plan_imag_MoyFin, plan_reel_normMoy, plan_imag_normMoy;

///fin exporation plan corrigé supredon

///######################################################exportation de centre######################################################
SAV(centre, 4*NXMAX*NYMAX, "/home/hui/maniptomo/centre.bin", UINT,"wb");
printf("ecriture de centre en 32 bits\nDimension centre : %d x%d\n",2*NXMAX,2*NYMAX);

delete[] centre, ampli_ref, mod_ref;
delete[] fft_reel, fft_imag, fft_reel_tmp, fft_imag_tmp;
delete[] distCal, phaseMod1pi, phaseunwrap, Edge, Pixel;
//getchar();
//libération memoire allouée pour les threads
void fftw_cleanup_threads(void);
return 0;
}

