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
//valeur controlant l'exclusion de certains centres
int xm0_limite=200; //centre maximum
int ym0_limite=200;
int rayon_inf=xm0_limite*xm0_limite+ym0_limite*ym0_limite;
//valeur du coin pour la découpe
Var2D coin={430,200};
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
char Dossier_acquiz[]="/home/hui/maniptomo/reflexion_acquisition/2012120601/";
char NumSession[]="i";
char NomHoloPart1[strlen(Dossier_acquiz)+strlen("i")]; //au max 9999 hologramme
concatener(Dossier_acquiz,NumSession,NomHoloPart1);
char Fichier_holo[strlen(NomHoloPart1)+strlen("XXXX_p0.pgm")+1];//Nom final


//chemin pour référence d'amplitude
char Dossier_ref[]="/home/hui/maniptomo/reflexion_acquisition/2012120403/";
char CheminModRef[strlen(Dossier_ref)+strlen("mod_ref.bmp")+1];
concatener(Dossier_ref,"mod_ref.bmp",CheminModRef);

/*
///############################repertoire transmission###########################################
char Dossier_acquiz[]="/home/hui/maniptomo/image_acquisition/session12012501/";
char NumSession[]="session12012501-";
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
premier_plan=0,
Num_Angle_final=299,
// int sans_signal=0;
NbAngle=Num_Angle_final-premier_plan+1,
SautAngle=1,
N=DIMX_CCD2*DIMY_CCD2;//nombre de pixel (pour le tableau 2D)
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
centres_exclus=0,
nb_proj=0;
Var2D dimDecoup={NXMAX,NXMAX};

///reservation de 4 tableaux (image 2D)  pour le phase shifting
unsigned char* holo1=new unsigned char[N];
unsigned char* holo2=new unsigned char[N];
unsigned char* holo3=new unsigned char[N];
unsigned char* holo4=new unsigned char[N];

///variable pour phase
double *phaseMod1pi=new double[N];
double *centre=new double[4*NXMAX*NYMAX];//pour mettre la position des centres translatés, on crée une variable 2D de la taille d'un plan apres tomo
double *d1 = new double[4*NXMAX*4*NXMAX];

///variable pour image de référence
double* ampli_ref=new double[N]; //reservation reference et "noir" camera pour une "tentative" de correction de la ref
unsigned char* mod_ref=new unsigned char[N]; //reservation reference et "noir" camera pour une "tentative" de correction de la ref

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
double *plan_reel_norm_shift=new double[4*NXMAX*4*NYMAX];
double *plan_imag_norm_shift=new double[4*NXMAX*4*NYMAX];
double *plan_reel_norm=new double[4*NXMAX*4*NYMAX];
double *plan_imag_norm=new double[4*NXMAX*4*NYMAX];
double *plan_module_norm=new double[4*NXMAX*4*NYMAX];


///variable pour synthèse ouverture 2D
//on agrandit le support en doublant largeur et longueur (facteur 4)
Var2D dimInit={2*NXMAX,2*NYMAX},dimFin={4*NXMAX,4*NYMAX};
//spectre repère centré
double *spectreC_reel=new double[dimFin.x*dimFin.y];
double *spectreC_imag=new double[dimFin.x*dimFin.y];
double *supRedon2D=new double[dimFin.x*dimFin.y];
//memset(supRedon2D, 0, dimFin.x*dimFin.y*8);

//spectre repère informatique
double  *spectreI_reel=new double[dimFin.x*dimFin.y];
double  *spectreI_imag=new double[dimFin.x*dimFin.y];
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
rempli_tableau(mod_ref,CheminModRef,coin,dimCCD);

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

///variable pour algo de déroulement de phase
    double *WrappedImage, *UnwrappedImage;
    int image_height=4*NXMAX,
    image_width=4*NXMAX;
	int i, j;
	int image_size = image_height * image_width;
	int two_image_size = 2 * image_size;
	int No_of_Edges = (image_width)*(image_height-1) + (image_width-1)*(image_height);
	WrappedImage = (double *) calloc(image_size, sizeof(double));
	UnwrappedImage = (double *) calloc(image_size, sizeof(double));
	PIXEL *pixel = (PIXEL *) calloc(image_size, sizeof(PIXEL));
	EDGE *edge = (EDGE *) calloc(No_of_Edges, sizeof(EDGE));;




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
        memset(qualiteHolo, 0, N*8);
		for(int pixel=0;pixel<N;pixel++)
		{
            ///qualité des hologrammes
            qualiteHolo[pixel]=2*sqrt((pow(holo4[pixel]-holo2[pixel],2)+pow(holo1[pixel]-holo3[pixel],2)))/(holo1[pixel]+holo2[pixel]+holo3[pixel]+holo4[pixel]);
            //cout<<"contraste pour chaque pixel: "<<qualiteHolo[pixel]<<endl;
        }
        SAV(qualiteHolo, N, "/home/hui/maniptomo/IDP/champ_2d/qualiteHolo.bin", FLOAT,"a+b");
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
            sommeQholo=sommeQholo+qualiteHolo[pixel];

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

        memset(plan_reel, 0, N*8);
        memset(plan_imag, 0, N*8);

        //if(qualiteHolo[pixelQ_max]>0.8)
        if(moyQholo>0.4)
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
        memset(plan_reel_shift, 0, N*8);
        memset(plan_imag_shift, 0, N*8);
        Var2D dimCCD={DIMX_CCD2,DIMY_CCD2},decalCCD={DIMX_CCD2/2,DIMY_CCD2/2};
		circshift2(plan_reel,plan_reel_shift, dimCCD,decalCCD);
		circshift2(plan_imag,plan_imag_shift, dimCCD,decalCCD);

///########################################TF 2D du front d'onde: passage du plan image au plan réciproque##########################
        memset(fft_reel_tmp, 0, N*8);
        memset(fft_imag_tmp, 0, N*8);
		TF2D(plan_reel_shift,plan_imag_shift,fft_reel_tmp,fft_imag_tmp,DIMX_CCD2,DIMY_CCD2);


///######################################################Découpage dans fourier à NXMAX#############################################
        //SAV(fft_reel_tmp, N, "/home/hui/maniptomo/IDP/champ_2d/TF_reel.bin", FLOAT,"a+b");
        memset(fft_reel, 0, dimInit.x*dimInit.y*8);
        memset(fft_imag, 0, dimInit.x*dimInit.y*8);
        decalCoupe(fft_reel,fft_imag,fft_reel_tmp,fft_imag_tmp,dimDecoup,dimCCD);
        //SAV(fft_reel, 4*NXMAX*NYMAX, "/home/hui/maniptomo/IDP/TF2d/TF_reel_apres.bin", DOUBLE,"a+b");
		//printf("  NXMAX : %i, NYMAX : %i\n",NXMAX,NYMAX);

///######################################################circshift 2D apres TF######################################################
        memset(fft_reel_shift, 0, dimInit.x*dimInit.y*8);
        memset(fft_imag_shift, 0, dimInit.x*dimInit.y*8);
        Var2D dim2D={2*NXMAX,2*NYMAX},decal2D={NXMAX,NYMAX};
        circshift2(fft_reel,fft_reel_shift, dim2D,decal2D);
	    circshift2(fft_imag,fft_imag_shift, dim2D,decal2D);

///######################################################calcul maximum de spéculaire###############################################
        memset(fft_module_shift, 0, dimInit.x*dimInit.y*8);
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
        memset(fft_reel_shift_norm, 0, dimInit.x*dimInit.y*8);
        memset(fft_imag_shift_norm, 0, dimInit.x*dimInit.y*8);
		for(int cpt=0;cpt<(4*NXMAX*NYMAX);cpt++)
			{
                fft_reel_shift_norm[cpt]=(fft_reel_shift[cpt]*max_part_reel+fft_imag_shift[cpt]*max_part_imag)/max_module;
                fft_imag_shift_norm[cpt]=(fft_imag_shift[cpt]*max_part_reel-fft_reel_shift[cpt]*max_part_imag)/max_module;
                //printf("fft_reel_shift_norm: %f, fft_imag_shift_norm: %f \n ",fft_reel_shift_norm[cpt], fft_imag_shift_norm[cpt]);
			}

///########################Exportation du champ complexe pour vérification##########################################################
memset(spectreC_reel, 0, dimFin.x*dimFin.y*8);
memset(spectreC_imag, 0, dimFin.x*dimFin.y*8);
memset(spectreI_reel, 0, dimFin.x*dimFin.y*8);
memset(spectreI_imag, 0, dimFin.x*dimFin.y*8);
changeDim2D(fft_reel_shift_norm, spectreC_reel, dimInit, dimFin);
changeDim2D(fft_imag_shift_norm, spectreC_imag, dimInit, dimFin);

///on redécale le spectre pour éliminer les franges d inclinaisons (-xc), on en profite pour revenir au zero informatique (+2*NXMAX)
Var2D decal{-xc+2*NXMAX,-yc+2*NXMAX};
circshift2(spectreC_reel,spectreI_reel,dimFin,decal);
circshift2(spectreC_imag,spectreI_imag,dimFin,decal);

///#####################################exportation spectre_2D######################################################################
SAV(spectreC_reel, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/spectre2D.bin", FLOAT,"a+b");
memset(spectreI_reelTot, 0, dimFin.x*dimFin.y*8);
memset(spectreI_imagTot, 0, dimFin.x*dimFin.y*8);
memset(supRedon2D, 0, dimFin.x*dimFin.y*8);

for(int cpt=0;cpt<dimFin.x*dimFin.y;cpt++)
		{
		    spectreI_reelTot[cpt]=spectreI_reelTot[cpt]+spectreI_reel[cpt];
		    spectreI_imagTot[cpt]=spectreI_imagTot[cpt]+spectreI_imag[cpt];
		    if(spectreI_reel[cpt]!=0)
			supRedon2D[cpt]=supRedon2D[cpt]++;
		}
//TF2D
memset(plan_reel_norm_shift, 0, dimFin.x*dimFin.y*8);
memset(plan_imag_norm_shift, 0, dimFin.x*dimFin.y*8);
TF2D_INV(spectreI_reel,spectreI_imag,plan_reel_norm_shift,plan_imag_norm_shift,4*NXMAX,4*NYMAX);

memset(plan_reel_norm, 0, dimFin.x*dimFin.y*8);
memset(plan_imag_norm, 0, dimFin.x*dimFin.y*8);
Var2D fftDecal{2*NXMAX,2*NXMAX};
circshift2(plan_reel_norm_shift,plan_reel_norm,dimFin,fftDecal);
circshift2(plan_imag_norm_shift,plan_imag_norm,dimFin,fftDecal);

///#####################################calcul module de champ2D####################################################################
    for(int cpt=0;cpt<(4*NXMAX*4*NYMAX);cpt++)
        {
            plan_module_norm[cpt]=pow(plan_reel_norm[cpt],2)+pow(plan_imag_norm[cpt],2);
        }

///#####################################déroulement de phase sur tous les angles####################################################
    phaseUnwrapping(phaseMod1pi, plan_reel_norm, plan_imag_norm, 4*NXMAX, 4*NXMAX, WrappedImage,
                    UnwrappedImage, pixel, edge, image_width, image_height, image_size, No_of_Edges);
    SAV(phaseMod1pi, 4*NXMAX*4*NXMAX, "/home/hui/maniptomo/IDP/champ_2d/phaseMod1pi.bin", FLOAT,"a+b");
    write_data("/home/hui/maniptomo/IDP/champ_2d/phaseUnwrapped.bin",UnwrappedImage,image_size);

    for(int cpt=0; cpt<4*NXMAX*4*NXMAX; cpt++)
    {
        d1[cpt] = UnwrappedImage[cpt]*lambda*pow(10,9)/(4*PI*n1);
        //cout<<"d="<<d<<endl;
    }
    ///exportation distance
    SAV(d1, 4*NXMAX*4*NXMAX, "/home/hui/maniptomo/IDP/champ_2d/distanceTAngles.bin", FLOAT,"a+b");

	//	printf("max_part_reel_norm : %f, max_part_imag_norm: %f,  max_part_module_norm %f \n ",max_part_reel_norm, max_part_imag_norm, max_part_module_norm);
///exportation champ_2D
SAV(plan_module_norm, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/champ2D.bin.bin", FLOAT,"a+b");


///fin exportation du champ complexe#
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

delete[] supRedonQ;
delete[] sommeP_imag, sommeP_reel;
delete[] qualiteHolo;
delete[] plan_I13, plan_I42;
delete[] masque;
delete[] holo1, holo2, holo3, holo4;
delete[] fft_reel_shift_norm;
delete[] fft_imag_shift_norm;
delete[] plan_reel;
delete[] plan_imag;
delete[] plan_reel_norm;
delete[] plan_imag_norm;
delete[] plan_module_norm;
delete[] plan_reel_norm_shift;
delete[] plan_imag_norm_shift;
delete[] d1;
delete[] fft_reel_shift, fft_imag_shift, fft_module_shift;
delete[] plan_reel_shift, plan_imag_shift;
delete[] spectreC_reel;
delete[] spectreC_imag;
delete[] spectreI_reel;
delete[] spectreI_imag;
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
        double  *plan_reel_normMoy=new double[4*NXMAX*4*NYMAX];
        double  *plan_imag_normMoy=new double[4*NXMAX*4*NYMAX];

    //SAV(spectreI_imagTot, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/spectre_2d_apress.bin", FLOAT,"a+b");

        TF2D_INV(spectreI_reelTot,spectreI_imagTot,plan_reel_normMoy,plan_imag_normMoy,4*NXMAX,4*NYMAX);

        double  *plan_reel_MoyFin=new double[4*NXMAX*4*NYMAX];
        double  *plan_imag_MoyFin=new double[4*NXMAX*4*NYMAX];

        Var2D fftDecal{2*NXMAX,2*NXMAX};

        circshift2(plan_imag_normMoy,plan_imag_MoyFin,dimFin,fftDecal);
        circshift2(plan_reel_normMoy,plan_reel_MoyFin,dimFin,fftDecal);

///#####################################exportation Moy_champ_2D##################################################################
    SAV(plan_imag_MoyFin, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/Moy_champ2D.bin", FLOAT,"a+b");

///#####################################déroulement de phase après synthèse2D#######################################################
    phaseUnwrapping(phaseMod1pi, plan_reel_MoyFin, plan_imag_MoyFin, 4*NXMAX, 4*NXMAX, WrappedImage,
                    UnwrappedImage, pixel, edge, image_width, image_height, image_size, No_of_Edges);
    write_data("/home/hui/maniptomo/IDP/champ_2d/phaseUnwrappedSyn2D.bin",UnwrappedImage,image_size);
    double *d2 = new double[4*NXMAX*4*NXMAX];
    for(int cpt=0; cpt<4*NXMAX*4*NXMAX; cpt++)
    {
        d2[cpt] = UnwrappedImage[cpt]*lambda*pow(10,9)/(4*PI*n1);
        //cout<<"d2="<<d<<endl;
    }
    ///exportation distance
    SAV(d2, 4*NXMAX*4*NXMAX, "/home/hui/maniptomo/IDP/champ_2d/distanceSyn2D.bin", FLOAT,"a+b");
    delete[] d2;

        delete[] spectreI_reelTot;
        delete[] spectreI_imagTot;
        delete[] supRedon2D;
        delete[] plan_reel_MoyFin;
        delete[] plan_imag_MoyFin;
        delete[] plan_reel_normMoy;
        delete[] plan_imag_normMoy;

///fin exporation plan corrigé supredon

///######################################################exportation de centre######################################################
    SAV(centre, 4*NXMAX*NYMAX, "/home/hui/maniptomo/centre.bin", UINT,"wb");
    printf("ecriture de centre en 32 bits\nDimension centre : %d x%d\n",2*NXMAX,2*NYMAX);
    delete[] centre;
	delete[] ampli_ref, mod_ref;
	delete[] fft_reel;
	delete[] fft_imag;
	delete[] fft_reel_tmp;
	delete[] fft_imag_tmp;
    delete[] phaseMod1pi;
    free(edge);
	free(pixel);
	free(UnwrappedImage);
	free(WrappedImage);
	//getchar();

//libération memoire allouée pour les threads
void fftw_cleanup_threads(void);
return 0;
}
