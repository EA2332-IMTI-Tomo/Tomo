#include "projet.h"

int main(int argc, char *argv[])
{

  Var2D dimCCD= {WINDOW_X, WINDOW_Y};
  Var2D decalCCD= {dimCCD.x/2,dimCCD.y/2};

  int DIMX_CCD2=dimCCD.x,  DIMY_CCD2=dimCCD.y;
  int NbPixROI2d=dimCCD.x*dimCCD.y;//nombre de pixel (pour le tableau 2D)

  //Forcer la dilmension de l'espace à une certaine taille, attention, ceci changera lataille des pixels tomo
  const int   dim_final=512,//peut être différent de 4*NXMAX, mais l'image final sera zoomée;
    N_tab=dim_final*dim_final*dim_final;//Nombre de pixel dans l'espace final (indice max du tableau 3D+1)

  
  // 
  IplImage *tampon_image_cv = cvCreateImage( cvSize(dimCCD.x, dimCCD.y), IPL_DEPTH_8U, 1); 


  ///--------------Chemin fichiers-------------------------------------
  char* input_dir=INPUT_DIR;
  string Dossier_acquiz=input_dir , CheminModRef=Dossier_acquiz+"mod_ref.bmp";
  string NumSession=RADIX_SESSION;
  string NomHoloPart1=Dossier_acquiz+NumSession;


  //valeur controlant l'exclusion de certains centres
  int xm0_limite=200,  ym0_limite=200;//centre maximum
  int rayon_inf=xm0_limite*xm0_limite+ym0_limite*ym0_limite;

  Var2D coin {CUT_ORIGIN_X,CUT_ORIGIN_Y}; //valeur du coin pour la découpe
  if(coin.x+dimCCD.x>IMG_DIMX || coin.y+dimCCD.y>IMG_DIMY)//largeur : 740, hauteur : 574.
    printf("!!!!!!!!!!!!Attention, la zone decoupee sort de l'image\n**\n****\n******\n");

  int jumeau_elimine=0, points_faux=0;
  int fftwThreadInit;
  fftwThreadInit=fftw_init_threads();
  int nthreads=4;
  printf("fftwThreadInit: %i\n",fftwThreadInit);

  ///--selections des hologrammes par numéro
  const int   premier_plan=1, Num_Angle_final= LAST_ANGLE_NUM,//
    NbAngle=Num_Angle_final-premier_plan+1, SautAngle=1;


  ///----------------------------CONSTANTES EXPERIMENTALES->PREPROC?--------------------------------------------------
  const float //indice huile,ouverture numerique de l'objectif? (celle du condenseur intervient sur la forme, la taille, du papillon),//theta_max
    n1=1.515, NA=1.40, theta=asin(NA/n1), lambda=632*pow(10,-9), G=100,	TpCam=8*pow(10,-6), //longueur d'onde, //grossissement telan+objectif,Taille pixel Caméra
    Rf=1,//1.2;=1/facteur grossissement
    tailleTheoPixelHolo=TpCam/(G/Rf)*pow(10,9);//Pour info
  double
    rayon=n1*Rf*TpCam*DIMX_CCD2/(G*lambda); //Rayon
  printf("theta : %f\n",theta);
  printf("Rayon %f \n",rayon);
  cout<<"############__--- Grandissement optique ---__###############################"<<endl;
  cout<<"Grandissement optique = "<<G/Rf<<" ->taille des pixels espace holo : "<<tailleTheoPixelHolo<<" nm"<<endl;
  const int
    NXMAX=round(n1*Rf*TpCam*DIMX_CCD2/(G*lambda)*NA/n1), NYMAX=NXMAX;//NXMAX=round(rayon*NA/n1)=round(rayon*sin theta_max);
  cout<<"NA/n1 : "<<NA/n1<<"NXMAX : "<<NA/n1*rayon<<endl;
  float
    tailleTheoPixelTomo=tailleTheoPixelHolo*DIMX_CCD2/(4*NXMAX);
  cout<<"Taille théorique pixel tomo sans zoom:"<<tailleTheoPixelTomo<<endl;
  cout<<endl<<"#############__--- Zoom numérique ---__###############################"<<endl<<endl;
  printf("NXMAX=%i,soit extension dans Fourier=4*NXMAX= %i \n",NXMAX,4*NXMAX);
  double zoom=double(dim_final)/double(4*n1*Rf*TpCam*DIMX_CCD2/(G*lambda)*NA/n1);
  if(zoom!=1) {
    cout<<"dim_final forcée à "<<dim_final<<" au lieu de "<<4*NXMAX<<"-->facteur de zoom="<<zoom<<endl;
    cout<<"nouvelle taille pixel tomo imposée par le zoom numérique= "<<tailleTheoPixelTomo/zoom<<" nm"<<endl;
    cout<<endl<<"###########################################"<<endl;
  }

  cout <<"Dim finale : " << round(4*zoom*n1*TpCam*DIMX_CCD2/(G*lambda)*NA/n1) << " pixels cube\n";
  int   centres_exclus=0,   nb_proj=0;
  //printf("Précision d'exportation: %i bits\n chemin :%s \n",8*sizeof(precision),Chemin);


  ///-----------------------------CALCULS CONSTANTES : éviter de les calculer dans la boucle--------------
  const int   dimVolX=round(dim_final),
    decalX=round(dim_final/2),decalY=round(dim_final/2),///demi longueur pour recalculer les indices
    dimPlanFinal=round(dim_final*dim_final), decalZ=round(dim_final/2);
  printf("dimx espace 3D, dimVolX: %i\n",dimVolX);
  fflush(stdout);
  cout.flush();

  //reservation de 4 tableaux (image 2D)  pour le phase shifting
  unsigned char* holo1=new unsigned char[NbPixROI2d];
  unsigned char* holo2=new unsigned char[NbPixROI2d];
  unsigned char* holo3=new unsigned char[NbPixROI2d];
  unsigned char* holo4=new unsigned char[NbPixROI2d];

  double* ampli_ref=new double[NbPixROI2d]; //reservation reference et "noir" camera pour une "tentative" de correction de la ref
  // 	unsigned char* noir_camera=new unsigned char[NbPixROI2d];
  double* masque=new double[NbPixROI2d]; //masque hamming

  nbCplx *plan=new nbCplx[NbPixROI2d], *plan_shift=new nbCplx[NbPixROI2d];
  ///variables utilisées pour la TF2D
  nbCplx *fft_tmp=new nbCplx[NbPixROI2d]; //Avant Crop à NXMAX;
  nbCplx *fft=new nbCplx[NbPixROI2d], *fft_shift=new nbCplx[4*NXMAX*NYMAX];//Après Crop à NXMAX;
  nbCplx *fft_shift_norm=new nbCplx[NbPixROI2d];
  double *fft_module_shift=new double[4*NXMAX*NYMAX];

  //double *qualite_phase_shift=new double[NbPixROI2d];
  //espace3D reel, imaginaire et support de redondance

  double *spectre3D_Re=new double[N_tab];//partie reelle du papillon dans l'espace reciproque
  double *spectre3D_Im=new double[N_tab];//partie imaginaire du papillon dans l'espace reciproque
  double *sup_redon=new double[N_tab];//pour normaliser par rapport au nombre de fois où l'on remplit la frequence
  double *centre=new double[4*NXMAX*NYMAX];//pour mettre la position des centres translatés, on crée une variable 2D de la taille d'un plan apres tomo


  ///-------------Réservation chrono----------------
  clock_t //temps en microseconde
    temps_depart= clock(), temps_initial=clock (), temps_final, temps_arrivee;
  float
    temps_cpu=0, temps_total=0;    /* temps total en secondes */


  ///-------------Chargement de la réference en module et calcul amplitude
  unsigned char* mod_ref=new unsigned char[NbPixROI2d]; //reservation reference et "noir" camera pour une "tentative" de correction de la ref
  // charger_image2D(mod_ref,CheminModRef,coin,dimCCD);

  //charge_decoupe_image_cv(mod_ref, CheminModRef, tampon_image_cv, coin.x, coin.y, dimCCD.x, dimCCD.y);

  for(int pixel=0; pixel<NbPixROI2d; pixel++) {
    ampli_ref[pixel]=1; //sqrt((double)mod_ref[pixel]);
  }


  ///---variable pour le masquage (fenêtre de tukey, "anti" gaussienne).
  float alpha=0.1;///coeff pour le masque de tuckey
  masque=tukey2D(DIMX_CCD2,DIMY_CCD2,alpha);
  const int t_mask=31; //taille totale du masque en pixel pour elimination objet jumeau;
  double* cache_jumeau2=new double[t_mask*t_mask];
  antigaussienne(cache_jumeau2,t_mask,round(t_mask*2),1,0);

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

  ///--------------------------------------BALAYAGE ANGULAIRE/////////////////////////
  char charAngle[4+1];
  double coef_IDP=0;

  for(int cpt_angle=premier_plan; cpt_angle<premier_plan+NbAngle; cpt_angle=cpt_angle+SautAngle) { //boucle sur tous les angles
    //cout<<"cpt_angle : "<<cpt_angle<<endl;
    if((cpt_angle-100*(cpt_angle/100))==0)
      printf("cpt_angle=%i\n",cpt_angle);
    ///Concaténer le numéro d'angle dans le nom
    sprintf(charAngle,"%i",cpt_angle);
    string imgExiste=Dossier_acquiz+NumSession+charAngle+"-001.bmp";
    string cheminHolo=Dossier_acquiz+NumSession+charAngle;

    ///-----------------------___ChARGE HOLO___---------------------------
    char Num_angle[4];
    sprintf(Num_angle,"record%i",cpt_angle);

    test_existence = fopen(imgExiste.c_str(), "rb");
    if(test_existence!=NULL) {
      fclose(test_existence);

      charge_decoupe_image_cv(holo1, cheminHolo + "-001" + ".bmp", tampon_image_cv, coin.x, coin.y, dimCCD.x, dimCCD.y);
      charge_decoupe_image_cv(holo2, cheminHolo + "-002" + ".bmp", tampon_image_cv, coin.x, coin.y, dimCCD.x, dimCCD.y);
      charge_decoupe_image_cv(holo3, cheminHolo + "-003" + ".bmp", tampon_image_cv, coin.x, coin.y, dimCCD.x, dimCCD.y);
      charge_decoupe_image_cv(holo4, cheminHolo + "-004" + ".bmp", tampon_image_cv, coin.x, coin.y, dimCCD.x, dimCCD.y);

      /*
      charger_image2D(holo1,cheminHolo+"-001"+".bmp", coin,dimCCD);
      charger_image2D(holo2,cheminHolo+"-002"+".bmp", coin,dimCCD);
      charger_image2D(holo3,cheminHolo+"-003"+".bmp", coin,dimCCD);
      charger_image2D(holo4,cheminHolo+"-004"+".bmp", coin,dimCCD);
      */

      ///------------------------___CALCUL CHAMP COMPLEXE PAR IDP___---------------------------------------------
      for(int pixel=0; pixel<NbPixROI2d; pixel++) {
	ampli_ref[pixel]=1;//sqrt((double)mod_ref[pixel]);
	plan[pixel].Re=((double)holo1[pixel]-(double)holo3[pixel])*(masque[pixel]/(double)ampli_ref[pixel]);
	plan[pixel].Im=((double)holo4[pixel]-(double)holo2[pixel])*(masque[pixel]/(double)ampli_ref[pixel]);
      }
      //SAV_Re(plan, NbPixROI2d, "/home/mat/tomo_test/plan_reel.bin", FLOAT,"a+b");
      /*///qualité du décalage phase
	double *qualiteIDP=new double[NbPixROI2d], moindreCarreIDP=0;
	for(int pixel=0;pixel<(DIMX_CCD2*DIMY_CCD2);pixel++)
	{
	qualiteIDP[pixel]=(double)holo4Re[pixel]-(double)holo3Re[pixel]-(double)holo1Re[pixel]+(double)holo2Re[pixel];
	moindreCarreIDP+=qualiteIDP[pixel]*qualiteIDP[pixel];
	}
	cout<<"moindre carré IDP="<<moindreCarreIDP<<endl;
	SAV(qualiteIDP, NbPixROI2d, "/home/mat/tomo_test/qualiteIDP/qualiteIDP", CHAR,"a+b");                                     */
      ///-------------------------___CIRCSHIFT, TF 2D, DECOUPE Nxmax, CIRCSHIFT___---------------------------------------------

      circshift2DCplx(plan, plan_shift,dimCCD,decalCCD);///-----Circshift avant TF2D
      //SAV_Re(plan_shift, NbPixROI2d, "/home/mat/tomo_test/plan_reel.bin", FLOAT,"a+b");
      TF2Dcplx(plan_shift,fft_tmp,dimCCD);
      Var2D dim2D= {2*NXMAX,2*NYMAX},decal2D= {NXMAX,NYMAX}, NMAX={NXMAX,NYMAX};
      decalCoupeCplx(fft, fft_tmp, NMAX, dimCCD);
      circshift2DCplx(fft, fft_shift,dim2D,decal2D);///semble bousiller fft??
      SAV_Re(fft_shift,4*NXMAX*NXMAX, "Tf2d_reel.bin", FLOAT,"a+b");

      ///-----------------------------------___CHERCHER SPÉCULAIRE=MAXIMUM___-----------------------------------
      //gloubiboulga function /!\ achtung, dont touch /!!\ Caramba, Né tochas

      int cpt_max=0;
      fft_module_shift[0]=pow(fft_shift[0].Re,2)+pow(fft_shift[0].Im,2);
      for(int cpt=1; cpt<(4*NXMAX*NYMAX); cpt++) {
	fft_module_shift[cpt]=pow(fft_shift[cpt].Re,2)+pow(fft_shift[cpt].Im,2);
	if(fft_module_shift[cpt]>fft_module_shift[cpt_max]) {
	  cpt_max=cpt;
	}
      }
      double max_part_reel = fft_shift[cpt_max].Re;///sauvegarde de la valeur cplx des  spéculaires
      //cout<<"angle "<<cpt_angle<<" ->max_part_reel = "<<max_part_reel<<endl;
      double max_part_imag = fft_shift[cpt_max].Im;
      double max_module = fft_module_shift[cpt_max];
      //SAV(fft_reel_shift, 4*NXMAX*NYMAX, "/home/mat/tomo_test/test_resultat_shift.bin", CHAR,"a+b");

      ///----------------------------Coordonnées dans l'espace 2D à partir de l'indice 1D------------
      int xmi=cpt_max%(2*NXMAX), ymi=cpt_max/(2*NYMAX);
      int xc=xmi-NXMAX, yc=ymi-NYMAX;
      int xmij=2*NXMAX-xmi, ymij=2*NYMAX-ymi;//coordonnée objet jumeaux
      int cptJumo_max=2*NYMAX*ymij+xmij;
      //cout<<"***rapport spéculaire/jumeau, angle "<<cpt_angle<<""<<endl;
      // if(abs(fft_reel_shift[cpt_max]/fft_reel_shift[cptJumo_max]) || abs(fft_imag_shift[cpt_max]/fft_imag_shift[cptJumo_max])<10)
      // cout<<"Le jumeau semble mal éliminé...R="<<fft_reel_shift[cpt_max]/fft_reel_shift[cptJumo_max]<<endl;

      ///----------------------------___VIRER JUMEAU + ORDRE ZÉRO (RÉSIDUELS)___------------------------------

      if((xc*xc+yc*yc)>9) { //35*35=1225;objet jumeau pas trop près de l'objet
	jumeau_elimine++;
	int t_mask_var=round(sqrt(xc*xc+yc*yc)/7)*2+1;
	// cout<<"tmask="<<t_mask_var<<endl;
	double* cache_jumeau3=new double[t_mask_var*t_mask_var];
	antigaussienne(cache_jumeau3,t_mask_var,round(t_mask_var*2),1,0);
	Var2D posJumeau={xmij,ymij}, posOrdreZero={NXMAX,NXMAX};
	multiplier_masque2Cplx(fft_shift, cache_jumeau3, NXMAX, t_mask_var, posJumeau);
	multiplier_masque2Cplx(fft_shift, cache_jumeau3, NXMAX, t_mask_var, posOrdreZero);
	delete[] cache_jumeau3;
      }
      //SAV(fft_reel_shift, 4*NXMAX*NYMAX, "/home/mat/tomo_test/TF2d_apres_masquage/fft_reel_shift.bin.bin", CHAR,"a+b");
      ///---------------------------___NORMALISER PAR SPÉCULAIRE=RECALER EN PHASE___--------------------------------------------------

      for(int cpt=0; cpt<(4*NXMAX*NYMAX); cpt++) {
	fft_shift_norm[cpt].Re=(fft_shift[cpt].Re*max_part_reel+fft_shift[cpt].Im*max_part_imag)/max_module;
	fft_shift_norm[cpt].Im=(fft_shift[cpt].Im*max_part_reel-fft_shift[cpt].Re*max_part_imag)/max_module;
      }
      //SAV(fft_imag_shift, 4*NXMAX*NYMAX, "/home/mat/tomo_test/test_resultat_imag_shift_norm_C.bin", CHAR,"a+b");
      ///----------------------------------___RÉTROPROJECTION 3D___------------------------------------

      int xm0=(xmi-NXMAX), ym0=(ymi-NYMAX) ;//coordonnée dans l'image2D centrée (xm0,ym0)=(0,0)=au centre de l'image

      if(xm0==0 && ym0==0)///indiquer les "blancs"
	printf("(xm0,ym0)=(0,0) dans le plan %i\n", cpt_angle);

      if((xm0*xm0+ym0*ym0)>rayon_inf) { //inutile?
	centres_exclus++;
      } else {
	centre[xmi*2*NXMAX+ymi]=cpt_angle;//sauvegarde des centres en coordonnées non centrée; on met le numero d'angle
	//pour vérifier les pb.
	Var2D posSpec={xm0,ym0}, NMAX={NXMAX,NYMAX};
	Var3D decal={decalX,decalY,decalZ};
	retroPropag(spectre3D_Re,spectre3D_Im, sup_redon, dim_final,fft_shift_norm, posSpec, decal, NMAX, rayon);
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
  delete[] plan, plan_shift;
  delete[] fft_tmp, fft, fft_shift, fft_module_shift, fft_shift_norm;
  delete[] masque;

  ///--------------------------------____AFFICHER QUELQUES INFOS IMPORTANTES

  printf("points_faux!! : %i,points_faux\n",points_faux);
  printf("jumeaux re-elimines : %i\n",jumeau_elimine);
  printf("temps_total proj : %f \n",temps_total);
  printf("nb_proj: %i\n",nb_proj);
  printf("Nombre de centre exclus : %i\n",centres_exclus);
  temps_final = clock ();
  temps_cpu = (temps_final - temps_initial) * 1e-6/nb_proj;
  printf("temps moyen pour 1 angle: %f\n",temps_cpu);
  temps_cpu = (temps_final - temps_initial) * 1e-6;
  printf("temps total pour %i angle(s): %f\n", nb_proj, temps_cpu);
  SAV(centre, 4*NXMAX*NYMAX, "OUT/centres.bin", UINT,"a+b");
  cout<<"centres exportés"<<endl;
  SAV(sup_redon, N_tab, "OUT/sup_redon_C.bin", UINT,"a+b");

  temps_initial=clock();//enclenchement du chronometre
  ///-----------------------___MOYENNAGE PAR SUP_REDON___----------------------------------------------------------
  for(int cpt=0; cpt<N_tab; cpt++) {
    if (sup_redon[cpt]==0) //remplace les 0 de sup_redon par des 1---> evite la division par 0
      sup_redon[cpt]=1;

    spectre3D_Re[cpt] = spectre3D_Re[cpt]/sup_redon[cpt];//normalisation par sup_redon
    spectre3D_Im[cpt] = spectre3D_Im[cpt]/sup_redon[cpt];
  }
  temps_final = clock ();
  temps_cpu = (temps_final - temps_initial) * 1e-6;
  printf("temps apres normalisation : %lf\n",temps_cpu);
  delete[] sup_redon;

  //SAV(spectre3D_Re, Ntab, "/home/mat/tomo_test/spectre3D_Re_norm_avant.bin", DOUBLE,"a+b");
  //interp3D(spectre3D_Re,  dimVolX,dimVolX, dimVolX);///interpolation dans Fourier
  //interp3D(spectre3D_Im,  dimVolX,dimVolX, dimVolX);

  //////////////////////////////ecriture de papillon binarisé
  double *papillon_masque=new double[N_tab];
  for(int compteur=0; compteur<N_tab; compteur++) {
    if(spectre3D_Re[compteur]==0)
      papillon_masque[compteur]=0;
    else
      papillon_masque[compteur]=1;
  }
  //SAV(papillon_masque, N_tab, "/home/mat/tomo_test/papillon_masque.bin", CHAR,"a+b");
  delete[] papillon_masque;
  ///-----------------------------------CIRCSHIFT3D, TF3D, SAUVEGARDE--------------------------------------------------

  printf("\n********circshift avant TF3D ************\n");
  printf("*******************************************\n");

  /////////////////////creation des tableaux qui recoivent le circshift
  double *spectre3D_Re_shift=new double[N_tab];//partie reelle du papillon dans l'espace reciproque
  double *spectre3D_Im_shift=new double[N_tab];//partie imaginaire du papillon dans l'espace reciproque

  temps_initial = clock();//enclenchement chronometre
  Var3D decal3D= {dimVolX/2,dimVolX/2,dimVolX/2}, dimFinal= {dimVolX,dimVolX,dimVolX};
  circshift3D2(spectre3D_Re, spectre3D_Re_shift, dimFinal,decal3D);//ancienne valeur, causant des problèmes d'arrondis
  circshift3D2(spectre3D_Im, spectre3D_Im_shift, dimFinal,decal3D);

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

  delete[] spectre3D_Re;//suppression des tableaux non circshiftés
  delete[] spectre3D_Im;

  ///---------------------------TF3D---------------------------------
  //TF3D(spectre3D_Re_shift,spectre3D_Im_shift,res_reel, res_imag);

  fftw_plan_with_nthreads(nthreads);
  int N3D=N_tab;
  //Déclaration des variables pour la FFT : entre,sortie et "fftplan"
  fftw_complex *in3D, *out3D;
  fftw_plan p3D;
  //Réservation memoire
  in3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3D);
  //out3D = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N3D);
  // p = fftw_plan_dft_2d(nx,ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  ///Récupération de l'image dans la partie reelle de l'entree
  for(int cpt=0; cpt<(N3D); cpt++) {
    in3D[cpt][0]=spectre3D_Re_shift[cpt];
    in3D[cpt][1]=spectre3D_Im_shift[cpt];
  }

  delete[] spectre3D_Re_shift,  spectre3D_Im_shift;//libertation de la memoire d'entree
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

  ///-------------------------------------CIRCSHIFT3D APRÈS TF3D-->OBJET 3D-----------------------

  double *final_reel_shift=new double[N_tab];//partie reelle du resultat final shifté (centré)
  double *final_imag_shift=new double[N_tab];//partie imaginaire du resultat final shifté (centré)

  circshift3D2(final_reel, final_reel_shift, dimFinal, decal3D);
  circshift3D2(final_imag, final_imag_shift, dimFinal, decal3D);

  delete[] final_reel, final_imag;

  ///---------------calcul du module carre
  double *valeur_module_shift=new double[N3D];

  for(int cpt=0; cpt<(N3D); cpt++) {
    valeur_module_shift[cpt] = final_reel_shift[cpt]*final_reel_shift[cpt]+final_imag_shift[cpt]*final_imag_shift[cpt];//calcul du module
  }
  temps_final = clock ();
  temps_cpu = (temps_final - temps_initial) * 1e-6;
  printf("temps apres circshift 3D final et calcul du module: %f\n",temps_cpu);
  temps_initial = clock();//enclenchement chronometre

  printf("*******************************************\n");
  printf("ecriture des résultats\n");
  printf("*******************************************\n");

  FILE* fichier_final_reel_shift = NULL;
  FILE* fichier_final_imag_shift = NULL;
  FILE* fichier_final_modul_shift = NULL;

  fichier_final_reel_shift = fopen("OUT/final_reel_shift.bin", "wb");
  fichier_final_imag_shift = fopen("OUT/final_imag_shift.bin", "wb");
  fichier_final_modul_shift = fopen("OUT/final_modul_shift.bin", "wb");

  float precision;

  for(int cpt=0; cpt<(N3D); cpt++) {
    precision=final_reel_shift[cpt];
    fwrite(&precision,sizeof(precision),1,fichier_final_reel_shift);

    precision=final_imag_shift[cpt];
    fwrite(&precision,sizeof(precision),1,fichier_final_imag_shift);

    // precision=valeur_module_shift[cpt];
    // fwrite(&precision,sizeof(precision),1,fichier_final_modul_shift);//ecriture du module
  }
  ecrire_rapport(NXMAX,rayon,Rf,DIMX_CCD2,coin.x,coin.y,sizeof(precision)*4,Dossier_acquiz,nb_proj,n1,NA,TpCam,G);

  delete[] ampli_ref, mod_ref;

  //libération memoire allouée pour les threads
  void fftw_cleanup_threads(void);

  fclose(fichier_final_reel_shift);
  fclose(fichier_final_imag_shift);
  fclose(fichier_final_modul_shift);
  printf("ecriture de final_reel_shift.bin,  final_imag_shift.bin, final_modul_shift : OK \n");
  //liberation des variables "resultat final"
  delete[] final_reel_shift, final_imag_shift;
  delete[] valeur_module_shift;

  temps_arrivee = clock ();
  temps_cpu = (temps_arrivee-temps_depart ) * 1e-6;
  printf("temps total: %f\n",temps_cpu);
  return 0;
}
