
void ecrire_rapport(int Nxmax,float rayon,float Rf, float K, int dimx_ccd2,int coin_x, int coin_y,double precision,char *chemin,int nb_proj,float n1,float NA,float Tp, int G);


// Interpolation3D : attend un volume "troué" et les taille en x,y,et z
void interp3D(double *volume_interp_3D, int taille_x,int taille_y,int taille_z);



//fonction de circshift "generale"
double *circshift(double *entree,int dimx,int dimy,int decal_x,int decal_y);


//ancienne fonction circshift3D uniquement pour decalage d'une demi dimension
void circshift3D(double *volume3D, double *volume3D_shift,int taille_x,int taille_y,int taille_z);


// circshift using memcpy, optimally fast
// cuts given volume in 8 subcubes and shifts these cubes inchanged:
// L1: 1 2  L2: 5 6     L'1: 8 7  L'2: 4 3 
//     3 4      7 8 =>       6 5       2 1
void 
circshift3D_memcpy(double *volIn, double *volOut, int dim_x, int dim_y, int dim_z);



double* 
circshift(double* entree, int dimx,int dimy,int decal_x,int decal_y)
{
  double *entree_shift;
  ARRAY_NEW(entree_shift, dimx * dimy, double);
  
  for(int xi=0;xi<decal_x;xi++)
    {
      for(int yi=0;yi<decal_y;yi++)
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
  delete[] entree_shift; //est ce utile dans une fonction?
}



void 
circshift3D(double *volume3D, double *volume3D_shift,int taille_x,int taille_y,int taille_z)
{
  //////attention il y a 8 cubes:
  ////// face superieur:  3 4     010 110
  //                      1 2     000 100

  //face inferieur :     7  8     011 111
  //		       5  6     001 101

  // 1   <-----> 8
  // 2   <-----> 7
  // 3   <-----> 6
  // 4   <-----> 5

  for(int xi=0;xi<taille_x/2;xi++)
    {
      for(int yi=0;yi<taille_y/2;yi++)
	{
	  for(int zi=0;zi<taille_z/2;zi++)
	    {

	      //cube 1 vers cube 8  //  000--->111
	      //test volume3D_shift[zi*taille_x*taille_y+yi*taille_x+xi]=1;
	      volume3D_shift[(zi+taille_x/2)*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi+taille_x/2]=volume3D[zi*taille_x*taille_y+yi*taille_x+xi];


	      //cube 8 vers cube 1
	      volume3D_shift[zi*taille_x*taille_y+yi*taille_x+xi]=volume3D[(zi+taille_x/2)*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi+taille_x/2];


	      //cube 2 vers cube 7  //  100 ---->  011
	      volume3D_shift[(zi+taille_x/2)*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi]=volume3D[zi*taille_x*taille_y+yi*taille_x+xi+taille_x/2];

	      //cube 7 vers cube 2
	      volume3D_shift[zi*taille_x*taille_y+yi*taille_x+xi+taille_x/2]=volume3D[(zi+taille_x/2)*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi];


	      //cube 3 vers cube 6 // 010 ---> 101
	      volume3D_shift[(zi+taille_x/2)*taille_x*taille_y+yi*taille_x+xi+taille_x/2]=volume3D[zi*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi];

	      //cube 6 vers cube 3
	      volume3D_shift[zi*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi]=volume3D[(zi+taille_x/2)*taille_x*taille_y+yi*taille_x+xi+taille_x/2];


	      //cube 4 vers cube 5 // 110 --->  001
	      volume3D_shift[(zi+taille_x/2)*taille_x*taille_y+yi*taille_x+xi]=volume3D[zi*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi+taille_x/2];

	      //cube 5 vers cube 4
	      volume3D_shift[zi*taille_x*taille_y+(yi+taille_y/2)*taille_x+xi+taille_x/2]=volume3D[(zi+taille_x/2)*taille_x*taille_y+yi*taille_x+xi];
	    }
	  
	}
    }
  
}





//version matthieu
void 
circshift3DM(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D)
{

  for(int xi=0;xi<dimFinal3D.x;xi++)
    {
      for(int yi=0;yi<dimFinal3D.y;yi++)
	{
	  for(int zi=0;zi<dimFinal3D.z;zi++)
	    {
	      ///variables décalés
	      short int x2=xi+decal3D.x;
	      short int y2=yi+decal3D.y;
	      short int z2=zi+decal3D.z;
	      ///compteur zi*taille_x*taille_y+yi*taille_x+xi
	      if(xi+decal3D.x>dimFinal3D.x)
		x2=(xi+decal3D.x)-dimFinal3D.x;
	      else
		{
		  if(xi+decal3D.x<0)
		    x2=dimFinal3D.x+(xi+decal3D.x);
		  //cout<<"coucou"<<endl;
		}
	      if(yi+decal3D.y>dimFinal3D.y)
		y2=(yi+decal3D.y)-dimFinal3D.y;
	      else
		{
		  if(yi+decal3D.y<0)
		    y2=dimFinal3D.y+(yi+decal3D.y);
		}
	      if(zi+decal3D.z>dimFinal3D.z)
		{
		  z2=(zi+decal3D.z)-dimFinal3D.z;
		}
	      else
		{
		  if(zi+decal3D.z<0)
		    z2=dimFinal3D.z+(zi+decal3D.z);
		}
	      int pixel=zi*dimFinal3D.x*dimFinal3D.y+yi*dimFinal3D.x+xi;
	      int pixel_decal=z2*dimFinal3D.x*dimFinal3D.y+y2*dimFinal3D.x+x2;
	      //         if(volume3D[pixel]==344444)
	      //         {
	      //           printf("xi+decal3D.x : %i,yi+decal3D.y : %i,zi+decal3D.z : %i\n",xi+decal3D.x,yi+decal3D.y,zi+decal3D.z);
	      //           printf("xi %i: ,yi : %i,zi %i: \n",xi,yi,zi);
	      //           printf("x2 %i: ,y2 : %i,z2 %i: \n",x2,y2,z2);
	      //           printf("volume3D_shift[pixel]: %f",volume3D_shift[pixel]);
	      //         }
	      volume3D_shift[pixel_decal]=volume3D[pixel];
	    }
	}
    }



}


void
circshift3D2(double *volume3D, double *volume3D_shift, int taille_x, int taille_y, int taille_z)
{
  Var3D vt{taille_x, taille_y, taille_z};
  Var3D vs{(int)(taille_x / 2), (int)(taille_y / 2), (int)(0)};
  circshift3DM(volume3D, volume3D_shift, vt, vs);  
}



void 
interp3D(double *volume_interp_3D, int taille_x,int taille_y,int taille_z)
{

  int z=0;
  // int cpt=x+y*dv0rf+z*taille_y*dv0rf;
  int z_min=taille_z;
  int z_max=0;

  // cout<< "valeur papillon : "<< volume_interp_3D[cpt] << endl;
  for (int x=0; x < taille_x; x++)//on balaye l'image, référentiel avec centre au milieu
    {
      for (int y=0; y<taille_y; y++)//on balaye l'image,
	{
	  z=0;
	  
	  while(z<taille_z) // pas de boucle for car l'indice est variable
	    {float valeur_test=volume_interp_3D[x+y*taille_x+z*taille_y*taille_x];
	      if(valeur_test!=0 && z<taille_z) //si valeur !=0, alors borne inférieure mais z_min ne peut valoir 19.
		{	
		  z_min=z; //on a trouvé z_min

		  z_max=z_min+1; //initialiser z_max au point suivant (sinon while jamais verifie)

		  while( z_max<taille_z && volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]==0) //soit trouver prochain z_max, soit fin de colonne
		    {	
		      z_max=z_max+1;
		      z=z_max;

		      if(z_max==taille_z-1)
			{
			  //printf("fin de colonne\n");
			}
		    }


		  if(z_max!=taille_z && z_max-z_min>1 && volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]!=0) //il faut au moins un trou pour interpoler
		    {
		      
		      //y=ax+b-> interpolation linéaire
		      double a=(volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]-volume_interp_3D[x+y*taille_x+z_min*taille_y*taille_x])/(z_max-z_min);
		      double b=volume_interp_3D[x+y*taille_x+z_max*taille_y*taille_x]-a*z_max;
		      for(int cpt_z=z_min+1; cpt_z<z_max; cpt_z++)

			{
			  volume_interp_3D[x+y*taille_x+cpt_z*taille_y*taille_x]=a*cpt_z+b;
			  //printf("interpolaion finie\n");
			}
		      
		      z=z_max-1; // nouveau compteur=ancienne borne sup (-1 car z++)
		      
		    }

		  // redémarrer le compteur de ligne à z_max
		}
	      z++;
	    }

	}
    }
}





void 
ecrire_rapport(int Nxmax,float rayon,float Rf, float K, int dimx_ccd2,int coin_x, int coin_y,double precision,char* chemin,int nb_proj,float n1,float NA,float Tp, int G)
{
  //////////////////////////////////////////////Ecriture du rapport de calcul
  time_t date;
  time(&date);

  FILE *fichier_rapport ;    /*  pointeur vers un type fichier     */
  /*  ouverture du fichier 'mon_fichier.txt' pour ecriture (w) en mode texte (t)  */
  fichier_rapport = fopen ( str_concat(g_OUTPUT_DIR, "/rapport_calcul.txt"), "wt") ;
  /*  en cas d'échec de l'ouverture, le pointeur est NULL: intercepter ce cas  */
  if (fichier_rapport == NULL){
    printf ("impossible de créer le fichier rapport_calcul.txt\n") ;
    exit (EXIT_FAILURE);
  }
  
  fprintf(fichier_rapport,"Date     : %s\nNxmax    : %i\n,Rayon : %f\nRf       : %f\nK : %f\ndimx_ccd : %d\ncoin_x   : %d\ncoin_y   : %d\nPrecision: %lu bits\nSession  : %s\nnb_proj  : %i\nindice n1: %f\nNA       : %f\nT_pixel  : %e\nG        : %i",ctime(&date),Nxmax,rayon,Rf,K,dimx_ccd2,coin_x,coin_y,(long unsigned int)8*sizeof(precision),chemin,nb_proj,n1,NA,Tp,G);
  /*  fermeture du fichier  */
  fclose (fichier_rapport) ;
}
