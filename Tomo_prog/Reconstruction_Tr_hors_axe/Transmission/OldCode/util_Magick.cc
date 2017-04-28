#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>

#include <Magick++.h>


using namespace std;
using namespace Magick;



#include "main.h"
#include "util.h"
#include "util_Magick.h"




struct mVar3D{ 
  int x; 
  int y; 
  int z;
};

typedef struct mVar3D Var3D;





unsigned char* 
charger_image2D(unsigned char* phasei,int numero,char *Chemin,int cpt_fichier,int coin_x,int coin_y,int taille_x,int taille_y)
{
  char CopieChemin[100];
  char FinNom[15];

  /*FILE* fichier_TF_real2 = NULL;*///nom bidon pour tester charger_image2D

  //sauvegarder chemin car sprintf écrase
  strcpy(CopieChemin,Chemin);
  sprintf(FinNom,"%ld-00%ld." INPUT_FORMAT, (long int) cpt_fichier, (long int) numero);
  //créer le chemin final
  strcat(CopieChemin, FinNom);
  phasei = rempli_tableau(CopieChemin, coin_x, coin_y,taille_x,taille_y);
  return phasei;
}



void
my_load_2D_image(unsigned char* array, char* filename, int coin_x,int coin_y,int taille_x,int taille_y)  
{
  g_ReadImage.read(filename);
  // à virer?
  g_ReadImage.type( GrayscaleType );
  g_ReadImage.compressType( NoCompression );

  g_ReadImage.getPixels(coin_x, coin_y, taille_x, taille_y);
  g_ReadImage.writePixels(GrayQuantum, array);
}



void
load_2D_image(unsigned char* array, int numero, char* path, int cpt_fichier, int coin_x,int coin_y,int taille_x,int taille_y)  
{
  char* CopieChemin;
  int taille_suffixe = 15;
  ARRAY_ALLOC(CopieChemin, strlen(path) + taille_suffixe, char);
  char FinNom[taille_suffixe];
  //sauvegarder chemin car sprintf écrase
  strcpy(CopieChemin, path);
  sprintf(FinNom,"%ld-00%ld." INPUT_FORMAT, (long int) cpt_fichier, (long int) numero);
  //créer le chemin final
  strcat(CopieChemin, FinNom);

  g_ReadImage.read(CopieChemin);
  // à virer?
  g_ReadImage.type( GrayscaleType );
  g_ReadImage.compressType( NoCompression );

  g_ReadImage.getPixels(coin_x, coin_y, taille_x, taille_y);
  g_ReadImage.writePixels(GrayQuantum, array);
}


unsigned char *
rempli_tableau(string path, int coinx, int coiny,int taille_x, int taille_y)
{
  /////// declaration des variables

  Image Monimage;
  //int i, currentImageWidth, currentImageHeight;
  unsigned char *finalArray;

  ////// chargement en memoire de l'image

  Monimage.read(path);


  ////// Mon image est N&B
  Monimage.type( GrayscaleType );
  Monimage.compressType(NoCompression);
  /*try{
    finalArray = new unsigned char[taille_x*taille_y];
    }
    catch (bad_alloc ex)
    {cerr << endl << ex.what(); exit(EXIT_FAILURE);}*/
  ARRAY_NEW(finalArray, taille_x * taille_y, unsigned char);
  ////// lecture de l'image
  Monimage.getPixels(coinx,coiny,taille_x,taille_y);
  ////// ecriture de l'image dans un tableau
  Monimage.writePixels(GrayQuantum, finalArray);
  //affichage pour contrôle
  return finalArray;
  // ** // delete[] finalArray;
}






