#ifndef __UTIL_MAGICK__
#define __UTIL_MAGICK__

#include <Magick++.h>

using namespace Magick;



extern Image g_ReadImage;
extern double* g_entree_shift;
extern char *g_OUTPUT_DIR;



//Fonction transferant les pixels de l'image 2D de taille (nx,ny) dans un tableau 1D de taille (nx*ny)
// nécessite la bibliotheque imageMagick

unsigned char *
rempli_tableau(string path, int coin_x, int coin_y,int taille_x,int taille_y);


//fonction chargeant une image 2D. Retourne un pointeur sur le tableau 1D qui la contient
//  Nécessite le nom du pointeur à remplir,
//le numéro du pas de phase shifting (1,2 3 ou 4), le chemin, et cpt_fichier
//pour boucler sur les 1000 projections

unsigned char* charger_image2D(unsigned char* phasei, int numero, char* chemin, int cpt_fichier, int coin_x, int coin_y, int taille_x, int taille_y);





// ecrit dans le tableau array le contenu de l'image stockée dans path et de taille spécifiée
// array doit déjà être alloué à la bonne taille
void
load_2D_image(unsigned char* array, int numero, char* path, int cpt_fichier, int coin_x,int coin_y,int taille_x,int taille_y);  




#endif
