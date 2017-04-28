#ifndef __UTIL_IMAGE__
#define __UTIL_IMAGE__



/// A TERME, vImg doit *TOUT* remplacer

#include "cv.h"

/* #include <Magick++.h> */
/* using namespace Magick; */


/* extern Image g_ReadImage; */
/* extern double* g_entree_shift; */
/* extern char *g_OUTPUT_DIR; */



//Fonction transferant les pixels de l'image 2D de taille (nx,ny) dans un tableau 1D de taille (nx*ny) déjà alloué
// ARRAY_ALLOC(array, taille_x * taille_y, unsigned char);
// nécessite la bibliotheque PGM // openCV selon
void
remplit_tableau_cv(unsigned char *array, char *path, size_t taille_x, size_t taille_y);



void
charge_decoupe_image_cv(unsigned char* dst_array, const char* path, IplImage* tmpimage, size_t window_edge_x, size_t window_edge_y, size_t taille_x, size_t taille_y);



void
normalize_vectra1(unsigned char* src, double* dst, size_t N);



// =============================================================================


void circshift2(double* entree, double* result, int dimX, int dimY, int centreX, int centreY);


#endif // __UTIL_IMAGE__




// ==============================================================================
// cimetière du code peu performant
// ==============================================================================


/*
// ecrit dans le tableau array le contenu de l'image stockée dans path et de taille spécifiée
// array doit déjà être alloué à la bonne taille
void
load_2D_image(unsigned char* dst_array, int numero, char* path, int cpt_fichier, int window_edge_x,int window_edge_y,int taille_x,int taille_y);  
*/


/*

// charge le contenu de l'image PGM (path) dans un tableau alloué dst_array (taille_x * taille_y)
// l'image subit une découpe: on fournir le coin supérieur gauche (coinx/y) et la taille de la fenêtre de découpe (taille x/y)
// il faut fournir un pointeur pour un tampon 
void
charge_decoupe_image_lpgm(unsigned char* dst_array, char* path, unsigned char *tmp_array, size_t window_edge_x, size_t window_edge_y, size_t taille_x, size_t taille_y);


// soient deux images allouées en row-major
// découpe une fenêtre de img_src et l'écrit dans img_dst
// il faut préciser la taille des deux images et la position du coin supérieur gauche de la fenêtre
template <typename T>
void 
crop_window(T* img_src, size_t src_dimx, size_t src_dimy, T* img_dst, size_t dst_dimx, size_t dst_dimy, size_t ul_corner_x, size_t ul_corner_y)
{
  assert(dst_dimx < src_dimx);
  assert(dst_dimy < src_dimy);

  assert((ul_corner_x + dst_dimx) <= src_dimx);
  assert((ul_corner_y + dst_dimy) <= src_dimy);

  T *start_src = img_src + (ul_corner_y * src_dimx);

  T *curs_src = start_src;
  size_t pos_dst = 0;
  size_t i, j;

  for (i = 0; i < dst_dimy; i++)
    {
      for (j = ul_corner_x; j < ul_corner_x + dst_dimx; j++)
	{
	  img_dst[pos_dst] = curs_src[j];
	  pos_dst++;
	}
      curs_src += src_dimx;
    }
}

*/
