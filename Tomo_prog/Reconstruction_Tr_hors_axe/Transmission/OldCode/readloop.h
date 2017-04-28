#ifndef __READLOOP__
#define __READLOOP__


#define GENERATED_FILENAME_LEN 20
#define CACHE_JUMEAU_LEN 100

#include "main.h"
#include "AIR_Volume.h"

extern size_t g_window_dim_x;
extern size_t g_window_dim_y;
extern char *g_MASKPATH;

// for debug only
extern size_t image_dim_x;
extern size_t image_dim_y;


// version qui lit un répertoire d'images BMP. Pour chaque angle, série de 4 images (phase 1/2/3/4) dont on déduit une image réelle et une image imaginaire 
//double* reel_arc, double* imag_arc

void 
readLoop_angleImages(AIR_Volume<RECON_TYPE> *V_RealPart, AIR_Volume<RECON_TYPE> *V_ImagPart, unsigned short int* sup_redon, int* centre, \
		     const int Nxmax, const int Nymax, int Nxmax_Rf,	\
		     const int window_edge_x, const int window_edge_y,		\
		     const size_t image_dim_x, const size_t image_dim_y,  \
		     const int xm0_limite, const int ym0_limite,	\
		     float rayon, float delta_zmax,			\
		     int angle_start, int angle_count, int angle_jump,	\
		     const char* images_radix);


void 
readLoop_binaryVolume(RECON_TYPE* reel_arc, RECON_TYPE* imag_arc, unsigned short int* sup_redon, int* centre, \
	 const int Nxmax, const int Nymax, int Nxmax_Rf,                \
	 const int window_edge_x, const int window_edge_y,				\
	 const int xm0_limite, const int ym0_limite,			\
	 float rayon, float delta_zmax, int angle_count,                    \
	 int image_dim_x, int image_dim_y,                                \
	 char* Chemin);


#endif // __READLOOP__
