#ifndef __PROJET__
#define __PROJET__

// résolution des images d'entrée
#define IMAGE_X 1312
#define IMAGE_Y 1082
// taille de la fenêtre utile
#define WINDOW_X 512
#define WINDOW_Y 512
#define LAST_ANGLE_NUM 10


#define INPUT_DIR "/ramdisk/session/"
#define RADIX_SESSION "session10032302-record"

#define IMG_DIMX 1280
#define IMG_DIMY 1082
#define CUT_ORIGIN_X 550
#define CUT_ORIGIN_Y 140

// ca devrait être 512, mais il semble que 256 suffise
// en fait, avec Nxmax de merde, on est réduit à 2x121 = 142
#define FICHIER_TF2D "Tf2d/Tf2d.bin"
//#define FICHIER_MAT2 "/home/mat/tomo_test/TF2d_apres_masquage/fft_reel_shift.bin"
#define FICHIER_MAT2 "OUT/fft_reel_shift.bin"
//#define FICHIER_MAT_CENTRE "/home/mat/tomo_test/centre.bin"
#define FICHIER_MAT_CENTRE "OUT/centre.bin"
//#define FICHIER_MAT_SUPRED "/home/mat/tomo_test/sup_redon_C.bin"
#define FICHIER_MAT_SUPRED "OUT/sup_redon_C.bin"
// #define FICHIER_MAT_PAP "/home/mat/papillon_masque.bin"
#define FICHIER_MAT_PAP "OUT/papillon_masque.bin"
// #define FICHIER_RESULT_R "/home/mat/tomo_test/final_reel_shift.bin"
#define FICHIER_RESULT_R "OUT/final_reel_shift.bin"
#define FICHIER_RESULT_I "OUT/final_imag_shift.bin"
#define FICHIER_RESULT_M "OUT/final_modul_shift.bin"
#define FICHIER_RAPPORT "OUT/rapport_calcul.txt"

/* --------------------------------------------------------------------------- */
// Types
/* --------------------------------------------------------------------------- */

typedef struct {
  double Re,Im;
}nbCplx;

typedef struct {
  int x,y;
}Var2D;

typedef struct {
  int x,y,z;
}Var3D;

typedef enum PRECISION {CHAR, UINT, INT, FLOAT, DOUBLE};
#endif



#include "cv.h"
#include "highgui.h"


#include "fonctions.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
// #include <Magick++.h> 
#include <fftw3.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <assert.h>
using namespace std;
//using namespace Magick; 

