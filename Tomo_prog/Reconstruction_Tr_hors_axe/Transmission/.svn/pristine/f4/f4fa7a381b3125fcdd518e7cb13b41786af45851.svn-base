#ifndef __MAIN_H__
#define __MAIN_H__


// macros jonathan
#include "macros.h"



/* -------------------------------------------------------------------------- */
// 
/* -------------------------------------------------------------------------- */



// liste des types de reconstruction possibles
enum ReconType {
  _offline = 0,                     // offline classique (on lit des fichiers sur disque en une fois)
  _offline_batch = 1,               // offline par paquets (déprécié)
  _online_batch = 2,                // online (gère une attente de fichiers, ok pour connexion microscope)
  _offline_movie = 3,               // on crée juste un film sur une coupe reconstruite à plusieurs étapes
  _offline_zoom = 4                 // comme offline classique, mais avec un zoom *2 centré dans fourier 
};



// liste des types de plans de coupe possibles
enum SlicePlanType {
  _xy = 0, 
  _xz = 1, 
  _yz = 2, 
};


typedef struct d_movieParams{
  size_t movie_on_slice;
  size_t movie_angle_gap;
  SlicePlanType movie_cut_view;
} movieParams;


// name given to output files inside the results folder
//#define OUTPUT_REAL_FILENAME "/final_reel_shift.bin"
//#define OUTPUT_IMAG_FILENAME "/final_imag_shift.bin"
//#define OUTPUT_REAL_FILENAME "Diatomee_R"
//#define OUTPUT_IMAG_FILENAME "Diatomee_I"

#define OUTPUT_MOD_FILENAME "/final_modul_shift.bin"
#define OUTPUT_BUTTERFLY_FILENAME "/papillon_masque.bin"
//#define OUTPUT_REDON_FILENAME "/sup_redon_C.bin"


#define OUTPUT_CENTER_FILENAME "/centre.bin"
#define DEBUG_REEL_FILENAME "/DEBUG_REAL_ARC.bin"
#define DEBUG_IMAG_FILENAME "/DEBUG_IMAG_ARC.bin"

#define RECON_TYPE float
// #define RECON_TYPE double 
#define OUTPUT_TYPE float
//#define SUPREDON_OUTPUT_TYPE unsigned short int
#define BYTE unsigned char


#define INPUT_TYPE float

#define INPUT_FORMAT "pgm"

// debug: if input is made of pictures, build a volume with Angle1real/Angle1imag pictures in sequence
#define SAVE_INPUT_VOLUME 0


// debug: sauve les volume (double) issus de la normalisation de sup-redon, prêts à paser à la TF3D (et avant le circshift)
//#define SAVE_FFT_VOLUME 1 
// un-definir pour non


#endif
