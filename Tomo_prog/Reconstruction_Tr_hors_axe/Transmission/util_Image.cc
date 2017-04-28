#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fstream>

#include "macros.h"
#include "pgmcode.h"
#include "util_Image.h"

#include "vectra.h"

#include "cv.h"
#include "highgui.h"


using namespace std;



#define CV_LOAD_IMAGE_COLOR       1
#define CV_LOAD_IMAGE_GRAYSCALE   0
#define CV_LOAD_IMAGE_UNCHANGED  -1


// SE REFERER DE PREFERENCE A vImg.hpp


void
remplit_tableau_cv(unsigned char *dst_array, char *path, size_t taille_x, size_t taille_y)
{
  ASSERT(dst_array);
   if (! vectra::file_exists_p(path))
    {
      cerr << endl << path << "not found";
      MSG_ASSERT( false , "remplit_tableau_cv: ne trouve pas le fichier source" );
    }

  IplImage* src = cvLoadImage( path, CV_LOAD_IMAGE_GRAYSCALE );    
  ASSERT((taille_x == src -> width) && (taille_y == src -> height));

  size_t size = taille_x * taille_y;
  memset(dst_array, 0, size * sizeof(uchar));

  uchar *data = ( uchar* ) src -> imageData; 
  memcpy(dst_array, data, size * sizeof(uchar));

  cvReleaseImage( &src );
}


// fonction locale
void
extract_subImage(unsigned char* src, unsigned char* dst, size_t src_dimx, size_t src_dimy, size_t edge_x, size_t edge_y, size_t dst_dimx, size_t dst_dimy)
{
  MSG_ASSERT( (dst_dimx + edge_x <= src_dimx) && (dst_dimy + edge_y <= src_dimy),
	      "découpe hors de l'image" );
  
  unsigned char* readPos = src + edge_y * src_dimx;
  unsigned char* writePos = dst;
  size_t nbytes = dst_dimx * sizeof(unsigned char);
 
  for (size_t l = 0; l < dst_dimy; l++)
    {
      memcpy(writePos, readPos + edge_x, nbytes);
      readPos += src_dimx;
      writePos += dst_dimx;
    }
}


void
charge_decoupe_image_cv(unsigned char* dst_array, const char* path, IplImage* tmpimage, size_t window_edge_x, size_t window_edge_y, size_t taille_x, size_t taille_y)
{
  IplImage* src = cvLoadImage( path, CV_LOAD_IMAGE_GRAYSCALE );

  if (! src) {
    cerr << endl << "error loading <" << path << ">";
    MSG_ASSERT(false, "cvLoadImage failed");
  }

  size_t size = taille_x * taille_y;
  
  uchar *data = ( uchar* ) src -> imageData; 
  uchar *datacut = ( uchar* ) tmpimage -> imageData; 

  extract_subImage(data, datacut, src -> width, src -> height, window_edge_x, window_edge_y, taille_x, taille_y);

  memcpy(dst_array, datacut, size * sizeof(uchar));

  cvReleaseImage( &src );
}


/*
  cvSetImageROI(src, cvRect(window_edge_x, window_edge_y, taille_x, taille_y));
  
  cvCopy(src, tmpimage);
  uchar *data = ( uchar* )tmpimage->imageData; 
  
  size_t size = taille_x * taille_y;
  memcpy(dst_array, data, size * sizeof(uchar));
*/


void
normalize_vectra1(unsigned char* src, double* dst, size_t N)
{
  unsigned char min = 512, max = 0;

  for(size_t pixel = 0; pixel < N; pixel++)
    {
      unsigned char lu = src[pixel];

      if (lu > max)
	max = lu;
      else if (lu < min)
	min = lu;
    }

  double range = max - min;
  double invrange = 1.0f / range;

#pragma omp parallel for num_threads( 3 ) 
  for(size_t pixel = 0; pixel < N; pixel++)
    {
      dst[pixel] = double(src[pixel] - min) * invrange;
    }
#pragma omp barrier

}


// ==============================================================================
// circshift non centré
// ==============================================================================


void circshift2(double* entree, double* result, int dimX, int dimY, int centreX, int centreY)
{

  //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
  //  decal.y=decal.y%dim.y;
  //decal.x=decal.x%dim.x;

  int decal_x = centreX % dimX;
  int decal_y = centreY % dimY;


  for(int yi = 0; yi < dimY; yi++)
    {

      int yi_decal = yi + decal_y;

      //si dépassement des dimensions
      if (yi + decal_y > dimY - 1)
	yi_decal = - dimY + yi + decal_y;

      //si décage négatif
      if (yi + decal_y < 0)
	yi_decal = dimY +yi + decal_y;

      int hauteur_decal = yi_decal * dimX;
      int hauteur = yi * dimX;
      for (int xi = 0; xi < dimX; xi++)
	{
	  int xi_decal= xi + decal_x;

	  if (xi + decal_x > dimX - 1)
	    xi_decal = -dimX + xi + decal_x;
	  if (xi + decal_x < 0)
	    xi_decal = dimX + xi + decal_x;

	  int pos1D_init = hauteur + xi;
	  int pos1D_shift = hauteur_decal + xi_decal; //nouvelle position après décalage
	  result [pos1D_shift] = entree[pos1D_init];

	}
    }
}



// ==============================================================================
// cimetière du code peu performant
// ==============================================================================


/*


  void
  remplit_tableau_lpgm(unsigned char *array, char *path, size_t taille_x, size_t taille_y)
  {
  assert(vectra::file_exists_p(path));
  assert(array);

  long src_x, src_y, maxc_x, startx;
  PGM_Read_Header(path, &src_x, &src_y, &maxc_x, &startx);
  
  assert((taille_x == src_x) && (taille_y == src_y));
  ARRAY_SET(array, taille_x * taille_y, 0);
  PGM_8b_Get_Data_Linear(path, src_x, src_y, startx, array);
  }


  void
  charge_decoupe_image_lpgm(unsigned char* dst_array, char* path, unsigned char *tmp_array, size_t window_edge_x, size_t window_edge_y, size_t taille_x, size_t taille_y)
  {
  long src_x, src_y, maxc_x, startx;
  PGM_Read_Header(path, &src_x, &src_y, &maxc_x, &startx);
  
  assert(tmp_array);

  remplit_tableau(tmp_array, path, src_x, src_y);
  crop_window<unsigned char>(tmp_array, src_x, src_y, \
  dst_array, taille_x, taille_y, window_edge_x, window_edge_y);
  
  }


*/
