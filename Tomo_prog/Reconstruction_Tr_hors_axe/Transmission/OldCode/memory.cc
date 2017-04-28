#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>
// #include <Magick++.h>

#include "cv.h"

using namespace std;
using namespace cv;

#include "main.h"
#include "memory.h"



// circshift sur place
// version performante à base de memcopy, mais qui marche, elle.
void
circshift_2D_cv(double* entree, size_t dimx, size_t dimy)
{
  size_t cx = dimx / 2;
  size_t cy = dimy / 2;

  // on habille entree avec une structure cv::Mat
  Mat mat_src(dimx, dimy, CV_64F, entree);

  Mat q0(mat_src, Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant 
  Mat q1(mat_src, Rect(cx, 0, cx, cy));  // Top-Right
  Mat q2(mat_src, Rect(0, cy, cx, cy));  // Bottom-Left
  Mat q3(mat_src, Rect(cx, cy, cx, cy)); // Bottom-Right

  Mat tmp;                           // swap quadrants (Top-Left with Bottom-Right)
  q0.copyTo(tmp);
  q3.copyTo(q0);
  tmp.copyTo(q3);

  q1.copyTo(tmp);                    // swap quadrant (Top-Right with Bottom-Left)
  q2.copyTo(q1);
  tmp.copyTo(q2);
}


// circshift déporté
void
circshift_2D_cv(double* entree, double* sortie, size_t dimx, size_t dimy)
{
  size_t cx = dimx / 2;
  size_t cy = dimy / 2;

  // on habille entree avec une structure cv::Mat
  Mat mat_src(dimx, dimy, CV_64F, entree);
  Mat mat_dst(dimx, dimy, CV_64F, sortie);

  Mat qs0(mat_src, Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant 
  Mat qs1(mat_src, Rect(cx, 0, cx, cy));  // Top-Right
  Mat qs2(mat_src, Rect(0, cy, cx, cy));  // Bottom-Left
  Mat qs3(mat_src, Rect(cx, cy, cx, cy)); // Bottom-Right

  Mat qd0(mat_dst, Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant 
  Mat qd1(mat_dst, Rect(cx, 0, cx, cy));  // Top-Right
  Mat qd2(mat_dst, Rect(0, cy, cx, cy));  // Bottom-Left
  Mat qd3(mat_dst, Rect(cx, cy, cx, cy)); // Bottom-Right


  Mat tmp;                           // swap quadrants (Top-Left with Bottom-Right)
  qs0.copyTo(tmp);
  qs3.copyTo(qd0);
  tmp.copyTo(qd3);

  qs1.copyTo(tmp);                    // swap quadrant (Top-Right with Bottom-Left)
  qs2.copyTo(qd1);
  tmp.copyTo(qd2);
}



// fonctionne enfin!  
void 
circshift_2D_mc(double *image2D, double *image2D_shift, int dim_x, int dim_y)
{
  int cx = floor( dim_x / 2 );
  int cy = floor( dim_y / 2 );

  // Creating a local structure with member functions enables us to capture an environment and then to invoke functions with a reduced set of parameters, alike lisp lambda functions
  struct Local{
    double* tab;
    int dimx;
    int dimy; 
    
    // returns address of array element     
    inline double* 
    T(int i, int j)
    {return (tab + j + i * dimx);}
    
 
    // source corner, destination corner, on same array or not. 
    // arrays are supposed to be of same dimension as Local, or smaller
    inline void
    SQUARECPY(double *CS, double *CD)
    {
      int j;
      const int cly = floor( dimy / 2);
      const int nbytes = cly * sizeof(double);
      int offset = 0;
      
      #pragma omp parallel for num_threads(4) schedule(static) private(offset)
      for (j = 0; j < cly; j++)
	{
	  offset = j * dimy;
	  memcpy(CD + offset, CS + offset, nbytes);
	  //avant: offset += dimy;
	}
    }
  };
 
   
  // declare lambda functions for source and destination cubes
  struct Local ts{image2D, dim_x, dim_y};
  struct Local td{image2D_shift, dim_x, dim_y};
  
  // pre-compute sub-cubes upperleft corners  
  //   // source cube  
  double* S1 = ts.T(0, 0);
  double* S2 = ts.T(cx, 0);
  double* S3 = ts.T(0, cy);
  double* S4 = ts.T(cx, cy);
  
  double* D1 = td.T(0, 0);
  double* D2 = td.T(cx, 0);
  double* D3 = td.T(0, cy);
  double* D4 = td.T(cx, cy);
  

  // performs block-copies  // omp ne sert à rine ici
  ts.SQUARECPY(S1, D4);
  ts.SQUARECPY(S2, D3); 
  ts.SQUARECPY(S3, D2); 
  ts.SQUARECPY(S4, D1); 
    
}







// fonctionne très bien, mais surtout prévu pour un décalage arbitraire. Un poil trop lent
void
circshift_ccd2(double* entree, double* entree_shifted, int dimx, int dimy, int decal_x,int decal_y)
{
  int xi, yi, pixel, pixel_shift;
  for(xi=0;xi<decal_x;xi++)
    {
      for(yi=0;yi<decal_y;yi++)
	{
	  pixel=yi*dimx+xi;
	  pixel_shift=(yi+decal_y)*dimx+xi+decal_x;
	  //1er quadrant vers 4 eme
	  entree_shifted[pixel_shift]=entree[pixel];
	  //4 eme quadrant vers 1er
	  entree_shifted[pixel]=entree[pixel_shift];
	  //2eme vers 3eme
	  entree_shifted[(yi+decal_y)*dimx+xi]=entree[pixel+decal_x];
	  //3eme vers 2eme
	  entree_shifted[pixel+decal_x]=entree[(yi+decal_y)*dimx+xi];
	}
    }
}




// efficace, mais désormais inutilisé (tâche dévouée à FFTW_Volume ou AIR_Volume)
void 
circshift3D_memcpy(double *volume3D, double *volume3D_shift, int dim_x, int dim_y, int dim_z)
{

  // Creating a local structure with member functions enables us to capture an environment and then to invoke functions with a reduced set of parameters, alike lisp lambda functions
  struct Local{
    double* tab;
    int dimx;
    int dimy; 
    int dimz;
    
    // returns adress of array element     
    inline double* 
    T(int i, int j, int k)
    {return (tab + k + dimz * (j + i * dimy));}
    
    // return the adress shift from array beginning to array element
    inline int
    T_shift(int i, int j, int k)
    {return (k + dimz * (j + i * dimy));}

    // source corner, destination corner, on same array or not. 
    // arrays are supposed to be of same dimension as Local, or smaller
    inline void
    CUBECPY(double *CS, double *CD)
    {
      int i, j;
      int clx = floor(dimx / 2);
      int cly = floor(dimy / 2);
      int clz = floor(dimz / 2);
      int nbytes = clz * sizeof(double);
      int shift;
      for (i = 0; i < clx; i++)
	for (j = 0; j < cly; j++)
	  {
	    shift = T_shift(i, j, 0);
	    memcpy(CD + shift, CS + shift, nbytes);
	  }
    }
  };
 
   
  // declare lambda functions for source and destination cubes
  struct Local ts{volume3D, dim_x, dim_y, dim_z};
  struct Local td{volume3D_shift, dim_x, dim_y, dim_z};
  
  // pre-compute sub-cubes upperleft corners  
  //   // source cube  
  double* S1 = ts.T(0, 0, 0);
  double* S2 = ts.T(0, floor(dim_y / 2), 0);
  double* S3 = ts.T(floor(dim_x / 2), 0, 0);
  double* S4 = ts.T(floor(dim_x / 2), floor(dim_y / 2), 0);
  double* S5 = ts.T(0, 0, floor(dim_z / 2));
  double* S6 = ts.T(0, floor(dim_y / 2), floor(dim_z / 2));
  double* S7 = ts.T(floor(dim_x / 2), 0, floor(dim_z / 2));
  double* S8 = ts.T(floor(dim_x / 2), floor(dim_y / 2), floor(dim_z / 2));
  //   // dest cube  
  double* D1 = td.T(0, 0, 0);
  double* D2 = td.T(0, floor(dim_y / 2), 0);
  double* D3 = td.T(floor(dim_x / 2), 0, 0);
  double* D4 = td.T(floor(dim_x / 2), floor(dim_y / 2), 0);
  double* D5 = td.T(0, 0, floor(dim_z / 2));
  double* D6 = td.T(0, floor(dim_y / 2), floor(dim_z / 2));
  double* D7 = td.T(floor(dim_x / 2), 0, floor(dim_z / 2));
  double* D8 = td.T(floor(dim_x / 2), floor(dim_y / 2), floor(dim_z / 2));


  // performs block-copies  
  ts.CUBECPY(S1, D8);
  ts.CUBECPY(S8, D1);
  ts.CUBECPY(S2, D7);
  ts.CUBECPY(S7, D2);
  ts.CUBECPY(S3, D6);
  ts.CUBECPY(S6, D3);
  ts.CUBECPY(S4, D5);
  ts.CUBECPY(S5, D4);
}



void 
zerofill_double(double *array, int card)
{
  for (int i = 0; i < card; i++)
    array[i] = 0.0f;
}


void
convert_double(INPUT_TYPE* src, double* dst, int cardinal)
{
  for (int i = 0; i < cardinal; i++)
    dst[i] = double( src[i] );
}

void
convert_float(double* src, OUTPUT_TYPE* dst, int cardinal)
{
  for (int i = 0; i < cardinal; i++)
    dst[i] = OUTPUT_TYPE( src[i] );
}


int
nonzero_p(INPUT_TYPE* src, int cardinal)
{
  int count = 0;
  for (int i = 0; i < cardinal; i++)
    if (src[i]) count++;

  return count;
}


int
nonzero_d(double* src, int cardinal)
{
  int c = 0, i = 0;

  for (i = 0; i < cardinal; i++)
    if (src[i] != 0) c++;
  
  return c;
}













// =============================================================================
// kode sematary
// =============================================================================



// effectue un circshift centré de l'image donnée de dimension dim_x, dim_y



/*
  void
  extract_subImage(double* src, double* dst, int src_dimx, int src_dimy, int edge_x, int edge_y, int dst_dimx, int dst_dimy)
  {
  double* readPos = src + edge_y * src_dimx;
  double* writePos = dst;
  int nbytes = dst_dimx * sizeof(double);

  for (int l = 0; l < dst_dimy; l++)
  {
  memcpy(writePos, readPos + edge_x, nbytes);
  readPos += src_dimx;
  writePos += dst_dimx;
  }
  }
*/



