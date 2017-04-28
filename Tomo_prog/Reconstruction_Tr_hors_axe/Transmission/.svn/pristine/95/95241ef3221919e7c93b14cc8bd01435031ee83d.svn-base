#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <assert.h>

using namespace std;


#include "../macros.h"

#define TYPE double



// circshift_ccd2(plan_reel, plan_reel_shift, g_dimx_ccd, g_dimy_ccd, g_dimx_ccd/2, g_dimy_ccd/2);



void 
circshift2D_memcpy(double *image2D, double *image2D_shift, int dim_x, int dim_y)
{

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
      int cly = floor(dimy / 2);
      int nbytes = cly * sizeof(double);
      int offset = 0;
      
      for (j = 0; j < cly; j++)
	{
	  memcpy(CD + offset, CS + offset, nbytes);
	  offset += cly;
	}
    }
  };
 
   
  // declare lambda functions for source and destination cubes
  struct Local ts{image2D, dim_x, dim_y};
  struct Local td{image2D_shift, dim_x, dim_y};
  
  // pre-compute sub-cubes upperleft corners  
  //   // source cube  
  double* S1 = ts.T(0, 0);
  double* S2 = ts.T(0, floor(dim_y / 2));
  double* S3 = ts.T(floor(dim_x / 2), 0);
  double* S4 = ts.T(floor(dim_x / 2), floor(dim_y / 2));
  
  double* D1 = td.T(0, 0);
  double* D2 = td.T(0, floor(dim_y / 2));
  double* D3 = td.T(floor(dim_x / 2), 0);
  double* D4 = td.T(floor(dim_x / 2), floor(dim_y / 2));
  

  // performs block-copies  
  ts.SQUARECPY(S1, D4);
  ts.SQUARECPY(S2, D3);
  ts.SQUARECPY(S3, D2);
  ts.SQUARECPY(S4, D1);
}




int main(int argc, char** argv)
{
  assert(argc >= 2);

  size_t cubeMiArete = atoi(argv[1]);
  bool display = argc > 2;
  size_t cubeArete = 2 * cubeMiArete;
  size_t size = pow((cubeMiArete * 2), 2);
  TYPE *src, *dst;
  int i;
  
  fprintf(stdout, "\n cube length: %d | array lenght: %d", cubeArete, size);
  ARRAY_ALLOC(src, size, TYPE);
  ARRAY_ALLOC(dst, size, TYPE);
  for (i = 0; i < size; i++)
    src[i] = (TYPE)i;
  

  circshift2D_memcpy(src, dst, cubeArete, cubeArete);

  if (display)
    {
      fprintf(stdout, "\n");
      
      //MATRIX_DISPLAY(src, cubeArete, cubeArete, d);
      for (i = 0; i < size; i++)
       	fprintf(stdout, "%d|", dst[i]);
      fprintf(stdout, "\n");

      //       fprintf(stdout, "%.0f|", dst[i]);
    }

  return EXIT_SUCCESS;
}
