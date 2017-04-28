#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <assert.h>

using namespace std;


#include "../macros.h"

#define TYPE unsigned short int


void 
circshift3D_memcpy(TYPE *volume3D, TYPE *volume3D_shift, int dim_x, int dim_y, int dim_z)
{
  struct Local{
    TYPE* tab;
    int dimx;
    int dimy; 
    int dimz;
    
    inline TYPE* 
    T(int i, int j, int k)
    {return (tab + k + dimz * (j + i * dimy));}
    
    inline int
    T_shift(int i, int j, int k)
    {return (k + dimz * (j + i * dimy));}

    // source corner, destination corner, on same array or not. 
    // arrays are supposed to be of same dimension as Local, or smaller
    inline void
    CUBECPY(TYPE *CS, TYPE *CD)
    {
      int i, j;
      int clx = floor(dimx / 2);
      int cly = floor(dimy / 2);
      int clz = floor(dimz / 2);
      int nbytes = clz * sizeof(TYPE);
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
  TYPE* S1 = ts.T(0, 0, 0);
  TYPE* S2 = ts.T(0, floor(dim_y / 2), 0);
  TYPE* S3 = ts.T(floor(dim_x / 2), 0, 0);
  TYPE* S4 = ts.T(floor(dim_x / 2), floor(dim_y / 2), 0);
  TYPE* S5 = ts.T(0, 0, floor(dim_z / 2));
  TYPE* S6 = ts.T(0, floor(dim_y / 2), floor(dim_z / 2));
  TYPE* S7 = ts.T(floor(dim_x / 2), 0, floor(dim_z / 2));
  TYPE* S8 = ts.T(floor(dim_x / 2), floor(dim_y / 2), floor(dim_z / 2));
  //   // dest cube  
  TYPE* D1 = td.T(0, 0, 0);
  TYPE* D2 = td.T(0, floor(dim_y / 2), 0);
  TYPE* D3 = td.T(floor(dim_x / 2), 0, 0);
  TYPE* D4 = td.T(floor(dim_x / 2), floor(dim_y / 2), 0);
  TYPE* D5 = td.T(0, 0, floor(dim_z / 2));
  TYPE* D6 = td.T(0, floor(dim_y / 2), floor(dim_z / 2));
  TYPE* D7 = td.T(floor(dim_x / 2), 0, floor(dim_z / 2));
  TYPE* D8 = td.T(floor(dim_x / 2), floor(dim_y / 2), floor(dim_z / 2));


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


int main(int argc, char** argv)
{
  assert(argc >= 2);

  int cubeMiArete = atoi(argv[1]);
  bool display = argc > 2;
  int cubeArete = 2 * cubeMiArete;
  int size = pow((cubeMiArete * 2), 3);
  TYPE *src, *dst;
  int i;
  
  fprintf(stdout, "\n cube length: %d | array lenght: %d", cubeArete, size);
  ARRAY_ALLOC(src, size, TYPE);
  ARRAY_ALLOC(dst, size, TYPE);
  for (i = 0; i < size; i++)
    src[i] = i;
  

  circshift3D_memcpy(src, dst, cubeArete, cubeArete, cubeArete);

  if (display)
    {
      fprintf(stdout, "\n");
      
      //       MATRIX_DISPLAY(src, cubeArete, cubeArete, d);
      for (i = 0; i < size; i++)
       	fprintf(stdout, "%d|", dst[i]);
      fprintf(stdout, "\n");

      //       fprintf(stdout, "%.0f|", dst[i]);
    }

  return EXIT_SUCCESS;
}
