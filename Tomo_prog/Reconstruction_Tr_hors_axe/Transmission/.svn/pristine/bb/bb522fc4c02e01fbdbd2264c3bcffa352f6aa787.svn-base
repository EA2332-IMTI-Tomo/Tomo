#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <macros.h>

#include <iostream>

using namespace std;


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
      int cly = floor( dimy / 2);
      int nbytes = cly * sizeof(double);
      int offset = 0;
      
      for (j = 0; j < cly; j++)
	{
	  memcpy(CD + offset, CS + offset, nbytes);
	  offset += dimy;
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
  

  // performs block-copies  
  ts.SQUARECPY(S1, D4);
  ts.SQUARECPY(S2, D3);
  ts.SQUARECPY(S3, D2);
  ts.SQUARECPY(S4, D1);
}




int
main(void)
{
  size_t wx = 8, wy = 8;

  ARRAY_DEC_ALLOC(data, wx * wy, double);
  ARRAY_DEC_ALLOC(data_s, wx * wy, double);

  size_t i, j;
  for (i = 0; i < wx * wy; i++)
    data[i] = i;


  cout << endl;
  
  for (j = 0; j < wy; j++)
    {
      for (i = 0; i < wx; i++)
	{
	  cout << data[j * wx + i] << " ";
	}
      cout << endl;
    }

  cout << endl;
  cout << endl;

  
  circshift_2D_mc(data, data_s, wx, wy);


  
  for (j = 0; j < wy; j++)
    {
      for (i = 0; i < wx; i++)
	{
	  cout << data_s[j * wx + i] << " ";
	}
      cout << endl;
    }

}
