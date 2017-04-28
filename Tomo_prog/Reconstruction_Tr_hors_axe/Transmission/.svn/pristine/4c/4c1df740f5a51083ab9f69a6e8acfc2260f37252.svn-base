#include <math.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>

using namespace std;

#include "../macros.h"


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



int
main(int argc, char** argv)
{
  int i;
  double* array, *dst;
  ARRAY_ALLOC(array, 30, double);
  ARRAY_ALLOC(dst, 6, double);
  
  for (i = 0; i < 30; i++)
    array[i] = i + 1;
    
  for (cout << endl, i = 0; i < 30; i++) cout << " " << array[i]; 
  
  extract_subImage(array, dst, 5, 6, 2, 3, 2, 3);
  
  for (cout << endl, i = 0; i < 6; i++) cout << " " << dst[i]; 

  
  return EXIT_SUCCESS;
}
