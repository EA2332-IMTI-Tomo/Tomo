#include <iostream>
using namespace std;

#include "cl1.h"



void 
cl1::do1()
{
  cout << "bonjour de cl1";

  
  double *a1, *a2;
  ARRAY_ALLOC(a1, x * y, double);  
  
  FFT -> import_from(a1, a1, false);
  FFT -> set_fourier_backward();
  FFT -> export_to(a1, a1, false);

  
  cl2_instance = new cl2(x, y);
  cl2_instance -> do2();
}
