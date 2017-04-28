#ifndef _CL2_
#define _CL2_


#include <iostream>
using namespace std;

#include "FFTW_Image.h"


class cl2{

 private:
  int x, y;
  FFTW_Image<double, fftw_complex> *FFT;
  
  

 public:
 cl2(int a, int b): x(a), y(b) 
  { 
    FFT = new FFTW_Image<double, fftw_complex> (x, y, 3);
    FFT -> allocate();
  }
  
  void do2();
  
};


#endif // _CL2_
