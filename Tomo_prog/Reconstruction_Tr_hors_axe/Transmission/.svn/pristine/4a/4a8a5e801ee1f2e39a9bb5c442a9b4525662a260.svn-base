#ifndef _CL1_
#define _CL1_


#include <iostream>
using namespace std;

#include "cl2.h"
#include "FFTW_Image.h"


class cl1{

 private:
  int x, y;
  FFTW_Image<double, fftw_complex> *FFT;
  cl2* cl2_instance;


 public:
 cl1(int a, int b): x(a), y(b) 
  { 
    FFT = new FFTW_Image<double, fftw_complex> (x, y, 3);
    FFT -> allocate();

    
  }
  
  void do1();

};


#endif // _CL1_
