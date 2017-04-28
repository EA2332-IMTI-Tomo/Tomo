#include <iostream>
using namespace std;

#include "macros.h"


#include "FFTW_Image.h"
#include "FFTW_Image.cc"


#include "cl1.h"
#include "cl2.h"




int
main(void)
{
 
  double* plan_reel;
  ARRAY_ALLOC(plan_reel, 256 * 256, double);
  double* plan_imag;
  ARRAY_ALLOC(plan_imag, 256 * 256, double);

  
  FFTW_Image<double, fftw_complex> *VF_Fourier2D;  
  VF_Fourier2D = new FFTW_Image<double, fftw_complex>(256, 256, 3);
  VF_Fourier2D -> allocate();
  
  VF_Fourier2D -> import_from(plan_reel, plan_imag, true);
  VF_Fourier2D -> set_fourier_forward();
  VF_Fourier2D -> export_to(plan_reel, plan_imag, true); 
  
  VF_Fourier2D -> unallocate();


  cl1 i1(256, 256);
  i1.do1();

  cl2 i2(512, 512);
  i2.do2();
  

  return 0;
}
