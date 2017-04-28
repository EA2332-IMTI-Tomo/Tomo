#include <iostream>
using namespace std;


#include <FFTW_Image.hpp>
// #include "test_incl1.h"
// #include "test_incl2.h"

void 
betises(FFTW_Image<double, fftw_complex> *VF)
{

  double* plan_reel;
  ARRAY_ALLOC(plan_reel, 256 * 256, double);
  double* plan_imag;
  ARRAY_ALLOC(plan_imag, 256 * 256, double);



  VF -> import_from(plan_reel, plan_imag, true);
  VF -> set_fourier_forward();
  VF -> export_to(plan_reel, plan_imag, true); 
  

  //   betises2(VF);
  
}
