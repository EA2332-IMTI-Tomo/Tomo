#include <iostream>
using namespace std;


#include "FFTW_Image.h"
#include "test_incl2.h"

void 
betises2(FFTW_Image<double, fftw_complex> *VF)
{

  double* plan_reel;
  ARRAY_ALLOC(plan_reel, 256 * 256, double);



  //VF -> import_from(plan_reel, plan_reel, true);
  VF -> set_fourier_forward();
  VF -> set_fourier_backward();
  VF -> export_to(plan_reel, plan_reel, true); 
  
}
