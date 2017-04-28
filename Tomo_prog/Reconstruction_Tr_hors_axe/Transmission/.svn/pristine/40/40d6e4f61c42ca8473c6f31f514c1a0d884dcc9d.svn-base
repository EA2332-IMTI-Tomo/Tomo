#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <assert.h>


using namespace std;


int 
main(void)
{

  double *reel_arc, *r2, *r3, *r4;


  try{
    reel_arc = new double[5000 * 50000000000];
    r2 = new double[5000 * 50000000000];
    r3 = new double[5000 * 50000000000];
    r4 = new double[5000 * 50000000000];  
  }
  catch(std::bad_alloc& ex)
    {cerr << endl << "memory allocation failure" << ex.what(); return EXIT_FAILURE; } 
  catch(const std::exception& ex)
    {cerr << endl << "GException" << ex.what(); return EXIT_FAILURE; } 
  catch(...)
    {cerr << endl << "other error"; return EXIT_FAILURE; } 
  

  return EXIT_SUCCESS;
}
