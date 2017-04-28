#include <stdlib.h>
#include <stdio.h>
#include <iostream>
using namespace std;


#include "macros.h"
#include "vArray.h"


int
main(int argc, char** argv)
{
  int size = 10;
  int* tab;

  if (argc - 1)
    size = atoi( argv[1] );
  
  ARRAY_ALLOC(tab, size, int);


  vArray<int> TAB(tab, size);

  TAB.fill(5);

 
  TAB[2] = 10;
  TAB[7] = 4;
  
  cout << TAB << endl;
  cout << TAB.mean() << endl;
  cout << TAB.sum() << endl;
  cout << TAB.contains(4) << endl;
  cout << TAB.contains(3) << endl;
  cout << TAB.max() << endl;

  int minval, maxval;
  TAB.get_min_max(minval, maxval);
  cout << "min=" << minval << " max=" << maxval << endl;
  cout << endl << TAB.position_min();

  TAB.Free();
  cout << endl << endl;


 
  vArray<float> TAB2(10);
  TAB2.fill(5.5f);

  float* tab2_data = TAB2.get_data();
  tab2_data[4] = 4.44;
 
  cout << TAB2 << endl;

  return EXIT_SUCCESS;
}
