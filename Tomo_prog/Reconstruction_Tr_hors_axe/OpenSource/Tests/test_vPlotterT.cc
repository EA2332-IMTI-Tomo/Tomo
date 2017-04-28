#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <string.h>
#include <fstream>



#include <vPlotterT.h>

using namespace std;

#define TYPE float

// =============================================================================



// fonctions pour la fleur normale

TYPE
flower_x(TYPE t)
{
  float c4t10 = 10 * cos(4 * t);
  return c4t10 * cos(t);
}


TYPE
flower_y(TYPE t)
{
  float c4t10 = 10 * cos(4 * t);
  return c4t10 * sin(t);
}


  

int
main(void)
{
  vPlotterT<float> P( &flower_x, &flower_y );
  P.set_sampling_regular(400);
  P.set_range(0, 2 * M_PI);


  for (int i = 0; i < 400; i++)
    {
      P.compute_next();
      cout << endl << P.get_x() << " " << P.get_y();
    }

  return 0;
}
