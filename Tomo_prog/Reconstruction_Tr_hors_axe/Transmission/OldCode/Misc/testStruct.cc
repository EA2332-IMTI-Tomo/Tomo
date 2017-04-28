#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>

using namespace std;


struct mVar3D{ 
  int x; 
  int y; 
  int z;
};

typedef struct mVar3D Var3D;


void 
main(int argc, char **argv)
{
  int taille_x = 1;
  int taille_y = 2;
  int taille_z = 3;

  struct mVar3D vt{taille_x, taille_y, taille_z};
}
