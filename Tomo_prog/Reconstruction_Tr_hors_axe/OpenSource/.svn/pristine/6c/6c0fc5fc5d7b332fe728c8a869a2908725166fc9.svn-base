#include <Array.h>

#include <iostream>
#include <string>
#include <math.h>



using namespace std;


int 
main(void)
{
  
  Array<int> t(5);
  t.fill(8);

  t.get_data()[0] = 5;
  t.set(1, 6);
  cout << t[0] << endl;


  cout << t;

  t.resize(10);

  cout << endl <<  t;
  cout << endl << t.mean() << " " << t.variance();

  t.derive();
  cout << endl << t;

  int i, size = 5;
  int* tab;
  assert(tab = (int*) calloc(size, sizeof(int)));
  for (i = 0; i < size; i++) tab[i] = i + 1;
  
  Array<int> t2(tab, size);

  int* mytab = t2.get_data();
  mytab[2] = 55;
  cout << endl << t2;

  /*
  // bug realloc
    
  int i, size = 5;
  int* tab;
  assert(tab = (int*) calloc(size, sizeof(int)));
  for (i = 0; i < size; i++) tab[i] = i + 1;
  
  for (i = 0, cout << endl; i < size; i++) cout << tab[i] << " ";

  size = 10;
  assert(tab = (int*) realloc(tab, size));
  for (i = 0, cout << endl; i < size; i++) cout << tab[i] << " ";

  //1 2 3 4 5 
  //1 2 3 4 16 0 139048 0 0 0
  // /users/these/bailleul/These/OpenSource/Tests>g++ --version
  //g++ (GCC) 3.2.2
  */

  return EXIT_SUCCESS;
}
