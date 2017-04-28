#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>


#include "macros.h"

#include "vVector.h"


using namespace std;


int 
main(void)
{
  vVector<size_t> A1(5, 5);

  size_t i;
  for (i = 0; i < 17; i++)
    A1.s_insert(i);

  //cout << 
  A1.display(cout);
  
  
  size_t* data = A1.get_data();

  cout << endl;
  for (i = 0; i < A1.get_pos_insert(); i++)
    cout << data[i] << " ";

  cout << endl;
  return EXIT_SUCCESS;
}
