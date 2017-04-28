#include <iostream>
#include "buf_vector.h"

int
main()
{

  buf_vector<int> vec(5);
  
  for (int i = 0; i < 7; i++)
    vec.push_back(i);


  buf_vector<int> toto(5);
  for (int i = 0; i < 7; i++)
    toto.my_push_back(i);
  
  

}
