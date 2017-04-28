#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>


#include "macros.h"
#include "Point.h"
#include "Buffer_vector.h"


#define TYPE int
//     Point<int>


using namespace std;


int 
main(void)
{
  Buffer_vector< TYPE > vec;
  int i;
  //   TYPE A(1, 2, 3);
  //   TYPE B(4, 5, 6);
  TYPE A = 200;
  TYPE B = 300;

  vec.resize_clusters(5);
  
  for (i = 0; i < 17; i++)
    {
       vec.add_element( i + 100 );
       //vec.add_element( TYPE(i + 100, 0, 0) );
       cout << endl << vec.ask_element(i);
    }


  cout << endl << endl ;
  for (i = 0; i < vec.get_elements(); i++)
    cout << endl << vec.ask_element(i);
  
  cout << endl << vec.get_elements() << " " << vec.get_capacity();


  Buffer_vector< TYPE > vec2(5);
  TYPE* array = vec2.grant_first_cluster();
  
  int size;
  for (i = 0, size = vec2.get_cluster_size(); i < size; i++)
    array[i] = i;
  
  vec2.add_element(123);
  vec2[2] = 22;
  cout << endl << vec2;

  TYPE swp;
  ARRAY_REVERSE(array, size, swp); 
  cout << endl << vec2;



  //seize buffer
  Buffer_vector< TYPE > vecs(5);
  cin >> vecs;
  cout << endl << vecs;


  Buffer_vector< TYPE > vecf(5);  
  ifstream infile ("buf.txt");
  if (infile.is_open())
    {
      infile >> vecf;
    }
  cout << endl << vecf;

  return EXIT_SUCCESS;
}
