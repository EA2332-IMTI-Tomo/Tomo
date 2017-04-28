#include <iostream>
#include <string>
#include <stdio.h>


#include "macros.h"


#include "Point.h"



using namespace std;



int 
main(void)
{
  int a = 5, b = 2;
  cout << a << b;

  Point<double> pt1(1.5, 1.2, 3.4);


  Point<int> toto(pt1);
  cout << "toto en int:" << toto;

  Point<double> pt2(-2.0, -2.0, 3.0);

  cout << endl << "!equals?" << pt1 << pt2 << (pt1 != pt2) << endl;
  cout << endl << "equals?" << pt1 << pt1 << (pt1 == pt1) << endl;
  cout << endl << pt1 << - pt1;


  return EXIT_SUCCESS;
}
