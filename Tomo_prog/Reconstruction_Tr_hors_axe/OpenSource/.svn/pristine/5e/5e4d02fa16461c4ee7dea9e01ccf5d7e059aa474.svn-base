#include <iostream>
#include <string>
#include <stdio.h>


#include "macros.h"


#include "Point.h"
#include "Vector.h"
#include "Plane.h"
#include "Segment.h"
#include "Triangle.h"


using namespace std;


template <class U>
void meswap(U& x, U& y)
{
  U tmp = x;
  x = y;
  y = tmp;
}


int 
main(void)
{
  int a = 5, b = 2;
  meswap(a, b);
  cout << a << b;

  Point<double> pt1(1.5, 1.2, 3.4);
  //   pt1.set_X(1.5);

  Point<int> toto(pt1);
  cout << "toto en int:" << toto;

  Point<double> pt2(-2.0, -2.0, 3.0);
  

  cout << endl << "!equals?" << pt1 << pt2 << (pt1 != pt2) << endl;
  cout << endl << "equals?" << pt1 << pt1 << (pt1 == pt1) << endl;
  cout << endl << pt1 
       << - pt1;


  Vector<double> v1(1.0, 1.0, 3.0);
  Vector<double> v2(-2.0, -2.0, 3.0);

  cout << endl << v1 << endl;

  cout << pt1 << endl;
  cout << pt1.distance_eucl(pt2) << endl;
  cout << " prod scalar:"; 
  cout << v1.product_scalar(v2) << " " << endl ;// << product_scalar(v1, v2) << endl;
  cout << v2.norm() << endl;
  
  cout << endl << v2;
  v2.normalize();

  cout << endl << v2 << v2.norm();
  


  // Plane tests
  Point<double> orig(0, 0, 0);
  Point<double> target(1, 0, -1);
  Plane<double> pl(orig, Vector<double>(0, 0, 1));
  cout << endl << target.distance_signed(pl);


  //segment
  cout << endl;

  Segment<double> S(Point<double>(-1, 0, 1), Point<double>(1, 0, -1));
  
  cout << endl << "parametric test: " << S.parametric(0.1) << endl;

  cout << endl << S.intersects_plane_p(pl)
       << " " <<  S.intersection_plane(pl); 
  
  
  
  Segment<double> S2(Point<double>(0, 0, 1), Point<double>(1, 1, 1));
  S2.verbose_mode_set(TRUE);
  cout << endl << S2.intersects_plane_p(pl) << endl; 
  //    << " " <<  S2.intersection_plane(pl); 
  cout << S2.length();

  //triangle
  Triangle<double> tr(Point<double>(-2, 0, 0), Point<double>(1, 1, 0), Point<double>(1, -1, 0));   
  Triangle<double> tr2(Point<double>(1, 0, 2), Point<double>(1, 0, -1), Point<double>(-2, 0, -1));   
  

  cout << endl << "Triangle intersection";
  cout << endl << tr2.intersection_plane(pl);
  cout << endl << tr2.intersection_plane(Plane<double>(Point<double>(0.0, 0.0, 1.0), Vector<double>(0, 0, 1)));

  cout << endl << tr2;



  Point<double> px(1, 2, 3);
  Point<double> py(1, 2, 3);
  Point<double> pz(1, 2, 3);
  py = px + py;
  Vector<double> vx(1, 2, 3);
  Vector<double> vy(1, 2, 3);
  Vector<double> vz(1, 2, 3);
  vx = vx * 3.2;
  vz = (vx * 3.2) + (vy - vz) ;
  

  pz = pz + vy;
  cout << endl << pz;

  return EXIT_SUCCESS;
}
