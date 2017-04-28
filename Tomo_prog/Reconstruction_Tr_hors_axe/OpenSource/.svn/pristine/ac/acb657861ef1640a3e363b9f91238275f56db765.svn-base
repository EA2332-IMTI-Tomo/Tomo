#ifndef _TRIANGLE_
#define _TRIANGLE_


/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2001 the GLOS development team (http://glos.sourceforge.net) */


#include <iostream>
#include <string>


#include "Point.h" // really used
#include "Plane.h"
#include "Segment.h"


using namespace std;


/* pre-declaration of classes to suit with following prototypes */
template <typename T>
class Triangle;

template <typename T> ostream& operator<< (ostream &o, const Triangle<T> &t);



template <class T>
class Triangle
{

 private:

  BOOL plane[3]; //indicates wether vertices ABC are above-TRUE/below-false given plane
  BOOL plane_set;
  

  void link()
    {
      // AHH! pourquoi on realloue?
      /*       this -> BC = Segment<T>(B, C); */
      AB.set(A, B);
      BC.set(B, C);
      AC.set(A, C);
      plane[0] = plane[1] = plane[2] = FALSE;
      plane_set = FALSE;
    }
      
 public: //protected:

  Point<T> A;
  Point<T> B;
  Point<T> C;
  
  Segment<T> AB;
  Segment<T> AC;
  Segment<T> BC;

  
 public:

  Triangle() 
    {link();}  
  
  Triangle(const Point<T> &A, const Point<T> &B, const Point<T> &C)
    {set(A, B, C);} // set call links
  
  ~Triangle() {;}

          
  void set(const Point<T> &A, const Point<T> &B, const Point<T> &C)
    {this -> A = A; this -> B = B; this -> C = C; link();}

  void reset(const Point<T> &A, const Point<T> &B, const Point<T> &C)
    {set(A, B, C);}
  

  //determines if current triangle intersects given plane. Result
  //undefined if parallel. Sets private values for avoiding further
  //recomputations.
  BOOL intersects_plane_p(const Plane<T> &pl);
  
  //determines the intersection of current triangle with given
  //plane. Asserts if no intersection.
  Segment<T> intersection_plane(const Plane<T> &pl);

  //returns a barycenter point
  Point<T> barycenter() const;
  //same without allocation
  void set_barycenter(Point<T> &b) const;
  
  // hauteur liant A à BC (projeté de A sur BC)
  double height() const; 
  double area() const;

  friend ostream& operator<< <> (ostream &o, const Triangle<T> &t);
      
};


#include "./Triangle.cc"


#endif /* _TRIANGLE_ */
