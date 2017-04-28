#ifndef _PLANE_
#define _PLANE_


#include <iostream>
#include <string>

#include "macros.h"


#include "Point.h"
#include "Vector.h"


using namespace std;


/* pre-declaration of classes to suit with following prototypes */
template <class T>
class Plane;


template <class T>
class Plane
{
 private:
  void set(const Point<T> &inPoint, const Vector<T> &normal)
    {this -> inPoint = inPoint; this -> normal = normal;}
  
 protected:
  Point<T> inPoint;
  Vector<T> normal;
  
 public:
  Plane()
    {set(Point<T>(), Vector<T>());}
  Plane(const Point<T> &inPoint, const Vector<T> &normal)
    {set(inPoint, normal);}

  void reset(const Point<T> &inPoint, const Vector<T> &normal)
    {set(inPoint, normal);}
  

  READWRITER(Point<T>, inPoint);
  READWRITER(Vector<T>, normal);
  
  ~Plane() {;}  

  friend ostream& operator<< (ostream &o, const Plane<T> &p);
  
};
  



#include "./Plane.cc"


#endif /* _PLANE_ */
