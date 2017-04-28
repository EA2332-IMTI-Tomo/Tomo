#ifndef _avSEGMENT_
#define _avSEGMENT_


/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2001 the GLOS development team (http://glos.sourceforge.net) */


#include <iostream>
#include <string>

#include "macros.h"
#include <assert.h> // a virer

#include "avPoint.h" // really used



using namespace std;


/* pre-declaration of classes to suit with following prototypes */
template <typename T>
class avSegment;

template <typename T> ostream& operator<< (ostream &o, const avSegment<T> &s);



// // !! check vtkHoleFill pour plein d'opérateurs


template <class T>
class avSegment
{
 private:
  
  double t_inter;
  BOOL t_inter_set;
  BOOL parallel;
  
  
  void cons() 
    {
      verbose_mode = FALSE; t_inter_set = FALSE;
      t_inter = 0; parallel = FALSE;
    }
  

 protected:
  avPoint<T> A;  //first pt
  avPoint<T> B;  //last pt
  avVector<T> vector;
  
  BOOL verbose_mode;
  BOOL strict_mode; // if true (default), the segmen
  

 public:
    avSegment() {cons();}
    avSegment(const avPoint<T> &A, const avPoint<T> &B)
      {this -> A = A; this -> B = B; vector = avVector<T>(B - A); cons();}
    avSegment(const avSegment<T> &S) //copy const.
      {this -> A = S.A; this -> B = S.B; this -> vector = S.vector; cons();}

    READER(avPoint<T>, A);
    READER(avPoint<T>, B);
    
    void set(const avPoint<T> &A, const avPoint<T> &B)
      {this -> A = A; this -> B = B; vector = avVector<T>(B - A); cons();}
  
    //todo: point and vector
      
    ~avSegment() {;}  


    //returns a segment of rounded coordinates
    avSegment<int> round_coords()
      {return avSegment<int>(A.round_coords(), B.round_coords());}	

    //length of a segment
    double length() const
      {return this -> vector.norm();}
	

    // return (create & valcpy) the point at middle of segment
    avPoint<T> middle() const;
    void set_middle(avPoint<T> &dest) const;
    // return (create & valcpy) the point at parametric value (t=0 -> A; t=1 -> B)
    avPoint<T> parametric(double t) const;

    // determines if segment and plane intersects. If not parallel,
    // value of t at wich intersection might occur is stored for
    // further computations
    BOOL intersects_plane_p(const Plane<T> &pl);
    
    // return the point at which current segment *intersects* given
    // plane: assert occur if no intersection exists. Will use
    // previous computed t value if available
    avPoint<T> intersection_plane(const Plane<T> &pl);
    
    /* due to internal use of private variables, above functions cannot be const */
    
    void verbose_mode_set(BOOL value)
      {verbose_mode = value;}

    //friend ostream& operator<< <> (ostream &o, const Segment<T> &s); // nvcc cale
};


#include "./avSegment.cc"


#endif /* _avSEGMENT_ */
