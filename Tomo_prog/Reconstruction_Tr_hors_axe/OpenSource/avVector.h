#ifndef _avVECTOR_
#define _avVECTOR_


/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2001 the GLOS development team (http://glos.sourceforge.net) */


#include <iostream>
#include <string>

#include "macros.h"




using namespace std;


/* pre-declaration of classes to suit with following prototypes */
template <typename T>
class avVector;


#include "avPoint.h" // really used



/* friend functions prototypes */
/* template <typename T> T product_scalar(const avVector<T> &v1, const avVector<T> &v2); */
template <typename T> avVector<T> operator+(const avVector<T> &v1, const avVector<T> &v2);
template <typename T> avVector<T> operator-(const avVector<T> &v1, const avVector<T> &v2);
template <typename T> avVector<T> operator-(const avVector<T> &v);
template <typename T> bool operator==(const avVector<T> &v1, const avVector<T> &v2);
template <typename T> bool operator!=(const avVector<T> &v1, const avVector<T> &v2);
template <typename T> ostream& operator<<  (ostream &o, const avVector<T> &v);


template <class T>
class avVector: public avPoint<T>
{
 private:  
  /* should not be necessary since present in mother 
  void set(T x, T y, T z)
    {avVector<T>::X = x; avVector<T>::Y = y; avVector<T>::Z = z;}
  void set(const avVector &v) 
    {set(v.X, v.Y, v.Z);}
  */
  
 protected:
  bool verbose_mode;

 public:
  avVector() 
    {this->set(0, 0, 0); this->is_valid = TRUE;}
  avVector(T x, T y, T z)
    {this -> set(x, y, z); this->is_valid = TRUE;}
  avVector(const avPoint<T> &p)
    {this -> set(p.get_X(), p.get_Y(), p.get_Z());}
  avVector(const avPoint<T> &p1, const avPoint<T> &p2)
    {this -> set(p2.get_X() - p1.get_X(), p2.get_Y() - p1.get_Y(), p2.get_Z() - p1.get_Z());
    this->is_valid = TRUE;}
  avVector(const avVector<T> &v) //copy constructor
    {this -> set(v.get_X(), v.get_Y(), v.get_Z()); this->is_valid = TRUE;}

  
  avPoint<T> toavPoint() const
    {return avPoint<T>(this -> X, this -> Y, this -> Z);}
  
  ~avVector() {;}


  // scalar product (dot product)
  T product_scalar(const avVector<T> &v) const;
  avVector<T> product_vector(const avVector<T> &v) const;


  double norm() const;
  void normalize();
  void scale(const double factor);
  
  
  void verbose_mode_set(bool value)
    {verbose_mode = value;}


  // friend T product_scalar(const avVector<T> &v1, const avVector<T> &v2);
  friend avVector<T> operator+ <> (const avVector<T> &v1, const avVector<T> &v2);
  friend avVector<T> operator- <> (const avVector<T> &v1, const avVector<T> &v2);
  friend avVector<T> operator- <> (const avVector<T> &v);
  friend bool operator== <> (const avVector<T> &v1, const avVector<T> &v2);
  friend bool operator!= <> (const avVector<T> &v1, const avVector<T> &v2);
  //friend ostream& operator<< <> (ostream &o, const avVector<T> &v); // nvcc cale
};


#include "./avVector.cc"


#endif /* _VECTOR_ */
