#ifndef _VECTOR_
#define _VECTOR_


/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2001 the GLOS development team (http://glos.sourceforge.net) */


#include <iostream>
#include <string>

#include "macros.h"


#include "Point.h" // really used



using namespace std;


/* pre-declaration of classes to suit with following prototypes */
template <class T>
class Vector;


template <class T>
class Vector
{
 private:
 protected:
 public:
  ~Vector() {;}  
};


#include "./Vector.cc"


#endif /* _VECTOR_ */
