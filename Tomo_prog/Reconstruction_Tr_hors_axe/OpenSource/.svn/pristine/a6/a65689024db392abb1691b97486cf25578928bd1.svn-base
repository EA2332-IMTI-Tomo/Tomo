#ifndef _VECTRA_MATRIX_
#define _VECTRA_MATRIX_


/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2001 the GLOS development team (http://glos.vectraproject.com) */


#include <iostream>
#include <string>
#include <math.h>

#include "macros.h"


using namespace std;



/* pre-declaration of classes to suit with following prototypes */
template <typename T> class vMatrix;

/* pre-declaration of friend functions, demanded by gcc3.4.* */
template <typename T> ostream& operator<< (ostream &o, const vMatrix<T> &m);


// =============================================================================


template <class T>
class vMatrix
{

 private:
  
  size_t nb_rows, nb_cols;
  T** data;


 public:
  
  vMatrix(size_t nb_rows, size_t nb_cols, T initval = 0);
  // creates a matrix object from an existing T** array. Data is not copied
  // WARNING: creates critical errors at destruction/deallocation if sizes not consistent 
  vMatrix(size_t nb_rows, size_t nb_cols, T** data);
  vMatrix(const vMatrix& m);
  ~vMatrix();
  
  
  // declares get_nb_rows() and get_nb_cols()
  READER(size_t, nb_rows);
  READER(size_t, nb_cols);

  // fills matrix with given value
  void fill(T value);
  // multiplies matrix elements by scalar
  void scalar_mult(T value);
  // map a 1-arg function to matrix elements
  void map(T (*map_fun)(T));
  
  // instance courante reçoit le produit de m1 {a, b} par m2 {b, c} 
  // ssi instance courante de dimension {a, c}
  void multiply(const vMatrix<T> &m1, const vMatrix<T> &m2);
  // compute result of mv * m: {a, a} * {a, b}
  // mv is square, full of zeros, and his diagonal contains values from vector {1, a}
  // this: {a, b}
  void multiply_from_diag(const vMatrix<T> &vector, const vMatrix<T> &m);
  void transpose(const vMatrix<T> &m);

  // override m[] to return row as array pointer 
  // and m[][] to return element as element pointer
  T* operator[] (int nlig) {return data[nlig];} 
  T* operator[] (int nlig) const {return data[nlig];} 
  
    
  void display(ostream &o) const;
  friend ostream& operator<< <> (ostream& o, const vMatrix<T> &m);  
  // piège: friend impose usage du <> vide, uniquement dans le .h
  bool display_limited;
  size_t display_limit;


  // opérateur de recopie à définir: par défaut, les deux matrices partageraient la même plage mémoire
  // n = m (ssi n et m de même dimension)
  void operator= (const vMatrix<T> &m);
  
};




// =============================================================================

#include "./vMatrix.cc"


#endif // _VECTRA_MATRIX_
