#ifndef _VECTRA_ARRAY_
#define _VECTRA_ARRAY_

/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net) */


#include <iostream>
#include <string>
#include <math.h>

#include "macros.h"
#include <assert.h> // a virer

using namespace std;


#define DEFAULT_SEPARATOR " "





/* pre-declaration of classes to suit with following prototypes */
template <typename T> class vArray;

/* pre-declaration of friend functions, demanded by gcc3.4.* */
template <typename T> ostream& operator<< (ostream &o, const vArray<T> &a);
template <typename T> void operator>> (istream &i, const vArray<T> &a);

//----------------------------------------------------------------------


// This class is here to perform basic operations upon a simple 1D array
// It is assumed that the user checked the size of the array before accessing it
// It can work on top of an allocated array. It can provide its array for external use.
// like an array, you need to unallocate by hand


template <class T>
class vArray
{

  // =============================================================================

 protected:

  T swp;

  int size;
  int pos_insert; // for add: position to insert elts incrementally
  bool DISPLAY_DATA_MODE;
  T* data;
  char* separator; // separator between elements when display
  
  
  void set_defaults() { }

  void alloc(int size) {
    ARRAY_ALLOC(data, size, T);
  } // incompatible avec realloc: data = new T(size); genere des erreurs indebogables avec 2 instances

  void fill_in(int start, int end, T value) {
    for (int idx = start; idx <= end; idx++) 
      data[idx] = value;
  }
  
  void copy_in(T* array, int length) {
    for (int idx = 0; idx <= length; idx++) 
      data[idx] = array[idx];
  }
   
   
  // =============================================================================
  // constructeurs, destructeurs
  // =============================================================================

 public:
    
  vArray() {
    size = 0; separator = DEFAULT_SEPARATOR; pos_insert = 0; DISPLAY_DATA_MODE = false;
  }

  vArray(int size) {
    alloc(size); this -> size = size; separator = DEFAULT_SEPARATOR; pos_insert = 0; DISPLAY_DATA_MODE = false;
  }

  // WARNING: the array is NOT copied, which is NOT a solution to use in general
  // the idea is to propose nice operations onto a standard array already allocated in a program
  vArray(T* existing_array, int size)
    {this -> size = size; data = existing_array; separator = DEFAULT_SEPARATOR; pos_insert = size; DISPLAY_DATA_MODE = false;}

  vArray(const vArray& a)
    {size = a.size; separator = a.separator; pos_insert = 0;//hoping = copies, otherwise sharing same string
      DISPLAY_DATA_MODE = false;
      alloc(a.size); copy_in(a.data, a.size);}
  
  // utiliser Free pour désallouer à la main
  ~vArray() {
  }

  vArray& operator=(const vArray &a);


  // =============================================================================
  // méthodes
  // =============================================================================

 public:
  

  READER(int, size);
  READER(T*, data);
  READWRITER(char*, separator);
  
  void Free();
  
  
  // set an array element without control
  // see also read operator[int] 
  void set(int pos, T elt) {
    data[pos] = elt;
  }

  void resize(int newsize);
  
  void reverse() {
    T swap;
    size_t stop = size / 2;

    for (size_t i = 0; i < stop; i++)
      {
	REVERSE_VALUES( data[i], data[ size - 1 - i], swap );
      }
  }

  //----------------------------------------------------
  // all these operations are performed on ALL the array
  // if you built your array with add(), please call

  // return position of max/minimum element
  int position_max() const;
  int position_min() const;

  // return max/minimum element
  T max() const {
    return data[position_max()];
  }
  T min() const {
    return data[position_min()];
  }
  void get_min_max(T &minval, T &maxval);


  // return position of given element, -1 if not found
  int position(T value) const;
  bool contains(T value) const {
    return (position(value) > -1);
  }
  int occurences(T value) const;
  T sum() const;
  double mean() const;
  double variance() const;
  double std_deviation() const;
  // computes finite differences onto given array from first element (set to 0) (side effect)
  void derive();
  void absolutize(); //(side effect)
  

  // maps map_fun onto current array
  void map(T (*map_fun)(T));
  
    
  void fill(T value);
  void fill_random();
  


  // friend ostream& operator<< <> (ostream &o, const vArray<T> &a); // nvcc cale
  //file read will modify vArray elements from 0 to last read element OR array size 
  friend void operator>> <> (istream &i, const vArray<T> &a);
  // for implementation, inspire from parsing in binarize
  T& operator[] (size_t i) const;


  ostream& display(ostream &o, size_t max = 150) const;


  



  //----------------------------------------------------
  // Ces méthodes sont prévues pour émuler le fonctionnement d'une File

  // in order to fill array from position 0
  // adds elt to last position.
  // reallocation performed if necessary
  void add(T elt);
  // works even if initial size is zero

  // return number of elements added. can be considered as number of
  // "array items" if you only entered elements though add.
  int get_added_count() const
  {return pos_insert;}

  // re-inits the "queue" in the sense of adding.
  // actually, set all entries to 0 and next add to position 0
  void reset()
  {fill_in(0, pos_insert - 1, (T)0); pos_insert = 0;}

  //in order to switch back to a 'classical' array, we have to
  //truncate the zeroes after last element and until actual end of
  //array (size).
  //this is performed in setting 'size' to a smaller value: little
  //memory is actually lost, but CPU overhead is minimal (realloc can
  //copy even for shrinking). 
  void to_classical()
  {size = pos_insert;}


};

  


#include "./vArray.cc"


#endif /* _ARRAY_ */
