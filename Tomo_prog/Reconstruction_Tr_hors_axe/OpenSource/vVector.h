#ifndef _VECTRA_MINIVECTOR_
#define _VECTRA_MINIVECTOR_

/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net) */


// minimal dyn-array class for handling reallocations and maintaining 1-D array compatibility
// at any point, you have a simple 1-D array and a size
// it just lets you reallocate size and/ or access and insert safely


#include <iostream>
#include <string>
#include <math.h>

#include "macros.h"


using namespace std;





template <class T>
class vVector
{
 private:

  T swp;

  size_t size;          // taille de data
  size_t extend_size;   // de combien on réallouera
  size_t pos_insert;    // si on utilise insert, index de la première case libre

  T* data;
  
  void self_alloc(size_t size);
  // réallocation de taille < ou >, peu safe
  void self_realloc(size_t new_size); 
  
 public:

  
  vVector(size_t init_size)
    {size = init_size; ; pos_insert = 0; extend_size = init_size; self_alloc(init_size);}
  vVector(size_t init_size, size_t extend_size)
    {size = init_size; ; pos_insert = 0; this -> extend_size = extend_size; self_alloc(init_size);}
  ~vVector() 
    { free(data); }


  READER(size_t, size);
  READER(size_t, pos_insert);
  READWRITER(size_t, extend_size);
  READER(T*, data);
  
  // si on utilise le mode sécure, extend est invoqué tout seul
  void extend(); // réallocation de oldsize + resize_clust


  // --------------------------------------------------
  // accès aléatoire

  // raw, sans sécurité. secure, avec test préalable

  inline void r_set(const size_t pos, const T elt)
  { data[pos] = elt; }

  inline T r_get(const size_t pos) const 
  { return data[pos]; }
  

  inline void s_set(const size_t pos, const T elt)
  { if (pos > size - 1) extend(); r_set(pos, elt); } 

  inline T s_get(const size_t pos) const
  { assert(pos < size); return data[pos]; } 

  
  // --------------------------------------------------
  // accès séquentiel

  inline void r_insert(const T elt)
  { data[pos_insert] = elt; pos_insert++; }
  
  inline void s_insert(const T elt)
  { if (pos_insert == size) extend(); 
    data[pos_insert] = elt; pos_insert++; }


  // --------------------------------------------------
  // affichage

  // par défaut, on affiche le tableau du début jusqu'au dernier élément inséré
  // si on a effectué des insertions aléatoires, c'est à nous de spécifier jusqu'où on va

  void display(ostream &o, size_t until = 0, const char* separator = " ")
  {
    if (! until)
      until = pos_insert;
    assert (until < size);

    for (size_t i = 0; i < until; i++)
      o << data[i] << separator;

    //return o;
  }
  
};

  


#include "./vVector.cc"


#endif 



// cf section Tests pour usage basique
