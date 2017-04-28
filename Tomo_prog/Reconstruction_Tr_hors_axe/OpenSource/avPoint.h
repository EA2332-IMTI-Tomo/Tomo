#ifndef _avPOINT_
#define _avPOINT_


/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2001 the GLOS development team (http://glos.sourceforge.net) */


/**
   there seems to be cases where flipping X/Y coordinates is useful,
   since VTK respects normal format but not that stupid SCT.
   
   you can: 
   -flip X/Y of current point (swap_xy())
   -flip input data read to the point (read_flip_xy to true)
   -flip output data written from the point (write_flip_xy to true)
*/


#include <iostream>
#include <string>
#include "macros.h"
#include "Tools.h"



using namespace std;


#define picture(PIC, POINT) \
PIC[POINT.get_X()][POINT.get_Y()]

#define volume(VOL, POINT) \
VOL[POINT.get_Z()][POINT.get_X()][POINT.get_Y()]


double sqr(double x);
float sqr(float x);


/* pre-declaration of classes to suit with following prototypes */
template <typename T>
class avPoint;
template <typename T> // necessary to avoid cycle
class avVector;
template <typename T> // necessary to avoid cycle
class Plane;



#include "avVector.h" // 2d/3d
#include "vArray.h"


/* friend functions prototypes */
/* essai gcc343 */
template <typename T> class avPoint; 
template <typename T> ostream& operator<< (ostream &o, const avPoint<T> &p);
template <typename T> void operator>>  (istream &i, avPoint<T> &p);
/* template <typename T> avPoint<T> operator+  (const avPoint<T> &p, const avVector<T> &v); */
template <typename T> avPoint<T> operator*  (const double &x, const avPoint<T> &p);
 


//----------------------------------------------------------------------

template <typename T>  //class pas propre
class avPoint{

  //provisoire
 public:
  T X, Y, Z;  
  ostream& display(ostream &o) const;
  ostream& display_nice(ostream &o) const;


 protected:
  //T X, Y, Z;
  bool is_valid;
  bool read_flip_xy;
  bool write_flip_xy;


 private:

  bool verbose_mode;
  // sert d'espace de stockage temporaire pour "retourner" une instance sous forme de tableau
  // si on veut une representation duale et surtout capable d'alterner X et Y efficacement, c'est la meilleure option
  T* data;
  avPoint<T>* neighbors; 
  vArray<double> *neighbor_distances;   

  void set_defaults()
    {is_valid = true; data = NULL; neighbors = NULL; read_flip_xy = write_flip_xy = verbose_mode = false;
      neighbor_distances = new vArray<double>(10);}
  

 public:
  
  avPoint() 
    {set_defaults(); set(0, 0, 0);}
  avPoint(T x, T y, T z)
    {set_defaults(); set(x, y, z);}
  avPoint(const avPoint<T> &p) // constructeur par copie, apellé partout où nécessaire
    {set_defaults(); set(p.X, p.Y, p.Z);}
  
  // pour construire une instance a partir d'une autre d'un type template different
  template <typename U> avPoint(const avPoint<U> &p)
    {set_defaults();set((T)(p.get_X()), (T)(p.get_Y()), (T)(p.get_Z()));} 

  avPoint(T* coords)
    {set_defaults(); set(coords[0], coords[1], coords[2]);}

  template <typename U> avPoint(U* coords)
    {set_defaults(); set(T(coords[0]), T(coords[1]), T(coords[2]));}
    
  
  // special: vtk read points revert x/y!
  //avPoint(T* coords, bool) {set(coords[1], coords[0], coords[2]); set_defaults();}
  

  ~avPoint();
      

  READWRITER(T, X); // T get_X();  set_X(T); declared 
  READWRITER(T, Y); 
  READWRITER(T, Z); 
  READER(bool, is_valid); // bool get_is_valid() 
  READWRITER(bool, verbose_mode);
  READWRITER(bool, read_flip_xy); 
  READWRITER(bool, write_flip_xy); 


  // this one is finally called in any case 
  // reverts x/y if asked
  void set(T x, T y, T z);
  void set(const avPoint &p) 
    {set(p.X, p.Y, p.Z);}
  void set(T* coords)
    {set(coords[0], coords[1], coords[2]);}

  // return pointer to a (possibly changing) array, private member of the class
  T* get_data();
  //T* get_data_vtk(); //actually, vtk needs double, with X/Y reverted
 
    
  void write_on(double* coords); 

  /*void write_on_vtk(double* coords) 
    {coords[0] = double(Y); coords[1] = double(X); coords[2] = double(Z);}*/

  /* return neighbor in freeman 8-connex (trigonometric) convention */
  avPoint<T> neighbor_8(int freeman_code) const;
  /* same in pandore (inv-trigo) convention */
  avPoint<T> neighbor_pan(int pandore_code) const;

  /* checks wether points are neigbors */
  bool neighbors_p(const avPoint<T> &Neigh) const;
  
  /* */
  bool in_bounds(T dim_X, T dim_Y, T dim_Z) const
    {return ((X < dim_X) && (Y < dim_Y) && (Z < dim_Z));}

  // euclidean positive distance to another point
  double distance_eucl(const avPoint<T> &p) const;
  
  // euclidean signed distance to a plane identified by a included
  // point and a normal vector (cf softsurfer.com)
  //double distance_signed(const avPoint<T> &p, const avVector<T> &n) const;
  double distance_signed(const Plane<T> &p) const;

  
  // return an avPoint<int> whose coordinates are separately rounded
  avPoint<int> round_coords() const
    {return avPoint<int>((int)round((double)X), (int)round((double)Y), (int)round((double)Z));}
  // return the closest avPoint<int> in 26-connexity
  avPoint<int> round_26() const;
  

  void swap_xy()
    {T tmp; REVERSE_VALUES(X, Y, tmp);}
  void flip()
    {X = -X; Y = -Y; Z = -Z;}


  inline void operator=(const avPoint<T> &p);

  
  bool operator== (const avPoint<T> &p) const;
  bool operator!= (const avPoint<T> &p) const
    {return (! (*this == p));}
  avPoint<T> operator+ (const avPoint<T> &p) const;
  avPoint<T>& operator+= (const avPoint<T> &p);
  avPoint<T> operator- (const avPoint<T> &p) const;
  avPoint<T>& operator-= (const avPoint<T> &p);
  avPoint<T> operator- () const; 
  avPoint<T> operator* (const avPoint<T> &p) const;
  avPoint<T>& operator*= (const avPoint<T> &p);
  avPoint<T> operator* (double t) const;
  avPoint<T>& operator*= (double t);
  avPoint<T> operator/ (double t) const;
  avPoint<T>& operator/= (double t);
  avPoint<T> operator+ (const avVector<T> &v) const;

  T& operator[] (int i);
  

  /* they are friend external functions, not methods: hence cannot be const */  
  // friend ostream& operator<< <> (ostream &o, const avPoint<T> &p);  // nvcc cale
  friend void operator>> <> (istream &i, avPoint<T> &p);
  // friend avPoint<T> operator+ <> (const avPoint<T> &p, const avVector<T> &v); 
  //friend avPoint<T> operator* <> (const double &x, const avPoint<T> &p);
  
};


#include "./avPoint.cc"


#endif /* _POINT_ */






/* #include "avVector.h"  */ // necessary to avoid cycle
/* #include "Plane.h"  */ // necessary to avoid cycle


//prototypes
//template <typename T> avPoint<T> operator-(const avPoint<T> &p1, const avPoint<T> &p2);
//template <typename T> avPoint<T> operator+(const avPoint<T> &p1, const avPoint<T> &p2); 
//template <typename T> bool operator==(const avPoint<T> &p1, const avPoint<T> &p2);
//template <typename T> bool operator!=(const avPoint<T> &p1, const avPoint<T> &p2);
//template <typename T> avPoint<T> convert(const avPoint<U> &p);

