#ifndef _avVOXEL_
#define _avVOXEL_


#include <iostream>
#include <string>
#include "macros.h"
#include "avPoint.h"


using namespace std;


/* a voxel point should consist of an object avPoint<int> for keeping
   information about its location completed (by inheritance) by a
   value of type T representing the color value of avPoint in the volume 
   
   more sophisticated version than SEG's Vol_avPoint
*/

//template <typename I>
//class Voxel;
//template <typename I> ostream& operator<< (ostream &o, const Voxel<I> &p);


template <typename T>
class avVoxel: public avPoint<int>
{
 protected:
  BOOL verbose_mode;
  
 private:
  T value;
  ostream& display_io(ostream &o) const;
  ostream& display_verbose(ostream &o) const;

  void init() {verbose_mode = FALSE;}

 public:
  Voxel(T init_val)
    {init(); this -> value = init_val;}
  Voxel(int x, int y, int z, T init_val)
    {init(); this -> set(x, y, z); this -> value = init_val;}
  Voxel(const avPoint<int> &p, T init_val)
    {init(); this -> set(p.X, p.Y, p.Z); this -> value = init_val;}

  ~avVoxel() {;}


  READWRITER(T, value); // I get_value(); I set_value(); declared

  avPoint<int> toPoint() const
    {return avPoint<int>(X, Y, Z);}

  inline void operator=(const avVoxel<T> &p);

  /* they are friend external functions, not methods: hence cannot be
     const */
  //friend ostream& operator<< <> (ostream &o, const Voxel<I> &p); nvcc cale

  void verbose_mode_set(BOOL value)
    {verbose_mode = value;}
};


#include "avVoxel.cc"


#endif /* _avVOXEL_ */

