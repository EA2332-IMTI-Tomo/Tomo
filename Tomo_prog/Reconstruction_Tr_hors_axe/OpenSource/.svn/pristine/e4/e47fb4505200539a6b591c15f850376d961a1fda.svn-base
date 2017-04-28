/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net) */



#ifndef _BUFFER_VECTOR_
#define _BUFFER_VECTOR_


/*** 
     Je pense qu'il faut flinguer cette classe et utiliser soit Array, soit une classe BArray a ecrire, qui se baserait sur des Arrays.
     A force de tout vouloir gerer ici, on se retrouve avec des erreurs incoercibles.
     Ceci dit, la classe est tres utile
*/


#include <iostream>
#include <string>

#include "macros.h"


using namespace std;



#define DEFAULT_CL_SIZE 20
#define DEFAULT_ID_SIZE 50
#define DEFAULT_IDR_SIZE 5
#define CONS_CALLED_VALUE 123




/* pre-declaration of classes to suit with following prototypes */
template <typename T>
class Buffer_vector;
template <typename T> ostream& operator<<  (ostream &o, Buffer_vector<T> &bv); 
template <typename T> void operator>> (istream &i, Buffer_vector<T> &bv);

template <class T>
class Buffer_vector
{
 private:
  int cluster_size;
  int clusters;
  int elements;
  int capacity;

  int index_capacity;
  int index_elements;
  int index_init_size;
  int index_extend_size;
  
  T** index;


  BOOL verbose_mode;
  int constructor_called;

  string separator;

  // precomputations
  int num_cluster, num_elt;

  // use of ask_element updates internal variables
  ostream& display(ostream &o); 
  ostream& display_nice(ostream &o);


  // constructor might not be called if user has the bad idea to
  // allocate in C mode an array of objects: methods will be callable
  // on these without constructor call. A disaster here!
  bool constructor_called_p()
    {return (constructor_called == CONS_CALLED_VALUE);}

  void set()
    {
      constructor_called = CONS_CALLED_VALUE;
      verbose_mode = FALSE;
      cluster_size = DEFAULT_CL_SIZE; 
      index_init_size = DEFAULT_ID_SIZE; index_extend_size = DEFAULT_IDR_SIZE;
      clusters = elements = capacity = index_elements = index_capacity = 0;  
      separator = " ";
    }
  

  // finds indexes of given element, assuming correctness of idx
  inline void locate_element(int idx);   
  // returns given element, assuming idx correct
  inline T& get_element(int idx); 
  inline void set_element(int idx, const T &elt); // sets...
  // if no more index room for new clusters, realloc index
  void extend_cluster_table();
  // grants direct access to asked cluster: proceed at your very own
  // risk - myvec::get_cluster(x) - if you use it. Possible errors:
  // cluster not existing, elements that will be erased by the
  // progress of the ones behind.
  T* get_cluster(int idx)
    {return index[idx];}
  
  
 protected:

 public:
  Buffer_vector()
    {set();}
  Buffer_vector(int cluster_size)
    {set(); this -> cluster_size = cluster_size;}
  Buffer_vector(int cluster_size, int index_init_size)
    {set(); this -> cluster_size = cluster_size; this -> index_init_size = index_init_size;}
  Buffer_vector(int cluster_size, int index_init_size, int index_extend_size)
    {set(); this -> cluster_size = cluster_size; this -> index_init_size = index_init_size; this -> index_extend_size = index_extend_size;}
  // efficiently builds an object instance from an array 
  Buffer_vector(T* elements_array, int nb_elts);      
  // copy constructor:
  Buffer_vector(const Buffer_vector &b)
    {cerr << "copy cons. not implemented!" << endl; exit(1);}


  // operateur pour affecter une instance existante aux valeurs d'une autre
  // a = b <=> a.operator=(b)
  Buffer_vector& operator=(const Buffer_vector &v)
    {cerr << "Buffer_vector: = not implemented!" << endl; exit(1);}


  
  ~Buffer_vector()
    {for (int i = 0; i < index_elements; i++) free(index[i]);
    free(index);}


  // adds a new cluster to the vector, whose elements are
  // zeros. interest is just to control when allocation might occur.
  void reserve_cluster();
  // if buffer is new, allocates a cluster in a single pass filled
  // with zero-elements. Their array is returned to you so as to let
  // you play directly with it. Interest is to emulate the efficiency
  // of a standard array.
  T* grant_first_cluster()
    {assert(elements == 0); reserve_cluster(); 
    elements += cluster_size; capacity = 0; return index[0];}
  
        

  READER(int, cluster_size);
  READER(int, elements);
  READER(int, capacity);
  READWRITER(BOOL, verbose_mode);
  READWRITER(string, separator);


  void set_cluster_size(int cluster_size)
    {assert(elements == 0); this -> cluster_size = cluster_size;}  

  BOOL cluster_resizeable_p()
    {return (elements == 0);} const

  void resize_clusters(int new_cluster_size)
    {assert(cluster_resizeable_p()); cluster_size = new_cluster_size;}

  void add_element(const T &elt);
  // indix-secure R/W access
  inline T& ask_element(int idx); //needs to update internal variables


  inline T& operator[] (int index);

  // write current instance into stream
  friend ostream& operator<< <> (ostream &o, Buffer_vector<T> &bv);  // use of view_elt -> update
  // read extra elements from an instance from stream (keyboard, existing file)
  friend void operator>> <> (istream &i, Buffer_vector<T> &bv);
  /* hint: if reading from keyboard, insertion(s) will occur after
     each newline, the procedure being terminated by ^D */
};


#include "./Buffer_vector.cc"


#endif /* _BUFFER_VECTOR_ */


// nota bene: inline works here since we have a template and that the
// code of the function is included in any calling file. inline is a
// kind of macro/template...
