#ifndef _BUFFERED_VECTOR_
#define _BUFFERED_VECTOR_


#include "macros.h"
#include <vector>
#include <iostream>


using namespace std;


template <class T>
class buf_vector: public vector<T>
{

 private:
  int cluster_size;
  ostream& display(ostream &o) const;

 public:
  buf_vector(int cluster_size)
    {
      buf_vector<T>::cluster_size = cluster_size;
      reserve(cluster_size);
      // cout << endl << cluster_size << "allocated" 
      //   << capacity() << " ";
    }

  /*
    buf_vector(const buf_vector<T> & other):vector<T>(other)
    {
      cluster_size = other.cluster_size;
      cout << " recopie() " << size() << " " << capacity() << endl;
      cout << " recopie() " << other.size() << " " << other.capacity() << endl;
    }
  */

  buf_vector()
    {cerr << "not instantiable without cluster size"; exit(1);} 
  

  inline void my_push_back(const T &val);
  friend ostream& operator<< <T> (ostream &o, const buf_vector<T> &v);
  

  READER(int, cluster_size);

};


#include "buf_vector.cc"


#endif _BUFFERED_VECTOR_
