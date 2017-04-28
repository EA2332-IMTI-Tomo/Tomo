#ifndef _VECTRA_TUPLE_
#define _VECTRA_TUPLE_

using namespace std;


template <class T> class vTuple;
template <class T> ostream& operator<< (ostream &o, const vTuple<T>& vt);


template <class T>
class vTuple{ 

 public:
  T car;
  T cdr;

 vTuple(T a, T b):car(a), cdr(b) {}
  
 void set(T a, T b)
   {car = a; cdr = b;}
 
 ostream& display(ostream &o) const
   { 
     o << "{ " << car << "; " << cdr << " }";
     return o;
   }
};


template <class T>
ostream&
operator<< (ostream &o, const vTuple<T>& vt) 
{
  vt.display(o);
  return o;
}


#endif // _VECTRA_TUPLE_
