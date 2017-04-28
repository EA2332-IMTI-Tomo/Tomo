// inclus par .h


template <class T>
T
avVector<T>::product_scalar(const avVector<T> &v) const
{
  return (this->X * v.X) + (this->Y * v.Y) + (this->Z * v.Z); 
}


//v_U ^ v_V : (yz' - zy') + (zx' - xz') + (xy' - yx')
template <class T>
avVector<T>
avVector<T>::product_vector(const avVector<T> &v) const
{
  avVector<T> R;
  avVector<T> C(*this);
  R.X = C.Y * v.Z - C.Z * v.Y;
  R.Y = C.Z * v.X - C.X * v.Z;
  R.Z = C.X * v.Y - C.Y * v.X;
  return R;
}
// on peut pas passer en ref: il FAUT retourner une copie de la
// variable allouee localement pour le resultat


/*
template <class T>
T
product_scalar(const avVector<T> &v1, const avVector<T> &v2)
{
  return (v1.X * v2.X) + (v1.Y * v2.Y) + (v1.Z * v2.Z);
}
*/


template <class T>
double
avVector<T>::norm() const
{
  return this->toavPoint().distance_eucl(avPoint<T>(0, 0, 0)); 
}


template <class T>
void
avVector<T>::scale(const double s)
{
  this->X = (T)(s * this->X);
  this->Y = (T)(s * this->Y);
  this->Z = (T)(s * this->Z);
}


template <class T>
void
avVector<T>::normalize()
{
  this -> scale(1 / this -> norm());
}


template <class T>
bool
operator==(const avVector<T> &v1, const avVector<T> &v2)
{
  return ((v1.X == v2.X) && (v1.Y == v2.Y) && (v1.Z == v2.Z));
}


template <class T>
bool
operator!=(const avVector<T> &v1, const avVector<T> &v2)
{
  return ! (v1 == v2);
}


template <class T>
avVector<T> 
operator+ (const avVector<T> &v1, const avVector<T> &v2)
{
  return avVector<T>(v1.X + v2.X, v1.Y + v2.Y, v1.Z + v2.Z);
}


template <class T>
avVector<T> 
operator- (const avVector<T> &v1, const avVector<T> &v2)
{
  return avVector<T>(v1.X - v2.X, v1.Y - v2.Y, v1.Z - v2.Z);
}


template <class T>
avVector<T> 
operator- (const avVector<T> &v)
{
  return avVector<T>(-v.X, -v.Y, -v.Z);
}


template <class T>
ostream& 
operator<<(ostream &o, const avVector<T> &p) 
{
  o << "avVector: ";  
  return p.display(o);
}


