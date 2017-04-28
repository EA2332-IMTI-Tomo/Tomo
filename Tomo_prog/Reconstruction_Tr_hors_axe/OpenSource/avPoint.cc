inline double
sqr(double x)
{
  return x * x; 
}

inline float
sqr(float x)
{
  return x * x; 
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>::~avPoint() 
{
  if (data != NULL) delete [] data;
  if (neighbors != NULL) delete [] neighbors;
  // deletion of neighbor_distances will be automatic (hope so)
}


//------------------------------------------------------------------------------

template <typename T>
void
avPoint<T>::set(T x, T y, T z)
{
  if (read_flip_xy) 
    {X = y; Y = x;}
  else 
    {X = x; Y = y;}
  Z = z;
}


//------------------------------------------------------------------------------

template <typename T>
T*
avPoint<T>::get_data()
{
  if (data == NULL) data = new T[3];
  if (write_flip_xy)
    {data[0] = Y; data[1] = X;} 
  else 
    {data[0] = X; data[1] = Y;} 
  data[2] = Z;
  return data;
}


//------------------------------------------------------------------------------

template <typename T>
void
avPoint<T>::write_on(double* coords) 
{
  if (write_flip_xy)
    {coords[0] = double(Y); coords[1] = double(X);} 
  else 
    {coords[0] = double(X); coords[1] = double(Y);} 
  coords[2] = double(Z);
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template <typename T>
avPoint<T>
avPoint<T>::neighbor_8(int freeman_code) const
{
  T d_X = 0, d_Y = 0;
  int code = ABS(freeman_code) % 8;
  
  switch (code)
    {
    case 0: d_X = 1; break;
    case 1: d_X = 1; d_Y = 1; break;
    case 2: d_Y = 1; break;
    case 3: d_X = -1; d_Y = 1; break;
    case 4: d_X = -1; break;
    case 5: d_X = -1; d_Y = -1; break;
    case 6: d_Y = -1; break;      
    case 7: d_X = 1; d_Y = -1; break;
    default:
      cerr << "invalid freeman code"; exit(1);
    }

  return avPoint<T>(X + d_X, Y + d_Y, Z);
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>
avPoint<T>::neighbor_pan(int pandore_code) const
{
  T d_X = 0, d_Y = 0;
  int code = (pandore_code % 8);
  
  switch (code)
    {
    case 2: d_X = 1; break;
    case 1: d_X = d_Y = 1; break;
    case 0: d_Y = 1; break;
    case 7: d_X = -1; d_Y = 1; break;
    case 6: d_X = -1; break;
    case 5: d_X = d_Y = -1; break;
    case 4: d_Y = -1; break;      
    case 3: d_X = 1; d_Y = -1; break;
    default:
      cerr << "invalid pandore code"; exit(1);
    }

  return avPoint<T>(X + d_X, Y + d_Y, Z);
}


//------------------------------------------------------------------------------

template <typename T>
bool
avPoint<T>::neighbors_p(const avPoint<T> &voisin) const
{
  return ((ABS(X - voisin.X) <= 1) 
	  && (ABS(Y - voisin.Y) <= 1) 
	  && (ABS(Z - voisin.Z) <= 1));
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template <typename T>
double
avPoint<T>::distance_eucl(const avPoint<T> &p) const
{
  return sqrt( sqr((float)(X - p.get_X())) 
	       + sqr((float)(Y - p.get_Y()))
	       + sqr((float)(Z - p.get_Z())));  
}


//------------------------------------------------------------------------------

template <typename T>
double
avPoint<T>::distance_signed(const Plane<T> &p) const
{
  avVector<T> N = p.get_normal();
  avPoint<T> B = p.get_inavPoint();

  avVector<T> ray(*this - B);

  return ( (double)(N.product_scalar(ray)) / N.norm());
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// local-use macro
#define CREATE_NEIGHBOR_PT(A, B, C) \
neighbors[i].set(A, B, C); \
i++;


template <typename T>
avPoint<int>
avPoint<T>::round_26() const
{
  int minpos = -1;
  int x, y, z;
  int i = 0;

  if (neighbors.get_size == 0) 
    {
      neighbors = new avPoint<T>(26); 
      neighbor_distances -> resize((int)26);
    }
  
  avPoint<T> rounded_central(*this);
  rounded_central.round_coords();
  x = rounded_central.get_X();
  y = rounded_central.get_Y();
  z = rounded_central.get_Z();

  // sets default-allocated avPoints<T> of Array with computed
  // coordinates of avPoint<int> neighbors
  AROUND_26_CONNEX(x, y, z, CREATE_NEIGHBOR_PT);
  
  for (i = 0; i < 26; i++)
    neighbor_distances[i] = rounded_central.distance_eucl(neighbors[i]);
  
  minpos = neighbor_distances -> position_min();
  return avPoint<int>(neighbors[minpos]);
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template <typename T>
ostream& 
operator<<(ostream &o, const avPoint<T> &p) 
{
  return (p.get_verbose_mode()) ? p.display_nice(o) : p.display(o);
}


//------------------------------------------------------------------------------

template <typename T>
void
operator>>(istream &i, avPoint<T> &p)
{
  T read_X, read_Y, read_Z;
  bool stream_over;

  i >> read_X;
  i >> read_Y;
  
  stream_over = i.eof();

  i >> read_Z;

  // returned point is incomplete:
  // stream probably contained end at newline  
  if (stream_over) 
    {
      read_X = read_Y = read_Z = 0;
      p.is_valid = FALSE;
    }
    
  p.X = read_X;
  p.Y = read_Y;
  p.Z = read_Z;
}


//------------------------------------------------------------------------------

template <typename T>
ostream&
avPoint<T>::display_nice(ostream &o) const
{
  o << endl << " ("
    << avPoint<T>::X << " ; " << avPoint<T>::Y << " ; " << avPoint<T>::Z 
    << ") ";
  return o;
}


//------------------------------------------------------------------------------

template <typename T>
ostream&
avPoint<T>::display(ostream &o) const
{
  //   o << avPoint<T>::X << " " << avPoint<T>::Y << " " << avPoint<T>::Z << " ";
  o << this -> X << " " << this -> Y << " " << this -> Z << " ";
  return o;
}


//------------------------------------------------------------------------------

template <typename T>
void
avPoint<T>::operator=(const avPoint<T> &p) {
  X = p.X;
  Y = p.Y;
  Z = p.Z;
}


//------------------------------------------------------------------------------

template <typename T>
bool
avPoint<T>::operator==(const avPoint<T> &p) const 
{
  return ( (X == p.X) && (Y == p.Y) && (Z == p.Z) );
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>
avPoint<T>::operator+(const avPoint<T> &p) const
{
  return avPoint<T>(X + p.X, Y + p.Y, Z + p.Z);
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>
avPoint<T>::operator-(const avPoint<T> &p) const
{
  return avPoint<T>(X - p.X, Y - p.Y, Z - p.Z);
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>
avPoint<T>::operator*(const avPoint<T> &p) const
{
  return avPoint<T>(X * p.X, Y * p.Y, Z * p.Z);
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>
avPoint<T>::operator*(double t) const
{
  return avPoint<T>((T)(X * t), (T)(Y * t), (T)(Z * t));
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>&
avPoint<T>::operator+=(const avPoint<T> &p)
{
  X += p.X, Y += p.Y, Z += p.Z;
  return *this;
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>&
avPoint<T>::operator-=(const avPoint<T> &p)
{
  X -= p.X, Y -= p.Y, Z -= p.Z;
  return *this;
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>&
avPoint<T>::operator*=(const avPoint<T> &p)
{
  X *= p.X, Y *= p.Y, Z *= p.Z;
  return *this;
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>&
avPoint<T>::operator*=(double t)
{
  X = (T)(X * t); Y = (T)(Y * t); Z = (T)(Z * t);
  return *this;
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>&
avPoint<T>::operator/=(double t)
{
  X = (T)(X / t); Y = (T)(Y / t); Z = (T)(Z / t);
  return *this;
}

//------------------------------------------------------------------------------

/*
template <typename T>
avPoint<T>
operator+(const avPoint<T> &p, const avVector<T> &v) 
{
  return avPoint<T>(p.X + v.X, p.Y + v.Y, p.Z + v.Z);
}
*/

template <typename T>
avPoint<T>
avPoint<T>::operator+(const avVector<T> &v) const
{
  return avPoint<T>(X + v.X, Y + v.Y, Z + v.Z);
}

//------------------------------------------------------------------------------

template <typename T>
avPoint<T>
operator*(const double &k, const avPoint<T> &p) 
{
  return avPoint<T>(k * p.X, k * p.Y, k * p.Z);
}


//------------------------------------------------------------------------------

template <typename T>
avPoint<T>
operator/(const double &k, const avPoint<T> &p) 
{
  return avPoint<T>(k / p.X, k / p.Y, k / p.Z);
}

//------------------------------------------------------------------------------

template <typename T>
avPoint<T>
avPoint<T>::operator-() const
{
  return avPoint<T>(- X, - Y, - Z);
}


//------------------------------------------------------------------------------

template <typename T>
T&
avPoint<T>::operator[] (int i)
{
  switch (i)
    {
    case 0: 
      return X;
    case 1:
      return Y;
    case 2:
      return Z;
    default:
      cerr << "out of bounds avPoint index" << endl;
      exit(1);
    }
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


/*
template <typename T>
avPoint<T> 
convert(const avPoint<U> &p)
{
  return avPoint<T>((T)p.get_X(), (T)p.get_Y(), (T)p.get_Z());
}
*/


/*
template <typename T>
T*
avPoint<T>::get_data_vtk()
{
  if (data == NULL) data = new T[3];
  data[0] = Y; 
  data[1] = X; 
  data[2] = Z;
  return data;
}
*/


/** EXEMPLE DE SPECIALISATION!!! */
/*
double
avPoint<int>::distance_eucl(const avPoint<int> &p) const
{
return sqrt( sqr((float)(X - p.get_X())) 
+ sqr((float)(Y - p.get_Y()))
+ sqr((float)(Z - p.get_Z())));  
}
*/
