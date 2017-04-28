template <class T>
BOOL
avSegment<T>::intersects_plane_p(const Plane<T> &pl)
{
  avVector<T> up(pl.get_inavPoint() - this -> A);
  avVector<T> down(this -> B - this -> A);
  
  double upterm = pl.get_normal().product_scalar(up);
  double downterm = pl.get_normal().product_scalar(down);
  
  verbose_mode = TRUE;
  t_inter_set = TRUE;
  if (downterm == 0.0) 
    {
      this -> parallel = TRUE;
      return FALSE;
    }

  this -> parallel = FALSE;
  this -> t_inter = upterm / downterm;
  
  return (BETWEEN(this -> t_inter, 0.0, 1.0)); 
}


template <class T>
avPoint<T>
avSegment<T>::intersection_plane(const Plane<T> &pl) 
{
  if (!(this -> t_inter_set))
    {
      BOOL inter = intersects_plane_p(pl);
      assert(inter);
    }
  else
    {assert (!(this -> parallel));}
  
  this -> t_inter_set = FALSE;
  return this -> parametric(this -> t_inter);
}


template <class T>
avPoint<T>
avSegment<T>::middle() const
{
  avPoint<T> mid(A + B);
  mid /= 2;
  return mid;
}


template <class T>
void
avSegment<T>::set_middle(avPoint<T> &dest) const
{
  dest.set(T(0), T(0), T(0));
  dest += this -> A;
  dest += this -> B;
  dest /= 2;
}



template <class T>
avPoint<T>
avSegment<T>::parametric(double t) const
{
  avVector<T> scaled_vec(this -> vector);
  scaled_vec.scale(t);
  return (this -> A + scaled_vec);
}


template <class T>
ostream& 
operator<<(ostream &o, const avSegment<T> &s) 
{
  o << "first pt: " << s.A << "second pt: " << s.B; 

  return o;
}
