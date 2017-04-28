template <typename T>
ostream& 
avVoxel<T>::display_verbose(ostream &o) const
{
  o << endl << "Position: ";
  avPoint<int>::display(o) << "VALUE: " << this -> value;
  return o;
}


template <typename T>
ostream& 
avVoxel<T>::display_io(ostream &o) const
{
  avPoint<int>::display(o) << " " << this -> value;
  return o;
}


template <typename T>
ostream&
operator<<(ostream &o, const avVoxel<T> &p) 
{
  if (p.verbose_mode)
    p.display_verbose(o);
  else
    p.display_io(o);

  return o; 
}


template <typename T>
void
avVoxel<T>::operator=(const avVoxel<T> &p) {
  this -> value = p.value;
  avPoint<int>::set(p.X, p.Y, p.Z);
}
