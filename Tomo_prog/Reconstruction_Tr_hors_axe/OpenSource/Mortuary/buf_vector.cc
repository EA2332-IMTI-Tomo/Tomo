template <class T>
void 
buf_vector<T>::my_push_back(const T &val)
{
  assert(capacity());

  if (size() == capacity())
    {
      reserve(size() + cluster_size);
      cerr << endl << "needed to realloc:" << size();
    }
  //   cout << "OK\n";
  vector<T>::push_back(val);
}


template <class T>
ostream&
buf_vector<T>::display(ostream &o) const
{
  vector<T>::iterator i;
  //o << endl << "size: " << size() << "content: ";
  
  i = begin();
  while (i != end())
    {
      o << *i << " ";
      i++;
    }
  
  return o;
}


template <class T>
ostream& 
operator<< (ostream &o, const buf_vector<T> &v)
{
  v.display(o);
  return o;
}
