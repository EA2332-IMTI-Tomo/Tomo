/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net) */


template <class T>
Buffer_vector<T>::Buffer_vector(T* elements_array, int nb_elts)
{
  // allocates index normally, allocate and attach only one cluster
  // sized to fit our data. Sets envt values so as to directly match
  // target configuration.
  constructor_called = CONS_CALLED_VALUE;
  this -> cluster_size = nb_elts;
  
  index_capacity = index_init_size = DEFAULT_ID_SIZE; 
  index_extend_size = DEFAULT_IDR_SIZE;
  
  clusters = 1; elements = nb_elts;
  capacity = 0; 
  assert(index_capacity);
  assert(index = (T**)calloc(index_init_size, sizeof(T*)));
  assert(index[0] = (T*)calloc(nb_elts, sizeof(T)));

  for (int i = 0; i < nb_elts; i++) 
    index[0][i] = elements_array[i];
  index_capacity--;
}


template <class T>
void
Buffer_vector<T>::reserve_cluster()
{
  assert(constructor_called_p());

  if (verbose_mode) {cout << "ic rc: " << index_capacity << endl;}
  if (index_capacity <= 0) 
    {
      if (index_elements <= 0) {assert(index = (T**)calloc(index_init_size, sizeof(T*)));}
      else {extend_cluster_table();}
    }
  
  // allocates another cluster and updates elt capacity
  assert(index[index_elements] = (T*)calloc(cluster_size, sizeof(T)));
  capacity += cluster_size;

  index_capacity--;
  index_elements++;
}


template <class T>
void
Buffer_vector<T>::extend_cluster_table()
{
  cerr << endl << "cluster table extended - possible error (bug)" << endl;
  assert(index = (T**)realloc(index, index_elements + index_extend_size));
  index_capacity += index_extend_size;
}


template <class T>
void
Buffer_vector<T>::add_element(const T &elt)
{
  assert(constructor_called_p());

  if (! capacity) // no more free slots: add a cluster
    {
      if (index_capacity == 0) // no more free space in array of clusters 
	if (index_elements == 0) // because the index was not yet allocated
	  {assert(index = (T**)calloc(index_init_size, sizeof(T*)));
	  index_capacity += index_init_size;}
	else extend_cluster_table();
      
      // allocates another cluster and updates elt capacity
      T* new_cluster;
      assert(new_cluster = (T*)calloc(cluster_size, sizeof(T)));
      capacity += cluster_size;

      // adds the element efficiently
      new_cluster[0] = elt;

      // updates index
      index[index_elements] = new_cluster;
      index_capacity--;
      index_elements++;
    }
  else
    {
      set_element(elements, elt);
    }
  
  capacity--;
  elements++;
}


template <class T>
void
Buffer_vector<T>::locate_element(int idx) 
{
  assert(constructor_called_p());

  num_cluster = idx / cluster_size;
  num_elt = idx % cluster_size;
}


template <class T>
T&
Buffer_vector<T>::get_element(int idx) 
{
  locate_element(idx);

  return index[num_cluster][num_elt];
}


template <class T>
void
Buffer_vector<T>::set_element(int idx, const T &elt) 
{
  locate_element(idx);

  index[num_cluster][num_elt] = elt;
}



template <class T>
T&
Buffer_vector<T>::ask_element(int idx) 
{
  assert(idx < elements);
  
  return get_element(idx);
}


template <class T>
ostream& 
operator<<(ostream &o, Buffer_vector<T> &bv)
{
  return (bv.get_verbose_mode()) ? bv.display_nice(o) : bv.display(o);
}


template <class T>
void
operator>>(istream &i, Buffer_vector<T> &bv)
{
  T read;
  i >> read;

  while(! i.eof())
    {
      //cout << "inserting: " << read; 
      bv.add_element(read);
      i >> read;
    }
}


template <class T>
ostream& 
Buffer_vector<T>::display_nice(ostream &o)
{ 
  o << "Buffered Vector: elements nb=" << elements;
  o << " Elements listing:";
  
  for (int i = 0; i < elements; i++)
    o << i << ":" << ask_element(i) << separator;
  
  return o;
}


template <class T>
ostream& 
Buffer_vector<T>::display(ostream &o) 
{
  for (int i = 0; i < elements; i++)
    o << ask_element(i) << separator;

 return o;
}


template <class T>
T& 
Buffer_vector<T>::operator[] (int index) 
{
  return ask_element(index);
}
