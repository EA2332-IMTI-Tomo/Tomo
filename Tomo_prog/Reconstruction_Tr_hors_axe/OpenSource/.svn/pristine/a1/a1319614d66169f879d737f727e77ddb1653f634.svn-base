/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net) */


// =============================================================================


template <class T>
void
vArray<T>::Free()
{
  free(data);
  size = pos_insert = 0;
}


// =============================================================================

template <class T>
void
vArray<T>::fill(T elt) 
{
  fill_in(0, size - 1, elt);
}

// =============================================================================


template <class T>
void
vArray<T>::add(T elt)
{
  if (pos_insert >= size)
    resize(2 * (size + 5));
  
  data[pos_insert] = elt;
  pos_insert++;
  DISPLAY_DATA_MODE = true;
}


// =============================================================================


template <class T>
void
vArray<T>::resize(int newsize)
{
  if (this -> size == 0)
    assert(data = (T*) calloc(newsize, sizeof(T)));
  else
    {
      //swp = data[size - 1];
      assert(data = (T*) realloc(data, newsize * sizeof(T))); 
      //data[size - 1] = swp;
      if (newsize > size) fill_in(size, newsize - 1, 0); 
    }
  this -> size = newsize;
}


// =============================================================================


template <class T>
int 
vArray<T>::position_max() const
{
  int max_pos;
  T max_val; 
  
  assert(size >= 0);

  max_pos = 0;
  max_val = data[max_pos];
  
  for(int idx = 1; idx < size; idx++)
    if(data[idx] > max_val)
      {
	max_pos = idx;
	max_val = data[max_pos];
      }

  return max_pos;
}


// =============================================================================

template <class T>
int 
vArray<T>::position_min() const
{
  int min_pos;
  T min_val; 
  
  assert(size >= 0);

  min_pos = 0;
  min_val = data[min_pos];
  
  for(int idx = 1; idx < size; idx++)
    if(data[idx] < min_val)
      {
	min_pos = idx;
	min_val = data[min_pos];
      }

  return min_pos;
}


// =============================================================================


template <class T>
void
vArray<T>::get_min_max(T &minval, T &maxval)
{
  T min_val, max_val;

  assert(size >= 0);

  min_val = max_val = data[0];
  
  for(size_t idx = 1; idx < size; idx++)
    {
      T val = data[idx];
      if (val < min_val)
	min_val = val;
      else if (val > max_val)
	max_val = val;
    }
     
  minval = min_val;
  maxval = max_val;
}

// =============================================================================


template <class T>
int 
vArray<T>::position(T value) const
{
  assert(size >= 0);
  
  for(int idx = 0; idx < size; idx++)
    if(data[idx] == value)
      return idx;

  return -1;
}

// =============================================================================

template <class T>
int 
vArray<T>::occurences(T value) const
{
  assert(size >= 0);
  int occ = 0;

  for(int idx = 0; idx < size; idx++)
    if(data[idx] == value)
      occ++;

  return occ;
}


// =============================================================================


template <class T>
T
vArray<T>::sum() const
{
  T sum = 0;

  for (int idx = 0; idx < size; idx++)
    sum += data[idx];

  return sum;
}


// =============================================================================

template <class T>
double
vArray<T>::mean() const
{
  return ((double)sum() / (double)size);
}

// =============================================================================

template <class T>
double
vArray<T>::variance() const
{
  double l_mean = mean();
  double sum_sqr = 0;

  for (int idx = 0; idx < size; idx++)
    sum_sqr += SQR(data[idx] - l_mean);
  
  return (sum_sqr / (double)size);
}


// =============================================================================

template <class T>
double
vArray<T>::std_deviation() const
{
  return sqrt(variance());
}

// =============================================================================

template <class T>
vArray<T>&
vArray<T>::operator=(const vArray<T> &a) 
{
  size = a.size; 
  separator = a.separator; //hoping = copies, otherwise sharing same string
  alloc(a.size); 
  copy_in(a.data, a.size);
  return *this;
}

// =============================================================================

template <class T>
ostream& 
operator<<(ostream &o, const vArray<T> &a) 
{
  int limit = a.get_size();
  
  for (int idx = 0; idx < limit; idx++)
    o << a.get_data()[idx] << a.get_separator();

  return o;
}

// =============================================================================

template <class T>
void
operator>>(istream &i, const vArray<T> &a)
{
  int n = 0;
  T read;
  i >> read;

  while((! i.eof()) && (n < a.size))
    {
      a.data[n] = read;
      n++;
      i >> read;
    }
}


// =============================================================================

template <class T>
T&
vArray<T>::operator[] (size_t i) const
{
  return data[i];
}

// =============================================================================

template <class T>
void
vArray<T>::derive()
{
  T prev = data[0], back;

  for (int idx = 1; idx < size; idx++)
    {
      back = data[idx];
      data[idx] -= prev;
      prev = back;
    }

  data[0] = 0;
}

// =============================================================================

template <class T>
void
vArray<T>::absolutize()
{
  for (int idx = 0; idx < size; idx++)
    data[idx] = ABS(data[idx]);
}

// =============================================================================


template <class T>
void 
vArray<T>::map(T (*map_fun)(T))
{
  for (int idx = 0; idx < size; idx++)
    data[idx] = map_fun(data[idx]);
}


// =============================================================================

// default: max = 150

template <class T>
ostream&
vArray<T>::display(ostream &o, size_t max) const
{
  o << endl << "vArray of " << size;
      
  for (size_t i = 0; i < max; i++)
    o << "i:" << i << "=" << data[i] << " | ";
  return o;
}


// =============================================================================


/*
  template <class T>
  void
  vArray<T>::a_add_elt(T *data, int size, T elt)
  {
  for (idx = 0; idx < size; idx++)
  data[idx] += elt;
  }

*/
