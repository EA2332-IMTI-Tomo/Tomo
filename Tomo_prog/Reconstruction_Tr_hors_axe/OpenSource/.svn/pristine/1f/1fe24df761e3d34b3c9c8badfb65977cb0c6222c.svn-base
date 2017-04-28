template <class T>
void 
vVector<T>::self_alloc(size_t size)
{
  ARRAY_ALLOC(data, size, T);
}


template <class T>
void 
vVector<T>::self_realloc(size_t new_size)
{
  ARRAY_REALLOC(data, new_size, T);
}


template <class T>
void 
vVector<T>::extend()
{
  assert(extend_size);
  size_t new_size = this -> size + this -> extend_size;
  ARRAY_REALLOC(data, new_size, T);
  this -> size = new_size;
}


