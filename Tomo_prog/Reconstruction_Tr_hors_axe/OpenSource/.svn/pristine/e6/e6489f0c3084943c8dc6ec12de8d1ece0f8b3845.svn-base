/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2001 the GLOS development team (http://glos.vectraproject.com) */


#ifndef _VECTRA_MATRIX_CC_
#define _VECTRA_MATRIX_CC_



template <class T>
vMatrix<T>::vMatrix(size_t nb_rows, size_t nb_cols, T initval)
{
  this -> nb_rows = nb_rows;
  this -> nb_cols= nb_cols;
  
  this -> display_limited = false;
  this -> display_limit = 30;

  MATRIX_ALLOC(data, nb_rows, nb_cols, T);
  
  this -> fill(initval);
}


template <class T>
vMatrix<T>::vMatrix(const vMatrix<T> &m)
{
  this -> nb_rows = m.nb_rows;
  this -> nb_cols= m.nb_cols;  
  this -> display_limited = m.display_limited;
  this -> display_limit = m.display_limit;
  
  MATRIX_ALLOC(data, nb_rows, nb_cols, T);
  
  size_t i, j;
  for (i = 0; i < nb_rows; i++)
    for (j = 0; j < nb_cols; j++)  // consider memcpy instead
      this -> data[i][j] = m.data[i][j];
}



template <class T>
vMatrix<T>::vMatrix(size_t nb_rows, size_t nb_cols, T** data)
{
  this -> nb_rows = nb_rows;
  this -> nb_cols= nb_cols;

  this -> display_limited = false;
  this -> display_limit = 30;

  this -> data = data;
}


template <class T>
vMatrix<T>::~vMatrix()
{
  MATRIX_FREE(data, nb_rows);
}


template <class T>
void
vMatrix<T>::fill(T value)
{
  size_t i, j;
  for (i = 0; i < nb_rows; i++)
    for (j = 0; j < nb_cols; j++)
      data[i][j] = value;
}


template <class T>
void
vMatrix<T>::scalar_mult(T value)
{
  size_t i, j;
  for (i = 0; i < nb_rows; i++)
    for (j = 0; j < nb_cols; j++)
      data[i][j] *= value;
}


template <class T>
void 
vMatrix<T>::map(T (*map_fun)(T))
{
  size_t i, j;
  for (i = 0; i < nb_rows; i++)
    for (j = 0; j < nb_cols; j++)
      data[i][j] = map_fun(data[i][j]);
}


template <class T>
void 
vMatrix<T>::multiply(const vMatrix<T> &m1, const vMatrix<T> &m2)
{
  size_t i, j, k;
  size_t m1r = m1.nb_rows, m2r = m2.nb_rows;
  size_t m1c = m1.nb_cols, m2c = m2.nb_cols;
 
  assert ( m1c == m2r &&			\
	   m1r == nb_rows && m2c == nb_cols);
  
  long double acc = 0;
  size_t vec_l = m1.nb_cols;

  for (i = 0; i < nb_rows; i++)
    for (j = 0; j < nb_cols; j++)
      {
	acc = 0;
	
	for (k = 0; k < vec_l ; k++)
	  acc += ((long double)m1[i][k] * (long double)m2[k][j]);
	
	data[i][j] = (T)acc;
      }
}


template <class T>
void 
vMatrix<T>::multiply_from_diag(const vMatrix<T> &vector, const vMatrix<T> &m)
{
  size_t i, j;

  assert (vector.nb_rows == 1);

  size_t m1c = vector.nb_cols;
  size_t m2r = m.nb_rows, m2c = m.nb_cols;
 
  assert ( m1c == m2r &&			\
	   m2r == nb_rows && m2c == nb_cols);
  
  for (i = 0; i < nb_rows; i++)
    for (j = 0; j < nb_cols; j++)
      data[i][j] = m[i][j] * vector[0][i];
}


template <class T>
void 
vMatrix<T>::transpose(const vMatrix<T> &m)
{
  assert( nb_rows == m.nb_cols && m.nb_rows == nb_cols);
  
  size_t i = 0, li =  m.nb_rows;
  size_t j = 0, lj =  m.nb_cols;

  for (i = 0; i < li; i++)
    for (j = 0; j < lj; j++)
      data[j][i] = m.data[i][j];
}


template <class T>
void
vMatrix<T>::display(ostream &o) const
{
  size_t i, j;
  for (i = 0; i < nb_rows; i++)
    {
      if (display_limited && (i > display_limit)) 
	{o << endl << "..." << endl; break;}
      o << endl;
      for (j = 0; j < nb_cols; j++) 
	if (display_limited && (j > display_limit)) 
	  {o << "...."; break;}
	else
	  {o << data[i][j] << " ";}
      //o << "|";
    }
}


// operateur de recopie: m = n (ssi mêmes tailles)
template <class T>
void
vMatrix<T>::operator= (const vMatrix<T> &m)
{
  assert(this -> nb_rows == m.nb_rows && \
	 this -> nb_cols == m.nb_cols);
  
  size_t i, j;
  for (i = 0; i < nb_rows; i++)
    for (j = 0; j < nb_cols; j++)  // consider memcpy instead
      this -> data[i][j] = m.data[i][j];
}


// =============================================================================


// friend function: cout = (cout << m)
template <class T>
ostream& 
operator<< (ostream& o, const vMatrix<T> &m)
{
  m.display(o);
  return o;
}


// =============================================================================


#endif /* _VECTRA_MATRIX_CC_ */





  // advanced features: read matrix from text file. Delimiter: space. One row per line
  // vMatrix(char* filename, size_t max_line_size = 1000, size_t offset_rows = 0, size_t offset_cols = 0);


/*
template <class T>
vMatrix<T>::vMatrix(char* filename, size_t max_line_size, size_t offset_rows, size_t offset_cols)
{
  ifstream infile (filename);
  assert (infile.is_open());
  
  char* read_line;
  ARRAY_ALLOC(read_line, max_line_size, char);
  
  
}
*/


