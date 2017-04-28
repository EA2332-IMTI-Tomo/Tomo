#include "pca_statlib.h"


template <class T>
double**
PCA<T>::normal_to_shifted(double** matrix, int dim_x)
{
  for (int i = 0; i < dim_x; i++)
    matrix[i] -= 1;

  return --matrix; //matrix - 1
}


template <class T>
double**
PCA<T>::shifted_to_normal(double** matrix, int dim_x)
{
  matrix++;

  for (int i = 0; i < dim_x; i++)
    matrix[i] += 1;

  return matrix; //matrix
}


template <class T>
void
PCA<T>::set_training_data(double** data_matrix)
{
  assert (vector_dimension);
  assert (vectors_count);
  reset();

  MATRIX_ALLOC(training_data, vectors_count, vector_dimension, double);

  int i, j;
  for (i = 0; i < vectors_count; i++)
    for (j = 0; j < vector_dimension; j++)
      training_data[i][j] = data_matrix[i][j];
      
  training_data_set = true;
}


template <class T>
void
PCA<T>::set_training_data(char* data_matrix_file)
{
  int i, j;
  T lu;
  FILE* fp;

  assert (vector_dimension);
  assert (vectors_count);
  
  if (training_data_set) reset();
  
  
  MATRIX_ALLOC(training_data, vectors_count, vector_dimension, double);
  MATRIX_FILE_TXT(training_data, vectors_count, vector_dimension, data_matrix_file, "lf");

  
  training_data_set = true;
}


template <class T>
void
PCA<T>::set_symsquare_data_matrix(double** matrix)
{
  symsquare_matrix = matrix;
  symsquare_data_set = true;
}


template <class T>
bool
PCA<T>::perform_covar_pca()
{
  bool success = false;
  int i, j;
  double** l_col_eigenvectors;

  assert(training_data_set);

  /* for some reason, if done just after PCA, erases eigenvalues! */
  // if (pca_performed) MATRIX_FREE(eigenvectors, vectors_count);
  MATRIX_ALLOC(eigenvectors, vector_dimension, vector_dimension, double);
  MATRIX_ALLOC(l_col_eigenvectors, vector_dimension, vector_dimension, double);
  ARRAY_ALLOC(eigenvalues, vector_dimension, double);


  if (! symsquare_data_set)
    assert (compute_data_covar_matrix());
  
      
  /* performs to statlib PCA. Data is transformed before and after
     since statlib uses shifted matrices and arrays. */
  // success = covariance_pca(vectors_count, vector_dimension, training_data, &eigenvalues, &col_eigenvectors);
  symsquare_matrix = normal_to_shifted(symsquare_matrix, vector_dimension);
  l_col_eigenvectors = normal_to_shifted(l_col_eigenvectors, vector_dimension);
  eigenvalues--;
  success = covariance_pca(vectors_count, vector_dimension, symsquare_matrix, eigenvalues, l_col_eigenvectors);
  symsquare_matrix = shifted_to_normal(symsquare_matrix, vector_dimension);
  l_col_eigenvectors = shifted_to_normal(l_col_eigenvectors, vector_dimension);
  eigenvalues++; /* shifted array! */

  if (success)
    {
      for (i = 0; i < vector_dimension; i++)
	for (j = 0; j < vector_dimension; j++)
	  eigenvectors[i][j] = l_col_eigenvectors[j][vector_dimension - 1 - i];
    }

  /*
    if (success)
    {
    int nb_vectors, nb_components;
    nb_components = nb_vectors = vector_dimension;

    symsquare_data_set = true;
      
    // "decodes" content of allocated eigenvalues/vectors into "civilized" form 
    for (i = nb_vectors; i > 0; i--) 
    for (j = 1; j <= nb_components; j++) 
    eigenvectors[vector_dimension - i][j - 1] = col_eigenvectors[j][i];
    }
  */

  MATRIX_FREE(l_col_eigenvectors, vector_dimension);
  return (pca_performed = success);
}


template <class T>
double
PCA<T>::nth_eigenvalue(int n) const
{
  assert(pca_performed);
  assert(BETWEEN(n, 0, vector_dimension - 1));
  return(eigenvalues[n]);  
}


template <class T>
double*
PCA<T>::nth_eigenvector(int n) const
{
  assert(pca_performed);
  assert(BETWEEN(n, 0, vector_dimension - 1));
  return(eigenvectors[n]);  
}


template <class T>
ostream& 
operator<<(ostream &o, const PCA<T> &p) 
{
  p.display(o);
  return o;
}


template <class T>
ostream&
PCA<T>::display(ostream &o) const
{
  int i, j;
  double* eigenvec;

  double eigensum = 0;
  for (i = 0; i < vector_dimension; i++)
    eigensum += nth_eigenvalue(i);


  o << endl <<  "#Eigenvalues (descending order):" << endl;
  for (i = 0; i < vector_dimension; i++)
    o << nth_eigenvalue(i) << " (" << (100 * nth_eigenvalue(i)) / eigensum << " %) | "  ; 
  
  o << endl << "#Eigenvectors (corresponding):";
  for (i = 0; i < vector_dimension; i++)
    {
      o << endl << "vector " << i << ": ";
      eigenvec = nth_eigenvector(i);

      for (j = 0; j < vector_dimension; j++)
	o << eigenvec[j] << " ";
    }
  

  return o;
}


template <class T>
void
PCA<T>::normalize_values()
{
  int i, j;
  double eval_sum = 0;
  double comp_sum = 0;

  for (i = 0; i < vector_dimension; i++)
    eval_sum +=  eigenvalues[i];
  
  for (i = 0; i < vector_dimension; i++)
    eigenvalues[i] /= eval_sum;

  for (i = 0; i < vector_dimension; i++)
    {
      for (comp_sum = 0, j = 0; j < vector_dimension; j++)
	comp_sum += eigenvectors[i][j];
      for (j = 0; j < vector_dimension; j++)
	eigenvectors[i][j] /= comp_sum;
    }
}



/* traduction de methode Cootes */
template <class T>
bool
PCA<T>::compute_data_covar_matrix()
{
  assert(training_data_set);
  
  double* dif_xmean;
  double* xmean;
  int i, j, j1, j2;


  MATRIX_ALLOC(symsquare_matrix, vector_dimension, vector_dimension, double);
  ARRAY_ALLOC(dif_xmean, vector_dimension, double); //dx_i = x_i - xm
  ARRAY_ALLOC(xmean, vector_dimension, double); // xm 

  
  /* compute mean shape instance */
  for (j = 0; j < vector_dimension; j++)
    {
      xmean[j] = 0.0;
      for (i = 0; i < vectors_count; i++)
	xmean[j] += training_data[i][j];
	
      xmean[j] /= (double)vectors_count;
    }
  
  
  int ii;
  //N: nb instances, S: final covariance matrix
  // S = (1/N) sum_{i=1}_{N}
  for (i = 0; i < vectors_count; i++)
    {
      //ARRAY_SET(dif_xmean, vector_dimension, 0.0);
      
      for (j = 0; j < vector_dimension; j++)
	dif_xmean[j] = training_data[i][j] - xmean[j];
      
      // S += dx_i * dx_i _T
      // matrice symetrique: on calcule la partie inf gauche
      for (ii = 0; ii < vector_dimension; ii++)
	for (j = 0; j <= ii; j++)
	  symsquare_matrix[ii][j] += (dif_xmean[ii] * dif_xmean[j]);
    }
  
  
  //fin: on divise par N et on remplit la partie droite
  double val;
  for (i = 0; i < vector_dimension; i++)
    for (j = 0; j <= i; j++)
      {
	symsquare_matrix[i][j] /= vectors_count;
	symsquare_matrix[j][i] = symsquare_matrix[i][j]; 
      }
  // S = (1/N) sum_{i=1}_{N} (dx_i * dx_i _T) implemente
  
  
  if (DEBUG > 0)
    {
      for (j1 = 0; j1 < vector_dimension; j1++)
	{
	  cout << endl;
	  for (j2 = 0; j2 < vector_dimension; j2++)
	    cout << "| " << symsquare_matrix[j1][j2];
	} 
    }


  free(xmean);
  free(dif_xmean);
  return (symsquare_data_set = true);
}



/* traduction de pca_statlib/covcol avec des indices civilises. 
   mauvaise idee, ca marche tres mal (l'original) */
/*
  template <class T>
  bool
  PCA<T>::compute_data_covar_matrix()
  {
  assert(training_data_set);
  
  double *col_mean;
  int i, j, j1, j2;


  MATRIX_ALLOC(symsquare_matrix, vector_dimension, vector_dimension, double);
  ARRAY_ALLOC(col_mean, vector_dimension, double);

  
  // compute mean of columns (mean value of component over training set)
  for (j = 0; j < vector_dimension; j++)
  {
  col_mean[j] = 0.0;
  for (i = 0; i < vectors_count; i++)
  col_mean[j] = col_mean[j] + training_data[i][j];
      
  col_mean[j] /= (double)vectors_count;
  }
  

  // Center the column vectors - avoid further redundant computations
  // for covariance matrix. Does not change anything from stat. point
  // of view, but alters original data. 
  for (i = 0; i < vectors_count; i++)
  for (j = 0; j < vector_dimension; j++)
  training_data[i][j] -= col_mean[j];

  training_data_centered = true;
  

  // Calculate the square, symetric covariance matrix. 
  for (j1 = 0; j1 < vector_dimension; j1++)
  for (j2 = j1; j2 < vector_dimension; j2++)
  {
  symsquare_matrix[j1][j2] = 0.0;
  for (i = 0; i < vectors_count; i++)
  symsquare_matrix[j1][j2] += training_data[i][j1] * training_data[i][j2];
	
  symsquare_matrix[j2][j1] = symsquare_matrix[j1][j2];
  }
  // cov[in][d] = sum_[k=0,in] (data[k][d] - col_mean[d]) 
  //  *  (data[k][in] - col_mean[in]
  // in: instances number, d: dimension 


  if (DEBUG > 0)
  {
  for (j1 = 0; j1 < vector_dimension; j1++)
  {
  cout << endl;
  for (j2 = 0; j2 < vector_dimension; j2++)
  cout << "| " << symsquare_matrix[j1][j2];
  } 
  }


  free(col_mean);
  return (symsquare_data_set = true);
  }
*/

