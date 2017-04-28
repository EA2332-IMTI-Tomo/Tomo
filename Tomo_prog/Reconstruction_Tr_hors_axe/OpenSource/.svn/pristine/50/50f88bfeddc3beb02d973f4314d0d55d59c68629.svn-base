#ifndef _PCA_
#define _PCA_


/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2001 the GLOS development team (http://glos.sourceforge.net) */
// underlying code provided by other sources


#include <iostream>
#include <string>

/* #include <fstream> */
#include "macros.h"


#define DEBUG 1

using namespace std;


template <class T>
class PCA{


 protected:

  ostream& display(ostream &o) const;
  ostream& display_nice(ostream &o) const;

  int vector_dimension, vectors_count;

  bool training_data_set, symsquare_data_set; // ie sym-square data
  bool dimension_set, vectors_count_set;
  bool pca_performed;
  bool training_data_centered;
  bool normalized_values;

  // matrix of original vectors training set.
  // Expressed in _statlib_ format! (mat[1][1] -> mat[n][n])
  double** training_data; 
  double** symsquare_matrix; // symmetric square matrix that can be derived
			// (generally by covariance) from training set
			// or specified directly.
  

  double* eigenvalues;
  double** eigenvectors;
  

 private:

  

  
  /* these methods transform matrices from standard <-> statlib
     (shifted) format */
  double** normal_to_shifted(double** matrix, int dim_x);
  double** shifted_to_normal(double** matrix, int dim_x);


  void constructor()
    {
      training_data_set = symsquare_data_set = dimension_set = normalized_values
	= vectors_count_set = pca_performed = training_data_centered = false;
    }

  void reset()
    {
      if (training_data_set) MATRIX_FREE(training_data, vectors_count + 1); // bastard-consed
      if (symsquare_data_set) MATRIX_FREE(symsquare_matrix, vector_dimension - 1);
      pca_performed = symsquare_data_set = false;
    }


 public:
  

  PCA()
    {constructor();}    
  PCA(int vectors_count, int vector_dimension)
    {
      constructor(); 
      this -> vector_dimension = vector_dimension; 
      this -> vectors_count = vectors_count;
      dimension_set = vectors_count_set = true;
    }
  ~PCA() {reset();} 
  
  
  READER(int, vector_dimension);
  READER(int, vectors_count);
  READER(double**, symsquare_matrix);
  READER(double**, training_data);


  void set_dimension(int vector_dimension) {
    this -> vector_dimension = vector_dimension; dimension_set = true;
    reset();
  }
  
  void set_instances_num(int vectors_count) {
    this -> vectors_count = vectors_count; vectors_count_set = true;
    reset();
  }
    
  
  // symsquare_data_matrix = covariance matrix (in the most general case)
  void set_symsquare_data_matrix(double** matrix);

  // compute covariance matrix from training set. Ignored by perform_covar_pca.
  bool compute_data_covar_matrix();
  
  
  /* attach the training data */
  void set_training_data(double** data_matrix); // matrix in standard form: double[nb_ins][dim]
  void set_training_data(char* data_matrix_file); // file: one instance per line
  
  /* computes PCA from covariance matrix inferred from training
     data. symsquare_matrix is ignored!!. */
  bool perform_covar_pca();
  /* relies on statlib code */
  
  /* compute PCA from given symsquare_matrix */
  /*   bool perform_pca();  */
  /* relies on ? code*/

  
  double nth_eigenvalue(int n) const;
  double* nth_eigenvector(int n) const;
  
  //once called, eigenvalues and eigenvectors are normalized
  //(eigenvals sum is one, such as eigenvector components sum)
  void normalize_values();
  
  friend ostream& operator<< <T> (ostream &o, const PCA<T> &p);
};



#include "PCA/PCA.cc"


#endif /* _PCA_ */
