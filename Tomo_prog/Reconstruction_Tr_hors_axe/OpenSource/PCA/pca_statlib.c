/* ***********************  Contents  ****************************************
   Principal Components Analysis: C, 638 lines. ****************************
   * Sample input data set (final 36 lines). *********************************
   ***************************************************************************/



/*********************************/
/* Principal Components Analysis */
/*********************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h> 

#include "pca_statlib.h"

#define SIGN(a, b) ( (b) < 0 ? -fabs(a) : fabs(a) )


/* internal functions prototypes */
void tqli(double d[], double e[], int n, double **z);
void tred2(double **a, int n, double *d, double *e);
double **matrix(int n, int m);
double *vector(int n);
void free_matrix(double **mat, int n, int m);
void free_vector(double* v, int n);
void erhand(char err_msg[]);
void scpcol(double **data, int n, int m, double **symmat);
void covcol(double **data, int n, int m, double **symmat);
void corcol(double **data, int n, int m, double **symmat);


int
covariance_pca(int nb_instances, int nb_components, double** tset_data, 
	       double* eigenvals, double** eigenvectors)
{
  FILE *stream;
  int  n, m,  i, j, k, k2;
  double **data, **symmat, **symmat2, *evals, *interm;

  double in_value;
  char option;
  

  /*********************************************************************
   Get from command line:
   input data file name, #rows, #cols, option.

   Open input file: fopen opens the file whose name is stored in the
   pointer argv[argc-1]; if unsuccessful, error message is printed to
   stderr.
  *********************************************************************/

  n = nb_instances;              /* # rows */
  m = nb_components;        /* # columns */

  /** on ne les laisse plus calculer la covariance, ils savent pas faire
   // n: rows, m: columns 
   data = tset_data;
   // matrix(n, m);  // Storage allocation for input data 

   symmat = matrix(m, m);  // Allocation of correlation (etc.) matrix 

   // Look at analysis option; branch in accordance with this. 

   //   printf("Analysis of variances-covariances chosen.\n"); 
   covcol(data, n, m, symmat);
   /// computation of covariance matrix from the training set 
   */

  //la matrice de covariance est deja calculee
  // copions la dans eigenvectors, 
  //symmat = matrix(m, m);
  for (i = 1; i <= m; i++)
    for (j = 1; j <= m; j++)
      eigenvectors[i][j] = tset_data[i][j];


  /*********************************************************************
    Eigen-reduction
  **********************************************************************/

  /* Allocate storage for dummy and new vectors. */
  evals = vector(m);     /* Storage alloc. for vector of eigenvalues */
  interm = vector(m);    /* Storage alloc. for 'intermediate' vector */


  /** computation of eigenvalues/vectors */
  tred2(eigenvectors, m, evals, interm);  /* Triangular decomposition */
  tqli(evals, interm, m, eigenvectors);   /* Reduction of sym. trid. matrix */
  /* evals now contains the eigenvalues,
     columns of symmat now contain the associated eigenvectors. */
  

  /** reordering in ascending order */
  /* if (*eigenvals) 
     free_vector(*eigenvals, m); */
  //*eigenvals = vector(m);     
  for (i = 1; i <= m; i++)
    {eigenvals[i] = evals[m - i + 1];}
  
  /* free(*eigenvectors); */
  //*eigenvectors = symmat;
  

  /* etait deja */
  /*   free_matrix(symmat, m, m); */
  /*   free_matrix(data, n, m);   */

  free_vector(evals, m);  
  free_vector(interm, m);  

  return L_TRUE;
}


/**  Variance-covariance matrix: creation  *****************************/



/** ca ne marche pas, ces abrutis!  */ 
/* Create m * m covariance matrix from given n * m data matrix. */
void covcol(double **data, int n, int m, double **symmat)
{
  double *mean;
  int i, j, j1, j2;

  /* Allocate storage for mean vector */

  mean = vector(m);

  /* Determine mean of column vectors of input data matrix */

  for (j = 1; j <= m; j++)
    {
      mean[j] = 0.0;
      for (i = 1; i <= n; i++)
        {
	  mean[j] += data[i][j];
        }
      mean[j] /= (double)n;
    }

  /*
    printf("\nMeans of column vectors:\n");
    for (j = 1; j <= m; j++)  {
    printf("%7.1f",mean[j]);  }   printf("\n");
  */

  /* Center the column vectors. */

  /* not in vtk either */
  for (i = 1; i <= n; i++)
    {
      for (j = 1; j <= m; j++)
	{
	  data[i][j] -= mean[j];
	}
    }


  /* Calculate the m * m covariance matrix. */
  for (j1 = 1; j1 <= m; j1++)
    {
      for (j2 = j1; j2 <= m; j2++)
        {
	  symmat[j1][j2] = 0.0;
	  for (i = 1; i <= n; i++)
            {
	      symmat[j1][j2] += data[i][j1] * data[i][j2];
            }
	  //symmat[j2][j1] /= (m - 1); /** comme chez cootes, et pas dans pca statlib de base */
	  symmat[j2][j1] = symmat[j1][j2];
        }
    }

  return;

}



/**  Error handler  **************************************************/


/* Error handler */
void erhand(char err_msg[])
{
  fprintf(stderr,"Run-time error:\n");
  fprintf(stderr,"%s\n", err_msg);
  fprintf(stderr,"Exiting to system.\n");
  exit(1);
}

/**  Allocation of vector storage  ***********************************/


/* Allocates a double vector with range [1..n]. */
double *vector(int n)
{

  double *v;

  v = (double *) malloc ((unsigned) n*sizeof(double));
  if (!v) erhand("Allocation failure in vector().");
  return v-1;

}

/**  Allocation of double matrix storage  *****************************/


/* Allocate a double matrix with range [1..n][1..m]. */
double **matrix(int n, int m)
{
  int i;
  double **mat;

  /* Allocate pointers to rows. */
  mat = (double **) malloc((unsigned) (n)*sizeof(double*));
  if (!mat) erhand("Allocation failure 1 in matrix().");
  /*   assert(mat); */
  
  mat -= 1;

  /* Allocate rows and set pointers to them. */
  for (i = 1; i <= n; i++)
    {
      mat[i] = (double *) malloc((unsigned) (m)*sizeof(double));
      if (!mat[i]) erhand("Allocation failure 2 in matrix().");
      /* assert(mat[i]); */

      mat[i] -= 1;
    }

  /* Return pointer to array of pointers to rows. */
  return mat;

}

/**  Deallocate vector storage  *********************************/


/* Free a double vector allocated by vector(). */
void free_vector(double* v, int n)
{
  free((char*) (v+1));
}

/**  Deallocate double matrix storage  ***************************/


/* Free a double matrix allocated by matrix(). */
void free_matrix(double **mat, int n, int m)
{
  int i;

  for (i = n; i >= 1; i--)
    {
      free ((char*) (mat[i]+1));
    }
  free ((char*) (mat+1));
}

/**  Reduce a real, symmetric matrix to a symmetric, tridiag. matrix. */


/* Householder reduction of matrix a to tridiagonal form.
   Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
   Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
   Springer-Verlag, 1976, pp. 489-494.
   W H Press et al., Numerical Recipes in C, Cambridge U P,
   1988, pp. 373-374.  */
void tred2(double **a, int n, double *d, double *e)
{
  int l, k, j, i;
  double scale, hh, h, g, f;

  for (i = n; i >= 2; i--)
    {
      l = i - 1;
      h = scale = 0.0;
      if (l > 1)
	{
	  for (k = 1; k <= l; k++)
	    scale += fabs(a[i][k]);
	  if (scale == 0.0)
	    e[i] = a[i][l];
	  else
	    {
	      for (k = 1; k <= l; k++)
		{
		  a[i][k] /= scale;
		  h += a[i][k] * a[i][k];
		}
	      f = a[i][l];
	      g = f>0 ? -sqrt(h) : sqrt(h);
	      e[i] = scale * g;
	      h -= f * g;
	      a[i][l] = f - g;
	      f = 0.0;
	      for (j = 1; j <= l; j++)
		{
		  a[j][i] = a[i][j]/h;
		  g = 0.0;
		  for (k = 1; k <= j; k++)
		    g += a[j][k] * a[i][k];
		  for (k = j+1; k <= l; k++)
		    g += a[k][j] * a[i][k];
		  e[j] = g / h;
		  f += e[j] * a[i][j];
		}
	      hh = f / (h + h);
	      for (j = 1; j <= l; j++)
		{
		  f = a[i][j];
		  e[j] = g = e[j] - hh * f;
		  for (k = 1; k <= j; k++)
		    a[j][k] -= (f * e[k] + g * a[i][k]);
		}
	    }
	}
      else
        e[i] = a[i][l];
      d[i] = h;
    }
  d[1] = 0.0;
  e[1] = 0.0;
  for (i = 1; i <= n; i++)
    {
      l = i - 1;
      if (d[i])
	{
	  for (j = 1; j <= l; j++)
	    {
	      g = 0.0;
	      for (k = 1; k <= l; k++)
		g += a[i][k] * a[k][j];
	      for (k = 1; k <= l; k++)
		a[k][j] -= g * a[k][i];
	    }
	}
      d[i] = a[i][i];
      a[i][i] = 1.0;
      for (j = 1; j <= l; j++)
	a[j][i] = a[i][j] = 0.0;
    }
}

/**  Tridiagonal QL algorithm -- Implicit  **********************/

void tqli(double d[], double e[], int n, double **z)
{
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;


  for (i = 2; i <= n; i++)
    e[i-1] = e[i];
  e[n] = 0.0;
  for (l = 1; l <= n; l++)
    {
      iter = 0;
      do
	{
	  for (m = l; m <= n-1; m++)
	    {
	      dd = fabs(d[m]) + fabs(d[m+1]);
	      if (fabs(e[m]) + dd == dd) break;
	    }
          if (m != l)
	    {
	      if (iter++ == 30) erhand("No convergence in TLQI.");
	      g = (d[l+1] - d[l]) / (2.0 * e[l]);
	      r = sqrt((g * g) + 1.0);
	      g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
	      s = c = 1.0;
	      p = 0.0;
	      for (i = m-1; i >= l; i--)
		{
		  f = s * e[i];
		  b = c * e[i];
		  if (fabs(f) >= fabs(g))
                    {
		      c = g / f;
		      r = sqrt((c * c) + 1.0);
		      e[i+1] = f * r;
		      c *= (s = 1.0/r);
                    }
		  else
                    {
		      s = f / g;
		      r = sqrt((s * s) + 1.0);
		      e[i+1] = g * r;
		      s *= (c = 1.0/r);
                    }
		  g = d[i+1] - p;
		  r = (d[i] - g) * s + 2.0 * c * b;
		  p = s * r;
		  d[i+1] = g + p;
		  g = c * r - b;
		  for (k = 1; k <= n; k++)
		    {
		      f = z[k][i+1];
		      z[k][i+1] = s * z[k][i] + c * f;
		      z[k][i] = c * z[k][i] - s * f;
		    }
		}
	      d[l] = d[l] - p;
	      e[l] = g;
	      e[m] = 0.0;
	    }
	}  while (m != l);
    }
}
