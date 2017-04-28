#ifndef __PCA_STATLIB__
#define __PCA_STATLIB__


/** This library can be considered as GPL? */

/** minor modifications (rough code to std gnu c-library) by
    J. Bailleul (bailleul@greyc.ismra.fr) */


  /*********************************************************************/
  /* Principal Components Analysis or the Karhunen-Loeve expansion is a
     classical method for dimensionality reduction or exploratory data
     analysis.  One reference among many is: F. Murtagh and A. Heck,
     Multivariate Data Analysis, Kluwer Academic, Dordrecht, 1987.

     Author:
     F. Murtagh
     Phone:        + 49 89 32006298 (work)
     + 49 89 965307 (home)
     Earn/Bitnet:  fionn@dgaeso51,  fim@dgaipp1s,  murtagh@stsci
     Span:         esomc1::fionn
     Internet:     murtagh@scivax.stsci.edu
   
     F. Murtagh, Munich, 6 June 1989                                   */   
  /*********************************************************************/


/* #include "current_type.h" */
#define L_TRUE 1

/* input: nb_instances, nb_components (number of components of
   each data instance), and tset_data of size nb_ins * inst_dim
   eigenvals and eigenvectors should be void pointers.

   output: eigenvals was allocated and filled with nb_components
   eigenvalues in descending orders.

   eigenvectors was allocated and filled with corresponding
   nb_components eigevectors (of size nb_components).

   -> now: eigenvalues and eigenvectors must be allocated (as shifted) before call
   will be set to desired values.


   returns EXIT_SUCCESS.
*/
int
covariance_pca(int nb_instances, int nb_components, double** tset_data, 
	       double* eigenvals, double** eigenvectors);



/** */
/*   printf("\nEigenvectors:\n"); */
/*
  printf("(First three; their definition in terms of original vbes.)\n");
  for (j = 1; j <= m; j++) {
  for (i = 1; i <= 3; i++)  {
  printf("%12.4f", symmat[j][m-i+1]);  }
  printf("\n");  }
*/


#endif /* __PCA_STATLIB__ */
