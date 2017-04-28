#ifndef __REFOCUS__
#define __REFOCUS__


///variable pour synthèse ouverture 2D
//on agrandit le support en doublant largeur et longueur (facteur 4)
//Var2D dimInit={2*NXMAX,2*NYMAX},dimFin={4*NXMAX,4*NYMAX};

typedef struct {
  int x,y;
}Var2D;

typedef struct {
  int x,y,z;
}Var3D;

typedef struct {
  double Re,Im;
}nbCplx;


void AutofocusGrad(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe, double *sortieIm, size_t dimX, size_t dimY, double rayon, size_t SpecPosX, size_t SpecPosY, int cpt_angle);

void 
AutofocusGrad(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe, double *sortieIm, Var2D dimInit, double rayon, Var2D CoordSpecI, int cpt_angle);


/*
// retourne le nombre de pts exclus
int 
retroPropag(double* spectre3D_Re, double* spectre3D_Im, double *sup_redon, \
	    int dim_final, double* Spectre2D_Re, double* Spectre2D_Im, \
	    Var2D posSpec, Var3D decal,Var2D NMAX, double rayon);
*/

void  
retroPropagSA(float deltaZ, double *fft_shift_normRe, double *fft_shift_normIm, \
	      double* spectrePropag_normRe, double* spectrePropag_normIm, \
	      Var2D dimSpctHolo, double rayon, Var2D CoordSpecI);


void 
Grad(double* entreeRe, double* entreeIm, double* Gradient, int taille, int NbPlan, int cptPlan);



/*
retroPropagSA
Plan_ds_Vol
Grad
Var2D Var3D
*/


#endif // __REFOCUS__
