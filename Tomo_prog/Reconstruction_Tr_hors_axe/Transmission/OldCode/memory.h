#ifndef __MEMORY__
#define __MEMORY__


void zerofill_double(double *array, int card);


// pour circhifts à décalage arbitraire uniquement. Trop lent
void
circshift_ccd2(double* entree, double* entree_shifted, int dimx, int dimy, int decal_x,int decal_y);

// circshift à décalage centré. On ne précise que la taille de l'image source (qui doit = dest)
// basé sur opencv pour effectuer des copies mémoire rapides
void
circshift_2D_cv(double* entree, size_t dimx, size_t dimy);

void
circshift_2D_cv(double* entree, double* sortie, size_t dimx, size_t dimy);

void 
circshift_2D_mc(double *image2D, double *image2D_shift, int dim_x, int dim_y);


// effectue un décalage circulaire du volume chargé en mémoire vers la destination spécifiée
// le cube est subdivisé en 8 cubes égaux et le contenu de chaque sous-cube est permuté 
// le volume est implémenté par un tableau séquentiel et la copie est effectuée
void 
circshift3D_memcpy(double *volume3D, double *volume3D_shift, int dim_x, int dim_y, int dim_z);
// déprécié si AIR_Volume/OpenCV présents


void
convert_double(INPUT_TYPE* src, double* dst, int cardinal);

void
convert_float(double* src, OUTPUT_TYPE* dst, int cardinal);

int
nonzero_d(double* src, int cardinal);

int
nonzero_p(INPUT_TYPE* src, int cardinal);

#endif // __MEMORY__



// void circshift2D_memcpy(double *image2D, double *image2D_shift, int dim_x, int dim_y);
// void extract_subImage(double* src, double* dst, int src_dimx, int src_dimy, int edge_x, int edge_y, int dst_dimx, int dst_dimy);
