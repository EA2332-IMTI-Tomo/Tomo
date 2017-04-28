#ifndef __VOLUME__COPY__
#define __VOLUME__COPY__


// TOUTE CETTE SECTION EST PLUS OU MOINS A FOUCHTRE A LA POUBELLE
// Se référer à AIR_Volume


// works for src > dst, src and dst even, and of course src and dst isocubes (we give edge lengthes)
template <typename T>
void
copy_3Dmatrix_in_bigger(T* src, T* dst, size_t src_edge, size_t dst_edge);


// conversely,
template <typename T>
void
copy_3Dmatrix_in_smaller(T* src, T* dst, size_t src_edge, size_t dst_edge);



// same with effective functions (2 volumes vs C2C structure)
/* #include <cuda.h> */
/* #include <cuda_runtime.h> */
#include <cufft.h>


void
copy_3Dmatrix_in_biggerC2C(double* srcR, double* srcI, cufftComplex* dst, size_t src_edge, size_t dst_edge);


void
copy_3DmatrixC2C_in_smaller(cufftComplex* big, double* smallR, double* smallI, size_t big_edge, size_t small_edge);




#endif // __VOLUME__COPY__
