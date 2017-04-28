#ifndef __CU_HOLO_PROCESS_KERNELS
#define __CU_HOLO_PROCESS_KERNELS


#include <iostream>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cutil_inline.h>

#include "macros.h"

using namespace std;




// tableau de blocs consécutifs
#define CUDAINDEXE				      \
  blockIdx.y * ( gridDim.x * blockDim.x ) +	      \
  blockIdx.x * blockDim.x +			      \
  threadIdx.x;					      \



#define CUDASYNCE				\
  cudaThreadSynchronize();			\


__global__ void
kernel_mapping_float(int l_xm0, int l_ym0, size_t p_Nxmax, size_t p_Nymax, \
		     float p_rayon, float p_delta_zmax, double zm0,	\
		     size_t c_dv0s2rf, size_t c_dv1s2rf, size_t c_dv0rf, size_t c_dv2s2rf, size_t c_dv1xdv0rf, \
		     unsigned short int *sup_redon_d, float *reel_arc_d, float *imag_arc_d, \
		     double *fft_reel_shift_norm_d, double *fft_imag_shift_norm_d);


__global__ void
kernel_crop2D_double(double *src_data, double *dst_data, size_t src_size, size_t dst_size);



__global__ void
kernel_norm_double(double *dst_data, double *src_data_x, double *src_data_y, size_t src_size);


__global__ void
kernel_normpeak_double(double *dst_data_r, double *dst_data_i, double *src_data_r, double *src_data_i, size_t src_size, double max_part_reel, double max_part_imag, double inv_max_mod);


__global__ void
kernel_applymask_double(double* data_r, double* data_i, size_t data_size_edge, \
			unsigned char* data_mask, size_t data_mask_edge, \
			size_t spot_x, size_t spot_y);


#endif //__CU_HOLO_PROCESS_KERNELS
