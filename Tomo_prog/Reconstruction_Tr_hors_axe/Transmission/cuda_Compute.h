#ifndef __CUDA_COMPUTE__
#define __CUDA_COMPUTE__

void
compute_gpu_online_batch(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
			 size_t window_edge_x, size_t window_edge_y, size_t image_dim_x, size_t image_dim_y, \
			 size_t xm0_limite, size_t ym0_limite, float rayon, float delta_zmax, \
			 size_t batch_size, \
			 const string &images_radix);



#endif //__CUDA_COMPUTE__
