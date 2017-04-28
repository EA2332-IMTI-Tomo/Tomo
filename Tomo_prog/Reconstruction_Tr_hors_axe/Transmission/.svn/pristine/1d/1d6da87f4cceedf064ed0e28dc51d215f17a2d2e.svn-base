#include "cu_Holo_Process_kernels.h"



// *****************************************************************************
// 3D mapping 
// *****************************************************************************


// -----------------------------------------------------------------------------
// kernel CUDA


__global__ void
kernel_mapping_float(int l_xm0, int l_ym0, size_t p_Nxmax, size_t p_Nymax, \
		     float p_rayon, float p_delta_zmax, double zm0,	\
		     size_t c_dv0s2rf, size_t c_dv1s2rf, size_t c_dv0rf, size_t c_dv2s2rf, size_t c_dv1xdv0rf, \
		     unsigned short int *sup_redon_d, float *reel_arc_d, float *imag_arc_d, \
		     double *fft_reel_shift_norm_d, double *fft_imag_shift_norm_d)
{
  int s = CUDAINDEXE;
  const size_t tabsize = 4 * p_Nxmax * p_Nymax;

  if (s < tabsize)
    {
      const size_t c_Nxmax2 = 2 * p_Nxmax, c_Nymax2 = 2 * p_Nymax;
      const float c_rayon_sqr = p_rayon * p_rayon;
      const float c_delta_zmax_sqr = p_delta_zmax * p_delta_zmax;
      const int c_scan_size = 2 * p_Nxmax * 2 * p_Nymax;

      double arg_z_arc, z_arc; 
      int k, i, j, _x, _y;
      
      i = s / c_Nxmax2;
      j = s % c_Nxmax2;
      _y = i - p_Nxmax;
      _x = j - p_Nymax;
    
      //argument sous la racine calcul dans ARC_centre
      arg_z_arc = c_rayon_sqr - _x*_x - _y*_y;

      //ne pas depasser l'ouverture numérique pour 1 hologramme
      if(arg_z_arc > c_delta_zmax_sqr)
	{
	  //altitude au centre du volume
	  z_arc = round(sqrt(arg_z_arc) - zm0);
	  
	  //indice du tableau 1D du volume 3D
	  k = round(( - l_xm0 + _x + c_dv0s2rf) + (- l_ym0 + _y + c_dv1s2rf) * c_dv0rf + \
		    (z_arc + c_dv2s2rf) * c_dv1xdv0rf);

	  sup_redon_d[k] += 1;//pour calculer le support
	  
	  reel_arc_d[k] += fft_reel_shift_norm_d[s];//pour calculer l'image
	  imag_arc_d[k] += fft_imag_shift_norm_d[s];//pour calculer l'image
	}
    }
}





// *****************************************************************************
// centered 2D crop
// *****************************************************************************



// -----------------------------------------------------------------------------
// kernel CUDA


__global__ void
kernel_crop2D_double(double *src_data, double *dst_data, size_t src_size, size_t dst_size)
{
  int s = CUDAINDEXE;
  size_t tabsize = src_size * src_size;

  size_t xs, ys, xd, yd;
  size_t xmin, xmax;
  size_t delta;

  //ASSERT( src_size >= dst_size);

  if (s < tabsize)
    {
      delta = (src_size - dst_size) / 2;

      ys = s / src_size;
      xs = s % src_size;

      xd = xs - delta;
      yd = ys - delta;

      if (xd >= 0 && xd < dst_size && \\
	  yd >= 0 && yd < dst_size)
	dst_data[ yd * dst_size + xd ] = src_data[ ys * src_size + xs ];
    }
    
}



// *****************************************************************************
// peak_normalize
// *****************************************************************************


// 	const size_t cpt_max = peak_normalize(c_NxNy4, fft_reel_shift_h, fft_imag_shift_h, \
// 					      fft_reel_shift_norm, fft_imag_shift_norm);

// -----------------------------------------------------------------------------
// kernel CUDA



__global__ void
kernel_norm_double(double *dst_data, double *src_data_x, double *src_data_y, size_t src_size)
{
  int s = CUDAINDEXE;
  if (s < src_size)
    {
      double X = src_data_x[s];      
      double Y = src_data_y[s];
      dst_data[s] = X * X + Y * Y; 
    }
}


__global__ void
kernel_normpeak_double(double *dst_data_r, double *dst_data_i, double *src_data_r, double *src_data_i, size_t src_size, double max_part_reel, double max_part_imag, double inv_max_mod)
{
  int s = CUDAINDEXE;
  if (s < src_size)
    {
      double REAL = src_data_r[s];
      double IMAG = src_data_i[s];
      dst_data_r[s] = (REAL * max_part_reel + IMAG * max_part_imag) * inv_max_mod;
      dst_data_i[s] = (IMAG * max_part_reel - REAL * max_part_imag) * inv_max_mod;
    }
}


// data_r/i : données e/s. 
// data_size_edge: côté du carré que le tableau r/i représente
// data_mask: tableau contenant le masque à appliquer
// data_mask_edge: on suppose le masque carré, on donne uniquement la valeur du côté
// spot: on applique le masque aux images r et i centré sur spot
__global__ void
kernel_applymask_double(double* data_r, double* data_i, size_t data_size_edge, \
			unsigned char* data_mask, size_t data_mask_edge, \
			size_t spot_x, size_t spot_y)
{
  const size_t data_mask_halfedge = data_mask_edge / 2;
  const size_t x_min = spot_x - data_mask_halfedge;
  const size_t x_max = spot_x + data_mask_halfedge;
  const size_t y_min = spot_y - data_mask_halfedge;
  const size_t y_max = spot_y + data_mask_halfedge;
  const double inv_255 = 1.0f / 255.0f;
  const size_t data_size = data_size_edge * data_size_edge;

  int s = CUDAINDEXE;
 
  if (s < data_size)
    {
      size_t x = s / data_size_edge;
      size_t y = s % data_size_edge;

      // si le voxel est dans la fenêtre
      if (( x > x_min ) && ( x < x_max ) && ( y > y_min ) && ( y < y_max ))
	{
	  size_t twin_index = data_mask_edge * ( y - spot_y + data_mask_halfedge ) + \
	    x - spot_x + data_mask_halfedge;
	  double twin_val_inv255 = data_mask[ twin_index ] * inv_255;
	  
	  data_r[s] *= twin_val_inv255;
	  data_i[s] *= twin_val_inv255;
	}
    }
}
