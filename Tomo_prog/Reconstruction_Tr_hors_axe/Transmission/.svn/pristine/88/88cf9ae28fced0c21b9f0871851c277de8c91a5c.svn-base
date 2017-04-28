#ifndef _CUDA_HOLO_PROCESS_
#define _CUDA_HOLO_PROCESS_


#include <stdio.h>
#include <iostream>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>


#include <thrust/reduce.h>
#include <thrust/sequence.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>


using namespace std;

#include "Holo_Process.h"
#include "cu_Holo_Process_kernels.h"


#include "cu_Array.h"
#include "cuFFT_Image.h"


// =============================================================================
//
// =============================================================================


class cu_Holo_Process : public Holo_Process
{
  // ==================================================
  dim3 dimGrid_window;  // dimension de la fenêtre dans l'image camera
  dim3 dimBlock_window;

  dim3 dimGrid_cropped; // dimension du crop dans fourier -> fréquences instrumentales
  dim3 dimBlock_cropped;

  // ==================================================
 public:
  double *fft_reel_shift_norm_d, *fft_imag_shift_norm_d;
  double *fft_reel_shift_norm_h, *fft_imag_shift_norm_h;

  double *plan_reel_h, *plan_imag_h;
  double *plan_reel_d, *plan_imag_d;

  double *fft_reel_shift_d, *fft_imag_shift_d;
  double *fft_reel_shift_h, *fft_imag_shift_h;
  
  cuFFT_Image<double, cufftDoubleComplex> *VF_Fourier2D_d;

  // image sur laquelle on calcule la norme uniquement pour trouver l'indice du max
  double* fft_norm_peak_d;

  unsigned char* cache_jumeau_d;
  

  // ==================================================
 public:

 cu_Holo_Process() : Holo_Process()
    { ; }
  
  
  ~cu_Holo_Process() 
    { 
      cudaFree(fft_norm_peak_d);
      cudaFree(fft_reel_shift_norm_d);
      cudaFree(fft_imag_shift_norm_d);
      cudaFree(fft_reel_shift_norm_h); cudaFree(fft_imag_shift_norm_d);
      cudaFree(plan_reel_d); cudaFree(plan_imag_d); 
      cudaFree(plan_reel_h); cudaFree(plan_imag_h); 
      cudaFree(fft_reel_shift_d); cudaFree(fft_imag_shift_d); 
      cudaFree(fft_reel_shift_h); cudaFree(fft_imag_shift_h); 
    }


  // *****************************************************************************
  //
  // *****************************************************************************
  

  bool set_prepare()
  {
    bool success = Holo_Process::set_prepare();
    
    CUDALLOC_HP(plan_reel_h, c_N, double);
    CUDALLOC_HP(plan_imag_h, c_N, double);

    CUDALLOC_HP(fft_reel_shift_h, c_N, double);
    CUDALLOC_HP(fft_imag_shift_h, c_N, double);
    

    CUDALLOC_GPU(plan_reel_d, c_N, double);
    CUDALLOC_GPU(plan_imag_d, c_N, double);


    const size_t c_NxNy4 = 4 * p_Nxmax * p_Nymax;
    
    CUDALLOC_GPU(fft_reel_shift_d, c_NxNy4, double);
    CUDALLOC_GPU(fft_imag_shift_d, c_NxNy4, double);

    CUDALLOC_HP(fft_reel_shift_norm_h, c_NxNy4, double);
    CUDALLOC_HP(fft_imag_shift_norm_h, c_NxNy4, double);
    
    CUDALLOC_GPU(fft_reel_shift_norm_d, c_NxNy4, double);
    CUDALLOC_GPU(fft_imag_shift_norm_d, c_NxNy4, double);
    
    CUDALLOC_GPU(fft_norm_peak_d, c_NxNy4, double);
    

    VF_Fourier2D_d = new cuFFT_Image<double, cufftDoubleComplex>(p_window_size, p_window_size);
    VF_Fourier2D_d -> allocate();
    
    // --------------------------------------------------
    CUDALLOC_GPU( cache_jumeau_d, c_mask_size + c_mask_size, unsigned char);
    CUDAPUSH( cache_jumeau, cache_jumeau_d, c_mask_size + c_mask_size, unsigned char );


    // --------------------------------------------------
    // gestion des grilles de calcul

    dim3_resizeGrid( dimGrid_cropped, c_NxNy4, 512, (1<<15));
    dim3_set(dimBlock_cropped, 512, 1, 1);

    dim3_resizeGrid( dimGrid_window, c_N , 512, (1<<15));
    dim3_set(dimBlock_window, 512, 1, 1);
    

    return success;
  }


  // *****************************************************************************
  //
  // *****************************************************************************

 
  bool set_output_volumes_gpu(RECON_TYPE* out_fourier_R_d, RECON_TYPE* out_fourier_I_d, unsigned short int* out_supredon_d)
  {
    cerr << endl << "gros hack";
    
    i_reel_arc = out_fourier_R_d;
    i_imag_arc = out_fourier_I_d;
    i_sup_redon = out_supredon_d;
     
    return i_prefs_output_vols = true;
     
  }


  // *****************************************************************************
  //
  // *****************************************************************************

  
  // seul le lancement de kernel change
  void angle_mapping_gpu(size_t pl_current_angle, int pl_xmi, int pl_ymi, double* in_fft_r, double* in_fft_i)
  {
    //coordonnée dans l'image2D centrée (xm0,ym0)=(0,0)=au centre de l'image
    const int l_xm0 = pl_xmi - p_Nxmax;
    //printf("xm0 %i \n", xm0);
    const int l_ym0 = pl_ymi - p_Nymax;

    if(l_xm0 == 0 && l_ym0 == 0)
      printf("(xm0,ym0)=(0,0) dans le plan %d\n", pl_current_angle);
      

    if(( sqr(l_xm0) + sqr(l_ym0)) > c_rayon_inf)
      {
	i_centres_exclus++; 
	MSG_ASSERT(false, "lol, encore jamais arrivé!!");
      }
    else
      {
	//sauvegarde des centres en coordonnées non centrée; on met le numero d'angle
	i_centre[pl_xmi * 2 * p_Nxmax + pl_ymi] = pl_current_angle;
	//pour vérifier les pb. A ouvrir en sizeof(int)=32 bits mais non signé

	//indice 3D
	// variables locales uniquement
	int k = 0;
	double arg_z_arc = 0;
	double z_arc = 0;
	double zm0;
	double zm0_carre = sqr(p_rayon) - sqr(l_xm0) - sqr(l_ym0);

	if(round(zm0_carre) > -1)
	  {
	    zm0 = sqrt( zm0_carre );
	    //printf("zm0:%i\n",zm0);
	    i_nb_proj++;

	    kernel_mapping_float<<<dimGrid_window, dimBlock_window>>> (l_xm0, l_ym0, p_Nxmax, p_Nymax, \
					       p_rayon, p_delta_zmax, zm0, \
					       c_dv0s2rf, c_dv1s2rf, c_dv0rf, c_dv2s2rf, c_dv1xdv0rf, \
					       i_sup_redon, i_reel_arc, i_imag_arc, \
					       fft_reel_shift_norm_d, fft_imag_shift_norm_d);
	    

	  }//fin if zm0>-1
      }//fin else xm0_limite

  }

    
  
  // *****************************************************************************
  //
  // *****************************************************************************
  
  size_t peak_normalize_gpu2(size_t cropped_size, double* in_fft_r, double* in_fft_i, \
			     double* out_fft_r, double* out_fft_i)
  {
    double max_part_reel, max_part_imag, max_module, max_mod_divisor;
    size_t cpt_max = 0;
    double max_val = 0; // valeur maximale du module des TFs, >0
    size_t c_NxNy4 = 2 * p_Nxmax * 2 * p_Nymax;

    
    kernel_norm_double<<<dimGrid_cropped, dimBlock_cropped>>>(fft_norm_peak_d, in_fft_r, in_fft_i, cropped_size);
    thrust::device_ptr<double> ptr_d ( fft_norm_peak_d );
    thrust::device_ptr<double> ptrmax = thrust::max_element( ptr_d, ptr_d + c_NxNy4 );

    cpt_max = ptrmax - ptr_d; //cropped_size / 2; 
    ASSERT( cpt_max < cropped_size );

    // cu_Array<double> array_in_fft_r( in_fft_r, cropped_size); // dès que c'est déclaré, ça pète 
    // cu_Array<double> array_in_fft_i( in_fft_i, cropped_size);
    
    //max_part_reel = array_in_fft_r.get(cpt_max);
    CUDAPULL( &max_part_reel, in_fft_r + cpt_max, 1, double);
    //max_part_imag = 0; // array_in_fft_i.get(cpt_max);
    CUDAPULL( &max_part_imag, in_fft_i + cpt_max, 1, double);
    max_module = sqr(max_part_reel) + sqr(max_part_imag);
    max_mod_divisor = 1.0f / max_module;

    kernel_normpeak_double<<<dimGrid_cropped, dimBlock_cropped>>>(out_fft_r, out_fft_i, in_fft_r, in_fft_i, cropped_size, max_part_reel, max_part_imag, max_mod_divisor);

    return cpt_max;
  }


  // *****************************************************************************
  //
  // *****************************************************************************


  void twin_eliminate_gpu(int pl_xmi, int pl_ymi, double* io_fft_r, double* io_fft_i, unsigned char* cache_jumeau_d)
  {
    const int l_xc = pl_xmi - p_Nxmax;
    const int l_yc = pl_ymi - p_Nymax;

  
      
    ////////////////////////////////////Virer l'objet jumeaux (0.001 à 5 s cumulé sur CPU)
    if(sqr(l_xc) + sqr(l_yc) > 900)//35*35=1225;objet jumeau pas trop près de l'objet
      { 
	i_jumeau_elimine++;

	const size_t c_Nxmax2 = 2 * p_Nxmax;
	const size_t c_Nymax2 = 2 * p_Nymax;
	const size_t c_NxNy4 = c_Nxmax2 * c_Nymax2;
		

	//coordonnée objet jumeaux
	const int l_xmij = c_Nxmax2 - pl_xmi;
	const int l_ymij = c_Nymax2 - pl_ymi;


	kernel_applymask_double<<<dimGrid_cropped, dimBlock_cropped>>> \ 
	  (io_fft_r, io_fft_i, c_Nymax2, cache_jumeau_d, c_mask_size, l_xmij, l_ymij);

		  
      }//fin if objet jumeau pas trop pres

  }
    
  
  // *****************************************************************************
  //
  // *****************************************************************************


  void batch_launch()
  {
    ASSERT(i_batch_ready);


    // constantes locales

    size_t batch_size = v_batch_angles.size();

    const size_t c_window_dim_x = p_window_size;
    const size_t c_window_dim_y = p_window_size;
    const size_t c_image_size = c_window_dim_x * c_window_dim_y;

    const size_t c_Nxmax2 = 2 * p_Nxmax, c_Nymax2 = 2 * p_Nymax;
    const size_t c_NxNy4 = 4 * p_Nxmax * p_Nymax;
    const size_t c_mask_halfsize = c_mask_size / 2;
    const double c_inv_255 = 1.0f / 255.0f;
    const float c_rayon_sqr = sqr(p_rayon);
    const float c_delta_zmax_sqr = sqr(p_delta_zmax);
    
    const size_t c_cropped_fourier_size = c_NxNy4;


    // peut-être à virer
    size_t i_grid_GridX = 2 * p_Nxmax;
    size_t i_grid_GridY = 2 * p_Nymax;
    size_t i_grid_ThreadsPerBlock = 512;
    ASSERT(p_window_size >= 2 * p_Nxmax);


    for (size_t i = 0; i < batch_size; i++) 
      {
	size_t current_angle = v_batch_angles[i];

	// affichage périodique progression
	if (! (current_angle % 50))
	  {
	    cout << endl << "processing: " << current_angle;
	  }
      
      

	//----------------------------------------------------------------------
	//////////////////////////phase shifting////////////////////////////////
	//----------------------------------------------------------------------

	phase_shifting(current_angle, plan_reel_h, plan_imag_h);
	//show_hdp();

	//----------------------------------------------------------------------
	/////////////// TF et crop dans les fréquences instrumentales
	//----------------------------------------------------------------------

	CUDAPUSH(plan_reel_h, plan_reel_d, c_image_size, double );
	CUDAPUSH(plan_imag_h, plan_imag_d, c_image_size, double );

	VF_Fourier2D_d -> import_from(plan_reel_d, plan_imag_d);
	VF_Fourier2D_d -> set_fourier_forward();
	VF_Fourier2D_d -> export_to(plan_reel_d, plan_imag_d); // re-circshifté
	
	
	kernel_crop2D_double<<<dimGrid_window, dimBlock_window>>>(plan_reel_d, fft_reel_shift_d, p_window_size, 2 * p_Nxmax );
	kernel_crop2D_double<<<dimGrid_window, dimBlock_window>>>(plan_imag_d, fft_imag_shift_d, p_window_size, 2 * p_Nxmax );
	CUDASYNCE;  

	// ----------------------------------------------------------------------
	// Recherche du maximum non centré  +
	// Normalisation par le pic central
	// ----------------------------------------------------------------------
	
	const size_t cpt_max = peak_normalize_gpu2(c_NxNy4, fft_reel_shift_d, fft_imag_shift_d, \
						   fft_reel_shift_norm_d, fft_imag_shift_norm_d);
	

	// ----------------------------------------------------------------------

	//Coordonnées dans l'espace 2D à partir de l'indice 1D: (xc,yc)=(0,0)=en haut à gauche de l'image
	const int l_xmi = cpt_max % c_Nxmax2;
	const int l_ymi = cpt_max / c_Nymax2;


	twin_eliminate_gpu(l_xmi, l_ymi, fft_reel_shift_norm_d, fft_imag_shift_norm_d, cache_jumeau_d);

	angle_mapping_gpu(current_angle, l_xmi, l_ymi, fft_reel_shift_norm_d, fft_imag_shift_norm_d);


	// ----------------------------------------------------------------------
	// if (show_hdp_p) cvWaitKey(1);
      
	i_angle_last_done = current_angle;
      }
    // FIN ITERATION  
    

    v_batch_angles.clear();
    ASSERT( ! v_batch_angles.size() );


    i_batch_ready = false;

  }
  

};


// =============================================================================
//
// =============================================================================


#endif //_CUDA_HOLO_PROCESS_

