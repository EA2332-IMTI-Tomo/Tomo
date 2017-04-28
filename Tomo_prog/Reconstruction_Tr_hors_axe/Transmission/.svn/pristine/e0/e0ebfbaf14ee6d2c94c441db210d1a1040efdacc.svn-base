#ifndef __COMPUTE
#define __COMPUTE


#define GENERATED_FILENAME_LEN 20
#define CACHE_JUMEAU_LEN 100

#include "recon_includes.h"
#include "AIR_Volume.h"

/* extern size_t g_window_dim_x; */
/* extern size_t g_window_dim_y; */
/* extern char *g_MASKPATH; */

// for debug only
/* extern size_t image_dim_x; */
/* extern size_t image_dim_y; */

/* extern std::map<string, bool> g_map_bools; */


/*
  void
  compute_cpu_offline(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
  size_t window_edge_x, size_t window_edge_y, size_t image_dim_x, size_t image_dim_y, \
  size_t xm0_limite, size_t ym0_limite, float rayon, float delta_zmax, \
  size_t angle_start, size_t angle_count, size_t angle_jump, \
  size_t offaxis_cx, size_t offaxis_cy,		\
  const string &images_radix, const char* modref_file, std::map<string, bool> my_map_bools);
*/



// lit des interférogrammes d'un coup entre angle_start et angle_last.
// soit sur disque (renseigner images_radix) ou sur mémoire (renseigner images_list), selon la valeur de read_from_mem

void
compute_cpu_offline(std::map<string, float> my_map_experiment,	\
		    const string &images_radix, unsigned char **images_list, \
		    std::map<string, bool> my_map_bools, \
		    std::map<string, size_t> my_map_dims, \
		    std::map<string, string> my_map_paths);


void
compute_cpu_offline_batch(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
			  float rayon, float delta_zmax,			\
			  const string &images_radix, const char* modref_file, \
			  std::map<string, bool> my_map_bools, \
			  std::map<string, size_t> my_map_dims, \
			  std::map<string, string> my_map_paths);


void
compute_cpu_online_batch(std::map<string, float> my_map_experiment,	\
			 const string &images_radix, unsigned char **images_list, \
			 std::map<string, bool> my_map_bools,  \
			 std::map<string, size_t> my_map_dims,	\
			 std::map<string, string> my_map_paths);



void
compute_cpu_offline_movie(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
			  float rayon, float delta_zmax,		\
			  SlicePlanType cut, size_t slice_saved, size_t slice_gap, \
			  const string &images_radix, \
			  std::map<string, bool> my_map_bools, \ 
			  std::map<string, size_t> my_map_dims, \
			  std::map<string, string> my_map_paths);


/* deprecated */
void
compute_cpu_zoom_offline(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
			 size_t window_edge_x, size_t window_edge_y, size_t image_dim_x, size_t image_dim_y, \
			 size_t xm0_limite, size_t ym0_limite, float rayon, float delta_zmax, \
			 size_t  angle_start, size_t angle_count, size_t angle_jump, \
			 bool off_axis_and_psh, size_t offaxis_cx, size_t offaxis_cy, \
			 bool autofocus_gradient,			\
			 const string &images_radix, bool proper_numbers, bool modref_asked, const char* modref_file);

/* test only */

/*
void
compute_memtest(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
		float rayon, float delta_zmax, 
		const string &images_radix, const char* modref_file,	\
		std::map<string, bool> my_map_bools, \ 
		std::map<string, size_t> my_map_dims, \
		std::map<string, string> my_map_paths);
*/
void
compute_memtest(std::map<string, float> my_map_experiment,	\
		const string &images_radix, 
		std::map<string, bool> my_map_bools,	\
		std::map<string, size_t> my_map_dims,	\
		std::map<string, string> my_map_paths);

#endif //__COMPUTE
