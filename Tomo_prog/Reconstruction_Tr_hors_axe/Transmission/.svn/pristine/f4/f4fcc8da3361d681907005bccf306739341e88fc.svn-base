// Ce fichier sert de patron pour un appel externe du programme

#include "cpu_launcher.h"

using namespace std;


// ============================================================================
// ============================================================================



void 
compute_cpu(const char* config_filename, unsigned char** interferograms)
{
  std::map<string, bool> recon_map_bools;
  std::map<string, size_t> recon_map_dims;
  std::map<string, string> recon_map_paths;
  std::map<string, float> recon_map_physics;
  // données expérimentales calculées
  std::map<string, float> recon_map_experiment;

  enum ReconType ReconMode = _offline; 

  init_vals(recon_map_dims, recon_map_bools, recon_map_paths, recon_map_physics);
  parse_file(config_filename, recon_map_dims, recon_map_bools, recon_map_paths, recon_map_physics);
  compute_vals(recon_map_dims, recon_map_physics, recon_map_experiment);

  
  recon_map_bools["read_from_mem"] = true;
  recon_map_paths["OUTPUT_DIR"] = "/tmp";

  compute_cpu_offline(recon_map_experiment, "/dev/null", interferograms,
		      recon_map_bools, recon_map_dims, recon_map_paths);
}
