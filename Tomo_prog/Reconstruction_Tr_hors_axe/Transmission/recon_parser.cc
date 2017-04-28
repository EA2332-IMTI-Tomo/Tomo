#include "macros.h"
#include "recon_parser.h"

#include "util.h"

#include <math.h>


#include <boost/filesystem.hpp>
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/progress.hpp"


using namespace std;

#define LINESIZE 2048

// ======================================================================
//
// ======================================================================


void
init_vals(std::map<string, size_t> &map_dims, \ 
	  std::map<string, bool> &map_bools, \
	  std::map<string, string> &map_paths, \
	  std::map<string, float> &map_physics)
{
	
  // ======================================================================
  // réglages interférogrammes

  // taille de l'image caméra, devinée si lecture fichiers, passée au programme si lecture mémoire
  map_dims["image_dim_x"] = 0;
  map_dims["image_dim_y"] = 0;
  // taille de la fenêtre de découpe ( -> résolution de l'image reconstruite)
  map_dims["window_dim_x"] = 0;
  map_dims["window_dim_y"] = 0;
  // position du coin de la fenêtre de découpe dans l'image caméra
  map_dims["window_edge_x"] = 0;
  map_dims["window_edge_y"] = 0;
  map_bools["window_autocenter"] = false;
  
  // limite de lecture des fichiers d'entrée
  map_dims["angle_start"] = 0; 
  map_dims["angle_end"] = 0; 
  map_dims["angle_jump"] = 0; 
  // nombre max d'erreurs de lecture fichiers (à mettre à beaucoup si on efface des holos à la main)
  map_dims["angle_error_limit"] = 100;

  // ======================================================================
  // grandeurs physiques

  map_physics["oil_index"] = 0; //indice de l'huile // n1
  map_physics["aperture"] = 0; //ouverture numerique de l'objectif? //NA
  // (celle du condenseur intervient sur la forme, la taille, du papillon)
  map_physics["lambda"] = 0; //632 * pow(10, -9); //longueur d'onde
  map_physics["factor"] = 0; //G; // grossissement telan+objectif
  map_physics["pixel_size"] = 0; //G; // grossissement telan+objectif
  map_physics["focal_ratio"] = 0; // 1.2; 
  // float p_FOCAL_RATIO=0.8;//1.7;//0.77;//1.3;// bonne valeur le 29/9; //rapport f1/f2



  // ======================================================================
  
  // valeur controlant l'exclusion de certains centres
  map_dims["xm0_limite"] = 0; //centre maximum
  map_dims["ym0_limite"] = 0; 


  // ======================================================================
  // gestion décalage de phase et/ou hors-axe

  // mode décalage de phase, mode hors-axe combiné avec le décalage, mode hors-axe tout court
  map_bools["off_axis_mode"] = false;
  map_bools["psh_mode"] = true;
  // centre du cercle contenant les fréquences objet en hors-axe
  map_dims["off_axis_circle_cx"] = 0;
  map_dims["off_axis_circle_cy"] = 0;  
  // en hors-axe, le rayon de découpe se mesure par show_fourier. Nxmax ne convient que rarement
  map_dims["off_axis_radius"] = 0;


  // ======================================================================

  // afficher un slicer basique pendant l'écriture des résultats
  map_bools["end_visu_mode"] = true;

  // correction plan central par image de référence
  map_bools["modref_correction_mode"] = false;// modref_asked_mode = false;
  //const char* modref_file;


  // ======================================================================
  // divers

  // by default, we read 4 phase images to infer real and imaginary images
  // if true, we just read these images in a binary volume: real1img1real2img2...
  map_bools["read_saved_volume"] = false;
  map_bools["autofocus_mode_debug"] = false;
  map_bools["autofocus_mode_show_vals"] = map_bools["autofocus_mode_show_pics"] = false;
  


  // mode autofocus
  map_bools["autofocus_mode"] = false;

  // numérotation correcte hologrammes %04d %03d
  map_bools["proper_numbers_mode"] = false;
  
  // nombre de threads pour calculer la TF sur CPU
  map_dims["fftw_threads"] = 2;

  map_bools["read_from_mem"] = false;
    
  map_bools["debug_record_centers"] = true;
  map_bools["debug_record_fourier_spots"] = false;
  map_bools["debug_record_fourier_windows"] = false;
  map_bools["debug_record_final_holograms"] = false;
  map_bools["nodebug"] = false; // shunter tous les debug pour callgrind

  // essai de normalisation post FFTW
  map_bools["fourier_normalize"] = false;

  // ======================================================================
  // chemins

  map_paths["MASKPATH"] = "/dev/null";	
  map_paths["INPUT_DIR"] = "/dev/null";	
  map_paths["OUTPUT_RADIX"] = "/dev/null";	
  map_paths["INFILE_RADIX"] = "/dev/null";	
  map_paths["OUTPUT_RADIX"] = "VolumeSpecimen";
  map_paths["OUTPUT_DIR"] = "./";
  map_paths["LOGFILE"] = "recon_log.txt";


  // ==========================================================
  // ajouts sauvages
  // ==========================================================
  map_bools["compute_module"] = false;
  map_bools["save_module"] = false;
  map_bools["save_sup_redon"] = false;
  map_bools["blind_mode"] = false; // logs limités, pas d'interactivité.


  map_paths["save_sup_redon_filename"] = "supRedon_Diatomee";

}


// ======================================================================
//
// ======================================================================



void
parse_file(const char* in_file, \ 
	   std::map<string, size_t> &map_dims, \ 
	   std::map<string, bool> &map_bools, \
	   std::map<string, string> &map_paths,
	   std::map<string, float> &map_physics)
{
  ifstream cfg_file(in_file);

  ASSERT(cfg_file.is_open());

  char read_line[LINESIZE];
  char *read_token, *variable_name, *value_name;

  while (! cfg_file.eof() )
    {	
      cfg_file.getline(read_line, LINESIZE, '\n');
      read_token = strtok(read_line, " ");
      if ((read_token == NULL) || (read_token[0] == '#')) 
	continue;

      variable_name = read_token;
      read_token = strtok(NULL, " ");
      
      if (read_token == NULL) {
	cerr << endl << "parse error: value expected after " << variable_name << " ; continue reading";
	continue;
      }
      else
	value_name = read_token;
      
      // strtok est destructif: si le token lu est une chaîne, il doit être recopié

      SWITCH_STR_BEGIN( variable_name )
	CASE_STR("WINDOW_DIMX") map_dims["window_dim_x"] = atoi(value_name);
      CASE_STR("CCD_DIMX") map_dims["window_dim_x"] = atoi(value_name);  // legacy
      CASE_STR("WINDOW_DIMY") map_dims["window_dim_y"] = atoi(value_name);
      CASE_STR("CCD_DIMY") map_dims["window_dim_y"] = atoi(value_name); // legacy
      CASE_STR("XM0_LIMIT") map_dims["xm0_limite"] = atoi(value_name);
      CASE_STR("YM0_LIMIT") map_dims["ym0_limite"] = atoi(value_name);
      CASE_STR("EDGE_X") map_dims["window_edge_x"] = atoi(value_name);
      CASE_STR("EDGE_Y") map_dims["window_edge_y"] = atoi(value_name);
      CASE_STR("EDGES_CENTER") map_bools["window_autocenter"] = (bool)atoi(value_name);
      CASE_STR("EDGE_CENTER") map_bools["window_autocenter"] = (bool)atoi(value_name); // car confusions
      CASE_STR("IMAGE_DIMX") map_dims["image_dim_x"] = atoi(value_name); // uniquement si lecture interférogrammes mémoire
      CASE_STR("IMAGE_DIMY") map_dims["image_dim_y"] = atoi(value_name);
      CASE_STR("FIRST_ANGLE") map_dims["angle_start"] = atoi(value_name);
      CASE_STR("FIRST_LANE") map_dims["angle_start"] = atoi(value_name); // legacy
      CASE_STR("FINAL_ANGLE") map_dims["angle_end"] = atoi(value_name);
      CASE_STR("ANGLE_ERROR_MAX") map_dims["angle_error_limit"] = atoi(value_name);
      CASE_STR("INC_ANGLE") map_dims["angle_jump"] = atoi(value_name);
      CASE_STR("OIL_INDEX") map_physics["oil_index"] = atof(value_name);
      CASE_STR("APERTURE") map_physics["aperture"] = atof(value_name);
      CASE_STR("LAMBDA") map_physics["lambda"] = atof(value_name);
      CASE_STR("FACTOR") map_physics["factor"] = atof(value_name);
      CASE_STR("PIXEL_SIZE") map_physics["pixel_size"] = atof(value_name);
      CASE_STR("RF_VAL") map_physics["focal_ratio"] = atof(value_name);
      CASE_STR("FFTW_THREADS") map_dims["fftw_threads"] = atoi(value_name);
#ifdef CUDA
      //CASE_STR("FFT3D_CUDA") g_fft3d_cuda = (bool)atoi(value_name);
#endif
      CASE_STR("READ_DUMPED_VOLUME") map_bools["read_saved_volume"] = (bool)atoi(value_name); 
      CASE_STR("MASKPATH") map_paths["MASKPATH"] = str_alloc_cpy(value_name);
      CASE_STR("OUTPUT_DIR") map_paths["OUTPUT_DIR"] = str_alloc_cpy(value_name); 
      CASE_STR("OUTPUT_RADIX") map_paths["OUTPUT_RADIX"] = str_alloc_cpy(value_name);
      CASE_STR("INPUT_DIR") map_paths["INPUT_DIR"] = str_alloc_cpy(value_name);
      CASE_STR("INFILE_RADIX") map_paths["INFILE_RADIX"] = str_alloc_cpy(value_name);

      CASE_STR("PHASE_SHIFT") map_bools["psh_mode"] = atoi(value_name);
      CASE_STR("OFF_AXIS") map_bools["off_axis_mode"] = atoi(value_name);
      CASE_STR("CIRCLE_CX") map_dims["off_axis_circle_cx"] = atoi(value_name);
      CASE_STR("CIRCLE_CY") map_dims["off_axis_circle_cy"] = atoi(value_name);
      CASE_STR("CIRCLE_R") map_dims["off_axis_radius"] = atoi(value_name);

      CASE_STR("AUTOFOCUS") map_bools["autofocus_mode"] = (bool)atoi(value_name);
      CASE_STR("AUTOFOCUS_DEBUG") map_bools["autofocus_mode_debug"] = (bool)atoi(value_name);
      CASE_STR("AUTOFOCUS_VALS") map_bools["autofocus_mode_show_vals"] = (bool)atoi(value_name);
      CASE_STR("AUTOFOCUS_PICS") map_bools["autofocus_mode_show_pics"] = (bool)atoi(value_name);

      CASE_STR("PROPER_NUMBERS") map_bools["proper_numbers_mode"] = (bool)atoi(value_name);
      CASE_STR("REF_CORRECTION") map_bools["modref_correction_mode"] = (bool)atoi(value_name);
      CASE_STR("REF_FILE") map_paths["modref_file"] = str_alloc_cpy(value_name);

      CASE_STR("FOURIER_NORMALIZE") map_bools["fourier_normalize"] = (bool)atoi(value_name); // +dummy val

      CASE_STR("NO_END_VISU") map_bools["end_visu_mode"] = (bool)atoi(value_name); //false; //+dummy val


      // file saved on output dir
      CASE_STR("SUPREDON_SAVENAME") { map_paths["save_sup_redon_filename"] = str_alloc_cpy(value_name); map_bools["save_sup_redon"] = true;}

      CASE_STR("FOURIER_SAVENAME") { map_paths["save_fourier_filename"] = str_alloc_cpy(value_name); map_bools["save_fourier"] = true;}
      // sauvegarde hologrammes recalés, gourmand. 
      CASE_STR("RECORD_HOLOGRAMS") { map_bools["debug_record_final_holograms"] = (bool)atoi(value_name);} 
      
      CASE_STR("ALL_DEBUG_ON") // + dummy argt
	{ 
	  map_bools["autofocus_mode_debug"] = map_bools["debug_record_fourier_spots"] = map_bools["debug_record_fourier_windows"] = map_bools["debug_record_final_holograms"] = map_bools["save_sup_redon"] = map_bools["save_fourier"] = true;
	  map_paths["save_sup_redon_filename"] = "debug_sup_redon";
	  map_paths["save_fourier"] = "debug_fourier";
	}
      CASE_STR("ALL_DEBUG_OFF") { // pas d'arguments!
	map_bools["autofocus_mode_debug"] = map_bools["debug_record_fourier_spots"] = map_bools["debug_record_fourier_windows"] = map_bools["debug_record_final_holograms"] = map_bools["save_sup_redon"] = map_bools["save_fourier"] =  false;
      }

      

      SWITCH_STR_ELSE_END( 
			  cerr << endl << "parse error skipped: "		\
			  << variable_name << " | " << value_name;  
			   )
	}
  cfg_file.close();
  
}


// ======================================================================
//
// ======================================================================


void
parse_argts(size_t l_argc, char** l_argv, \ 
	    std::map<string, size_t> &map_dims, \ 
	    std::map<string, bool> &map_bools, \
	    std::map<string, string> &map_paths,
	    std::map<string, float> &map_physics,\
	    movieParams &movie_settings,\
	    enum ReconType &ReconMode)
{
  char read_line[LINESIZE];
  char *read_token, *variable_name, *value_name;

  while (l_argc > 0)
    {
      variable_name = l_argv[0];
      if (!(l_argc - 1)) break;
      value_name = l_argv[1];

      SWITCH_STR_BEGIN( variable_name )
	CASE_STR("-window_dimx") map_dims["window_dim_x"] = atoi(value_name);
      CASE_STR("-window_dimy") map_dims["g_window_dim_y"] = atoi(value_name);
      CASE_STR("-xm0_limit") map_dims["xm0_limite"] = atoi(value_name);
      CASE_STR("-ym0_limit") map_dims["ym0_limite"] = atoi(value_name);
      CASE_STR("-edge_x") map_dims["window_edge_x"] = atoi(value_name);
      CASE_STR("-edge_y") map_dims["window_edge_y"] = atoi(value_name);
      CASE_STR("-edges_center") map_bools["window_autocenter"] = (bool)atoi(value_name);
      CASE_STR("-image_dim_x") map_dims["image_dim_x"] = atoi(value_name); // uniquement si lecture interférogrammes mémoire
      CASE_STR("-image_dim_y") map_dims["image_dim_y"] = atoi(value_name);
      CASE_STR("-first_angle") map_dims["angle_start"] = atoi(value_name);
      CASE_STR("-final_angle") map_dims["angle_end"] = atoi(value_name);
      CASE_STR("-inc_angle") map_dims["angle_jump"] = atoi(value_name);
      CASE_STR("-err_angle") map_dims["angle_error_limit"] = atoi(value_name);
      CASE_STR("-oil_index") map_physics["oil_index"] = atof(value_name);
      CASE_STR("-aperture") map_physics["aperture"] = atof(value_name);
      CASE_STR("-lambda") map_physics["lambda"] = atof(value_name);
      CASE_STR("-factor") map_physics["factor"] = atof(value_name);
      CASE_STR("-pixel_size") map_physics["pixel_size"] = atof(value_name);
      CASE_STR("-rf_val") map_physics["focal_ratio"] = atof(value_name);
      CASE_STR("-fftw_threads") map_dims["fftw_threads"] = atoi(value_name);
#ifdef CUDA
      // CASE_STR("-cuda_3d") g_fft3d_cuda = (bool)atoi(value_name);
#endif
      CASE_STR("-masks") map_paths["MASKPATH"] = value_name;
      CASE_STR("-input_dir") map_paths["INPUT_DIR"] = value_name;
      CASE_STR("-input") map_paths["INPUT_DIR"] = value_name; // legacy
      CASE_STR("-output_dir") map_paths["OUTPUT_DIR"] = value_name;
      CASE_STR("-output_name") map_paths["OUTPUT_RADIX"] = value_name;
      CASE_STR("-in_radix") map_paths["INFILE_RADIX"] = value_name;
      CASE_STR("-realtime_mode") ReconMode = _online_batch; //  
      CASE_STR("-psh") map_bools["psh_mode"] = (bool)atoi(value_name);
      CASE_STR("-off_axis") map_bools["off_axis_mode"] = (bool)atoi(value_name);
      CASE_STR("-circle_cx") map_dims["off_axis_circle_cx"] = atoi(value_name);
      CASE_STR("-circle_cy") map_dims["off_axis_circle_cy"] = atoi(value_name);
      CASE_STR("-circle_r") map_dims["off_axis_radius"] = atoi(value_name);
      CASE_STR("-autofocus") map_bools["autofocus_mode"] = true; // autofocus par gradient
      CASE_STR("-proper_numbers_mode") map_bools["proper_numbers_mode"] = true; // les numéros d'hologrammes normalisés en nb chiffres
      CASE_STR("-movie_on_slice") {ReconMode = _offline_movie; movie_settings.movie_on_slice = atoi(value_name);}
      CASE_STR("-movie_every_angle") {ReconMode = _offline_movie; movie_settings.movie_angle_gap = atoi(value_name);}
      CASE_STR("-movie_cut_plane") {ReconMode = _offline_movie; movie_settings.movie_cut_view = (SlicePlanType) atoi(value_name);}
      CASE_STR("-nv") map_bools["end_visu_mode"] = false;
      CASE_STR("-fourier_normalize") map_bools["fourier_normalize"] = true;
      CASE_STR("-mod_ref") { map_paths["modref_file"] = value_name; map_bools["modref_correction_mode"] = true; }

      SWITCH_STR_ELSE_END(
			  cerr << "critical parse error met" << endl; 
			  //usage(argc, argv);
			  fflush(stdout); cout.flush();
			  exit (EXIT_FAILURE);
			  );

      l_argv += 2; l_argc -= 2;
    }

}


// ======================================================================
// capter dès le lancement toutes les erreurs détectables.
// ======================================================================


void
check_vals(std::map<string, size_t> &map_dims, std::map<string, bool> &map_bools, std::map<string, string> &map_paths, std::map<string, float> &map_experiment)
{
  // =============================================================================
  // A) vérifier que la fenêtre de découpe reste dans l'image // FINI
  
  /*
  // NB! il peut encore y avoir autocentrage map_bools["window_autocenter"]

  if((map_dims["window_edge_x"] + map_dims["window_dim_x"] > map_dims["image_dim_x"]) || \
     (map_dims["window_edge_y"] + map_dims["window_dim_y"] > map_dims["image_dim_y"]))
    {
      cerr << endl << map_dims["window_edge_x"] << "+" << map_dims["window_dim_x"] << " > " <<  map_dims["image_dim_x"] << " OR"; 
      cerr << map_dims["window_edge_y"] << "+" << map_dims["window_dim_y"] << " > " <<  map_dims["image_dim_y"];
      MSG_ASSERT(false, "recon_parser: erreur critique: la zone decoupee sort de l'image\n");
    }
  // largeur : 740, hauteur : 574 // 1312 * 1082
  */
  
  // NON PLUS MAINTENANT. Transféré dans Holo_Process, pas besoin de doublons, surtout maintenant
  // que la taille de l'image est déduite du fichier image

  
  // =============================================================================
  // B) vérifier que le répertoire de sortie est inscriptible
  // c'est plus simple que pour l'interface: le répertoire DOIT être présent et vide et inscriptible
  // c'est à TomOsiris de fournir un tel répertoire, et on s'en fiche de ce qu'il y a dessus

  // doublons de code avec vectradialog_util.cc:probeImagesOutDir

  namespace bf = boost::filesystem;
  bf::path outDir( map_paths["OUTPUT_DIR"] );
  bf::path parentPath( outDir.parent_path() );


  if (! (bf::exists( outDir ) && bf::is_directory( outDir )))
    {
      cerr << endl << "recon_parser: invalid output dir: " << outDir; 
      MSG_ASSERT(false, "invalid outDir");
    }

  try
    {
      bf::path testpath( parentPath / "toto_dir");
      testpath += "_rr"; // insérer ici chaine aléatoire
      bool create_dir = bf::create_directory( testpath );    
      bool delete_dir = bf::remove( testpath );  
    
      if (! create_dir) 
	cerr << endl << "recon_parser: Can't write in given directory" << testpath;
      else 
	cout << endl << "test écriture dans répertoire parent ok avec:" << testpath;		
      
      MSG_ASSERT( create_dir && delete_dir, "critical failure when testing outdir writing");

    }
  catch(bf::filesystem_error e) 
    {
      cerr << e.what(); 
      MSG_ASSERT( false, "critical failure when testing outdir writing");
    }


  // =============================================================================
  // B) vérifier que les masques sont disponibles // BIENTOT FINI

  bf::path maskDir( map_paths["MASKPATH"] );
  MSG_ASSERT( (bf::exists( maskDir ) && bf::is_directory( maskDir )), \ 
	      "input directory (masks) incorrect");
  
  //string maskName("kmasq_"); maskName += "x"; maskName += ".pgm";

  /*  // Fini, on le génère maintenant
  ARRAY_DEC_ALLOC(maskName, 512, char);
  sprintf(maskName, "kmasq_%dx%d_30.pgm", map_dims["window_dim_x"], map_dims["window_dim_y"] );
  bf::path maskPath( maskDir / maskName );
  
  MSG_ASSERT( bf::is_regular_file(maskPath), "Mask file not available");
  */

  // =============================================================================
  // B) vérifier que les données sont disponibles
  
  bf::path inDir( map_paths["INPUT_DIR"] );

  MSG_ASSERT( (bf::exists( inDir ) && bf::is_directory( inDir )), \ 
	      "input directory (interferograms) incorrect");
  
  // penser à lister le répertoire avec .pgm et input_radix, vérifier fichiers lisibles


  // =============================================================================
  // B) calculer la place nécessaire sur disque

  size_t angles_count = (map_dims["angle_end"] - map_dims["angle_start"]) / map_dims["angle_jump"];
  MSG_ASSERT( angles_count, "recon_parser: aucun angle ne sera balayé");

  size_t cube_edge = (size_t) map_experiment["cube_edge"];
  // 512^3 *2 (complex) *4 (32 bits)
  size_t reconstruction_volume_size = cube_edge * cube_edge * cube_edge * 2 * 4;
  size_t Nxmax = (size_t) map_experiment["nxmax"];
  size_t Nymax = (size_t) map_experiment["nxmax"]; // "volontaire"
  size_t debug_slice_size = 4 * Nxmax * Nymax * 8; // 64 bits = 4o
  size_t debug_slice_stack = debug_slice_size * angles_count;

  MSG_ASSERT( debug_slice_size, "recon_parser: erreur avec Nxmax ou Nymax");

  size_t stacks_count = 0;
  stacks_count += (size_t) map_bools["debug_record_fourier_spots"]; 
  stacks_count += (size_t) map_bools["debug_record_fourier_windows"]; 
  stacks_count += (size_t) map_bools["debug_record_final_holograms"];  // pb: on ne le sait qu'a posteriori dans mon prog
  size_t mandatory_size = reconstruction_volume_size + stacks_count * debug_slice_stack;

  if (map_bools["save_sup_redon"]) mandatory_size += cube_edge * cube_edge * cube_edge * 1;


  // et vérifier qu'il y a assez de place (nb: on omet de vérifier si on réécrit)
  size_t space_avail_b;

  try{
    bf::space_info infoDisk;
    infoDisk = bf::space(parentPath);
    space_avail_b = size_t(infoDisk.available);
  }
  catch(bf::filesystem_error e) {
    cerr << endl << "recon_parser: failed to check available space on: <" << outDir << ">";
  }

  MSG_ASSERT( space_avail_b > mandatory_size, "recon_parser: not enough space for reconstruction and debug volumes"); 
  float giga_b = 1024 * 1024 * 1024;
  cout << endl << float(mandatory_size) / float(giga_b) << " Gb will be written (" << float( space_avail_b ) / float(giga_b) << " available)";
  

  // =============================================================================
  // fini
  
  cout << endl << "basic tests passed";
  cout.flush(); cerr.flush();
}


// ======================================================================
//
// ======================================================================


void
compute_vals(std::map<string, size_t> &map_dims, \ 
	     std::map<string, float> &map_physics,	\
	     std::map<string, float> &map_experiment)
{
  float p_LAMBDA = map_physics["lambda"];
  float p_FACTORg = map_physics["factor"];
  float p_APERTURE = map_physics["aperture"];
  float p_PIXEL_SIZE = map_physics["pixel_size"];
  float p_FOCAL_RATIO = map_physics["focal_ratio"];  
  float p_OIL_INDEX = map_physics["oil_index"];  
  

  //short int precision;
  float precision;

  size_t angle_count = map_dims["angle_end"] - map_dims["angle_start"] + 1;

  float theta = asin(map_physics["aperture"]/map_physics["oil_index"]);//theta_max
  printf("theta : %f\n",theta);

  // taille d'un pixel de l'image reconstruite, probablement à un facteur genre 512/2xNxMax près
  float pixel_size_recon = p_FOCAL_RATIO * p_PIXEL_SIZE /  p_FACTORg ;
  cout << endl << "taille d'un pixel (non zoomé) sur l'image reconstruite (nm):" << pixel_size_recon * 1000000000 << endl;
  map_experiment["pixel_size_recon"] = pixel_size_recon;

  /*
    # grossissement telan+objectif 
    FACTOR 100
    # Rapport focale 0.83333 avant, 1 pour experience fom 2013, 0.5 pour ictms 13
    RF_VAL 0.5
    
    # Taille des pixels: (sur la CCD camera)
    #11.2e-6 ancienne cam, 8e-6 now 
    PIXEL_SIZE 0.0000080
    
    #taille d'un pixel sur l'image reconstruite:
    # au 100x:
    #8 microns / 100 = 80nm | avec le RF à 0.5 => 40nm. Oui, c'est petit, mais c'est pas la résolution. 
  */

  //   char *images_radix;
  //   ARRAY_ALLOC(images_radix, strlen(g_INPUT_DIR) + strlen(g_INFILE_RADIX) + 10, char);
  
  // Facteur d'echelle
  float K_SCALE = p_LAMBDA * p_FACTORg/(2 * p_APERTURE * p_PIXEL_SIZE * p_FOCAL_RATIO); //idem HA
  printf("K_SCALE : %f\n", K_SCALE);
  float PIXEL_SIZE2 = p_PIXEL_SIZE;//*K_SCALE; //Taille pixel apres mise à l'echelle important si RF<>1

  int Ntx = map_dims["window_dim_x"];//K_SCALE; //Nombre total de pixel en x
  /*int Nty=Ntx;*/ //(support carré)

  // rayon = le rayon de la sphère dans laquelle on découpe une calotte != rayon calotte
  // le rayon calotte correspond à Nxmax, et en HA pur, 
  // il faut le (Nxmax) remplacer par le rayon mesuré sur show_fourier

  // idem en code HA mat, 
  float rayon = p_OIL_INDEX * p_FOCAL_RATIO * PIXEL_SIZE2 * map_dims["window_dim_x"]/(p_FACTORg * p_LAMBDA); //Rayon
  //int rayon=round(p_OIL_INDEX*p_FOCAL_RATIO*p_PIXEL_SIZE*Ntx/(G*p_LAMBDA));ici un commentaire archéologique

  float delta_zmax = rayon * cos(theta);//float?

  //int rayon=100;

  //dimension du support dans Fourier
  //Nxmax=round(rayon*p_APERTURE/p_OIL_INDEX)=round(rayon*sin theta_max);
  double resultat_inter = p_OIL_INDEX * PIXEL_SIZE2 * map_dims["window_dim_x"] / (p_FACTORg * p_LAMBDA) * p_APERTURE/p_OIL_INDEX;
  const int Nxmax = round( p_FOCAL_RATIO * resultat_inter );
  // NXMAX=round(n1*Rf*Tps*DIMX_CCD2/(G*lambda)*NA/n1)
  // NXMAX Mat HA: round(n1 * Rf * Tps * DIMX_CCD2 / (G*lambda) * NA/n1): pareil

  /// ajouts ici pour passage en cubique
  const size_t dim_final = 512;
  double zoom = double(dim_final) / (4.0f * p_FOCAL_RATIO * resultat_inter);
  //double zoom=double(dim_final)/double(4*n1*Rf*Tps*DIMX_CCD2/(G*lambda)*NA/n1);

  // dimension de l'espace final
  // ici, arbitrairement, on short-cute le calcul de dimension pour forcer 512^3 
  // (4 * 128)^3 = 64 * 128 * 128 * 128, cf N_tab
  int Nxmax_Rf = 128; //round( zoom * resultat_inter );
  //NXMAX_Rf=round(zoom*n1*Tps*DIMX_CCD2/(G*lambda)*NA/n1);
  cout << endl << "dimension de sortie: " << Nxmax_Rf; cout.flush();

  const int Nymax = Nxmax;

  printf("\nRayon %f \n",rayon);
  printf("\ndelta_zmax %f \n",delta_zmax);
  printf("Nxmax=%i,4*Nxmax=%i \n",Nxmax,4*Nxmax);

  // Taille de l'arête d'un volume destination (real ou imag). A été forcé à 512
  size_t cube_edge = 4 * Nxmax_Rf;
  //N_tab taille du tableau 1D dans lequel on stocke un volume (real ou imag) en row-major order
  int N_tab = 64 * Nxmax_Rf * Nxmax_Rf * Nxmax_Rf;
  int cube_nbvoxels = N_tab; /// alias

  printf("N_tab (indice espace 3D)= %i \n",N_tab);
  cout <<"Espace final : " << 4 * Nxmax_Rf << " pixels cube\n";

  
  const int dv0= round(4 * Nxmax); //dimension totale en x,y
  const int dv0rf= round(4 * Nxmax_Rf);
  //const int cube_edge = dv0rf; /// alias. 4 * 128 (ou 127)


  map_experiment["cube_edge"] = (float)cube_edge;
  map_experiment["nxmax"] = (float)Nxmax;
  map_experiment["nxmax_rf"] = (float)Nxmax_Rf;
  map_experiment["rayon"] = rayon;
  map_experiment["delta_zmax"] = delta_zmax;
}
