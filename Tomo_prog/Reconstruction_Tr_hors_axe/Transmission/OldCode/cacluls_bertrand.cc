#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>
#include <assert.h>
#include <exception>


using namespace std;

#include "main.h"
#include "util.h"
//#include "readloop.h"
#include "memory.h"
#include "IO.h"

#include "fourier.h"
#include "fourier_shifted.h"

#include "AIR_Volume.h"



#define LINESIZE 512


/* -------------------------------------------------------------------------- */
// 
/* -------------------------------------------------------------------------- */



template <typename U>
U
fun_double(U numeric_var);



template <typename U>
U
fun_double(U num_val)
{
  return 2 * num_val;
}


template <typename U>
U
fun_add(U left, U right)
{
  return left + right;
}



/* -------------------------------------------------------------------------- */
// 
/* -------------------------------------------------------------------------- */


inline double
nozero(double i)
{
  return ((i == 0) ? 1 : i);
}



/* -------------------------------------------------------------------------- */
// Usage 
/* -------------------------------------------------------------------------- */


static void
usage(int argc, char **argv) 
{

  if ((argc - 1) == 0)
    {
      printf("Tests hypothèse Bertrand: calculs de post-processing sur l'image module et le papillon binarisé\n");
      printf("on relit l'image module produite par le logiciel de reconstruction en Transmission, de même que le papillon\n");
      printf("on produit une seule image résultat module améliorée");
      printf("Paramètres: on reprend exactement le fichier de configuration utilisé pour la reconstruction, sans modification.\n");
      printf("CUDA n'est pas géré du fait de volumes de taille doublée\n");
      printf("\n les paramètres ci-dessous sont rapellés à titre indicatif \n ");
      printf("\n \n ********************** ");


      printf("Values for I/O directories and processing parameters are read from mandatory ./config.txt file and can be redefined as follow:\n");
      printf("mandatory: <file> configuration file \n");
      printf("opt: -input <dir> directory holding input images \n");
      printf("opt: -output <dir> directory where results are written \n");
      printf("opt: -masks <dir> directory holding masks \n");
      printf("opt: -in_radix <name> input images filename radix (program adds image number + extension) \n");
      printf("opt: -fftw_threads <int> number of threads to activate when computings ffts with fftw software library\n");      
      printf("opt: -image_dimx(y) <int> to define the size of the input images\n");      
      printf("opt: -ccd_dimx(y) <int> to define the size of the frame to process on aforementioned images \n");      
      printf("opt: -edge_x(y) <int> to define the upper left corner position of aforementioned frame \n");      
      printf("opt: -x(y)m0_limit <int> to set exclusion for some centers \n");      
      printf("opt: -first_lane <int> \n");      
      printf("opt: -final_angle <int> \n");      
      printf("opt: -inc_angle <int> \n");      
      printf("opt: -oil_index <float> \n");      
      printf("opt: -aperture <float> \n");  
      printf("opt: -lambda <float> to set wavelength \n");  
      printf("opt: -factor <int> pour définir grossissement du telan+objectif \n");  
      printf("opt: -pixel_size <float> \n");  
      printf("opt: -rf_val <float> \n");  
      
      exit(EXIT_FAILURE);
    }
}


/* -------------------------------------------------------------------------- */
// 
/* -------------------------------------------------------------------------- */



int main(int argc, char *argv[])
{
  usage(argc, argv);
  

  //* int image_dimx, image_dimy; temporairement global
  // valeur controlant l'exclusion de certains centres
  int xm0_limite; //centre maximum
  int ym0_limite;
  //valeur du coin pour la découpe
  int coin_x, coin_y;
  int premier_plan;
  // limite de lecture des fichiers d'entrée
  int Num_Angle_final; 
  // int sans_signal=0;
  int SautAngle;
  float p_OIL_INDEX;	//indice de l'huile // n1
  float p_APERTURE;	//ouverture numerique de l'objectif? //NA
  // (celle du condenseur intervient sur la forme, la taille, du papillon)
  float p_LAMBDA; //632 * pow(10, -9); //longueur d'onde
  int p_FACTORg; //G; // grossissement telan+objectif
  float p_PIXEL_SIZE; // 11.2 * pow(10, -6); //Taille des pixels
  float p_FOCAL_RATIO; // 1.2; 
  // float p_FOCAL_RATIO=0.8;//1.7;//0.77;//1.3;// bonne valeur le 29/9; //rapport f1/

  
  // ==========================================================
  // 
  // ==========================================================

  
  // résolution de la caméra ccd 
  // résolution de l'image reconstruite
  int g_window_dim_x;
  int g_window_dim_y;

  // nombre de threads employés sous fftw en TF3D
  int g_fftw_threads;

  char* g_MASKPATH;
  char* g_OUTPUT_DIR;
  char* g_INPUT_DIR;
  char* g_INFILE_RADIX; 

  
  // à des fins de déboguage
  int image_dimx, image_dimy;
  bool read_saved_volume;


  
  // ==========================================================
  // set parameters alike parsed configuration file
  // ==========================================================


  ifstream cfg_file (argv[1]);
  argc -= 2; argv += 2;

  assert (cfg_file.is_open());

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
	CASE_STR("CCD_DIMX") g_window_dim_x = atoi(value_name);
      CASE_STR("CCD_DIMY") g_window_dim_y = atoi(value_name);
      CASE_STR("XM0_LIMIT") xm0_limite = atoi(value_name);
      CASE_STR("YM0_LIMIT") ym0_limite = atoi(value_name);
      CASE_STR("EDGE_X") coin_x = atoi(value_name);
      CASE_STR("EDGE_Y") coin_y = atoi(value_name);
      CASE_STR("IMAGE_DIMX") image_dimx = atoi(value_name);
      CASE_STR("IMAGE_DIMY") image_dimy = atoi(value_name);
      CASE_STR("FIRST_LANE") premier_plan = atoi(value_name);
      CASE_STR("FINAL_ANGLE") Num_Angle_final = atoi(value_name);
      CASE_STR("INC_ANGLE") SautAngle = atoi(value_name);
      CASE_STR("OIL_INDEX") p_OIL_INDEX = atof(value_name);
      CASE_STR("APERTURE") p_APERTURE = atof(value_name);
      CASE_STR("LAMBDA") p_LAMBDA = atof(value_name);
      CASE_STR("FACTOR") p_FACTORg = atoi(value_name);
      CASE_STR("PIXEL_SIZE") p_PIXEL_SIZE = atof(value_name);
      CASE_STR("RF_VAL") p_FOCAL_RATIO = atof(value_name);
      CASE_STR("FFTW_THREADS") g_fftw_threads = atoi(value_name);
      CASE_STR("READ_DUMPED_VOLUME") read_saved_volume = (bool)atoi(value_name); 
      CASE_STR("MASKPATH") g_MASKPATH = str_alloc_cpy(value_name);
      CASE_STR("OUTPUT_DIR") g_OUTPUT_DIR = str_alloc_cpy(value_name);
      CASE_STR("INPUT_DIR") g_INPUT_DIR = str_alloc_cpy(value_name);
      CASE_STR("INFILE_RADIX") g_INFILE_RADIX = str_alloc_cpy(value_name);
      SWITCH_STR_ELSE_END( 
			  cerr << endl << "parse error: "		\
			  << variable_name << " | " << value_name;  
			   )
	}
  cfg_file.close();
  

  // ==========================================================
  // override parameters by command-line arguments
  // ==========================================================

  
  while (argc > 0)
    {
      variable_name = argv[0];
      if (!(argc - 1)) break;
      value_name = argv[1];

      SWITCH_STR_BEGIN( variable_name )
	CASE_STR("-ccd_dimx") g_window_dim_x = atoi(value_name);
      CASE_STR("-ccd_dimy") g_window_dim_y = atoi(value_name);
      CASE_STR("-xm0_limit") xm0_limite = atoi(value_name);
      CASE_STR("-ym0_limit") ym0_limite = atoi(value_name);
      CASE_STR("-edge_x") coin_x = atoi(value_name);
      CASE_STR("-edge_y") coin_y = atoi(value_name);
      CASE_STR("-image_dimx") image_dimx = atoi(value_name);
      CASE_STR("-image_dimy") image_dimy = atoi(value_name);
      CASE_STR("-first_lane") premier_plan = atoi(value_name);
      CASE_STR("-final_angle") Num_Angle_final = atoi(value_name);
      CASE_STR("-inc_angle") SautAngle = atoi(value_name);
      CASE_STR("-oil_index") p_OIL_INDEX = atof(value_name);
      CASE_STR("-aperture") p_APERTURE = atof(value_name);
      CASE_STR("-lambda") p_LAMBDA = atof(value_name);
      CASE_STR("-factor") p_FACTORg = atoi(value_name);
      CASE_STR("-pixel_size") p_PIXEL_SIZE = atof(value_name);
      CASE_STR("-rf_val") p_FOCAL_RATIO = atof(value_name);
      CASE_STR("-fftw_threads") g_fftw_threads = atoi(value_name);
      CASE_STR("-masks") g_MASKPATH = value_name;
      CASE_STR("-input") g_INPUT_DIR = value_name;
      CASE_STR("-output") g_OUTPUT_DIR = value_name;
      CASE_STR("-in_radix") g_INFILE_RADIX = value_name;
      SWITCH_STR_ELSE_END(
			  cerr << endl; 
			  usage(argc, argv);
			  exit (EXIT_FAILURE);
			  );

      argv += 2; argc -= 2;
    }



    
  // ==========================================================
  // perform computations depending from parameter values
  // ==========================================================


  if(coin_x + g_window_dim_x > image_dimx || coin_y + g_window_dim_y > image_dimy)
    printf("Attention, la zone decoupee sort de l'image\n");
  //largeur : 740, hauteur : 574.

  //short int precision;
  float precision;

  int NbAngle = Num_Angle_final - premier_plan + 1;

  // coordonnees du coin haut gauche dans l'image originale (748x570) 
  // pour l'image recadrée --------

  float theta = asin(p_APERTURE/p_OIL_INDEX);//theta_max
  printf("theta : %f\n",theta);

  // chemin de lecture fichiers entrée + radix de fichier (+ numfic.bmp)
  //   char Chemin[] = INPUT_DIR "/" INFILE_RADIX ;
  char *Chemin = str_concat(str_concat(g_INPUT_DIR, (char *)"/"),
			    g_INFILE_RADIX);
  
  // Facteur d'echelle
  float K_SCALE = p_LAMBDA * p_FACTORg/(2 * p_APERTURE * p_PIXEL_SIZE * p_FOCAL_RATIO);
  printf("K_SCALE : %f\n", K_SCALE);
  float PIXEL_SIZE2 = p_PIXEL_SIZE;//*K_SCALE; //Taille pixel apres mise à l'echelle important si RF<>1

  int Ntx = g_window_dim_x;//K_SCALE; //Nombre total de pixel en x
  /*int Nty=Ntx;*/ //(support carré)

  float rayon = p_OIL_INDEX * p_FOCAL_RATIO * PIXEL_SIZE2 * g_window_dim_x/(p_FACTORg * p_LAMBDA); //Rayon

  float delta_zmax = rayon * cos(theta);//float?
  //int rayon=round(p_OIL_INDEX*p_FOCAL_RATIO*p_PIXEL_SIZE*Ntx/(G*p_LAMBDA));ici un commentaire archéologique
  //int rayon=100;

  //dimension du support dans Fourier
  //Nxmax=round(rayon*p_APERTURE/p_OIL_INDEX)=round(rayon*sin theta_max);
  const int Nxmax = round(p_OIL_INDEX * p_FOCAL_RATIO * PIXEL_SIZE2 * g_window_dim_x / (p_FACTORg * p_LAMBDA) * p_APERTURE/p_OIL_INDEX);
  int Nxmax_Rf = round(p_OIL_INDEX * PIXEL_SIZE2 * g_window_dim_x / (p_FACTORg * p_LAMBDA) * p_APERTURE/p_OIL_INDEX);
  //int Nxmax_Rf = round(p_OIL_INDEX * PIXEL_SIZE2 * g_window_dim_x / (G*p_LAMBDA) * p_APERTURE/p_OIL_INDEX) + 1; // trick to test 2^n

  //   float test_nxmax_Rf = p_OIL_INDEX * PIXEL_SIZE2 * g_window_dim_x / (p_FACTORg * p_LAMBDA) * p_APERTURE/p_OIL_INDEX;
  //   printf("testnxmax_Rf %f \n", test_nxmax_Rf);
  //   printf("round(testnxmax_Rf) %f \n",round(test_nxmax_Rf));
  // 
  const int Nymax = Nxmax;

  printf("Rayon %f \n",rayon);
  printf("Nxmax=%i,4*Nxmax=%i \n",Nxmax,4*Nxmax);
  //Indice crée pour balayer le tableau final 1D (informatique) 
  //contenant le volume final de données en 3D
  int N_tab=64*Nxmax_Rf*Nxmax_Rf*Nxmax_Rf;
  int cube_nbvoxels = N_tab; /// alias

  printf("N_tab (indice espace 3D)= %i \n",N_tab);
  cout <<"Espace final : " << 4 * Nxmax_Rf << " pixels cube\n";

  
  const int dv0=4*Nxmax; //dimension totale en x,y
  const int dv0rf=4*Nxmax_Rf;
  const int cube_edge = dv0rf; /// alias. 4 * 128 (ou 127)






  // ***************************************************************************
  //
  // Commencement des calculs
  //
  // ***************************************************************************

  size_t i;
  


  
  // initialize multithreaded fftw 
  if (g_fftw_threads > 1)
    assert(fftw_init_threads());
  

  
  /*


    AIR_Volume<unsigned short int> V_sup_redon("/ramdisk/sup_redon_C.hdr");
    V_sup_redon.set_data_linear_mode(true);
    V_sup_redon.update_definitions(); 
    V_sup_redon.allocate();
    V_sup_redon.update_data();
    V_sup_redon.binarize_data(0);

    //V_sup_redon.map(&fun_double);
    //V_sup_redon.map(&fun_double);
    //V_sup_redon.double_values();
    //V_sup_redon.double_values();
    V_sup_redon.scalar_mult(5);
    V_sup_redon.fill(8);

    V_sup_redon.change_files("/ramdisk/sup_redon_C_binary");
    V_sup_redon.write_data();
    V_sup_redon.write_hdr();
  


    AIR_Volume<float> V_sup_redonD(508, 508, 508, 1);
    //V_sup_redonD.set_bits(32); 
    V_sup_redonD.set_data_linear_mode(true);
    V_sup_redonD.allocate();
    V_sup_redonD.copy_from(V_sup_redon);
    V_sup_redonD.change_files("/ramdisk/float_redon_bin");
    V_sup_redonD.write_hdr();
    V_sup_redonD.write_data();



    AIR_Volume<unsigned short int> IntVol(508, 508, 508, 1);
    //IntVol.set_bits(16);  // traitre si on oublie (assert dans write_hdr) 
    IntVol.set_data_linear_mode(true);  // traitre si on oublie (asserts plus tard) 
    IntVol.allocate();
    float* float_data = V_sup_redonD.get_data_linear();
    //ARRAY_ALLOC(float_data, 508 * 508 * 508, float);
    IntVol.copy_from(float_data);
    IntVol.change_files("/ramdisk/sizet_redon_bin");
    IntVol.write_data();
    IntVol.write_hdr();

  

    IntVol.fill(5);
    IntVol+=2;
    IntVol*=2; // vaut 14
    IntVol.map2_star(&fun_add, V_sup_redon); // on ajoute celui à 8 -> 22

    V_sup_redon = IntVol;
    V_sup_redon.change_files("/ramdisk/toto14"); //non, 22!
    V_sup_redon.write_hdr();
    V_sup_redon.write_data();
  
  
    //unsigned short int *sup_redon;
    //ARRAY_ALLOC(sup_redon, N_tab, unsigned short int);
  

    */





  // on binarise sup_redon
  //AIR_Volume<double> V_sup_redon_binaire();

  /*
    try{
    sup_redon_binaire = new double[N_tab]; 
    }
    catch(std::bad_alloc& ex)
    {cerr << endl << "memory allocation failure" << ex.what(); exit(EXIT_FAILURE); } 
    catch(...)
    {cerr << endl << "other error"; exit(EXIT_FAILURE); } 
  */
  
  /*
    unsigned short int val_sr;
    for(int cpt = 0; cpt < N_tab; cpt++)
    {
    val_sr = V_sup_redon[i];
    V_sup_redon_binaire[i] = (double)((val_sr > 0) ? 1 : 0);
    }
  */

  //TF3D_cuda_r2c_backward()


  return(EXIT_SUCCESS);
}



