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
      printf("chemin d'un volume AIR requis");
      printf("écrit des merdes dans /ramdisk");
      
      exit(EXIT_FAILURE);
    }
}


/* -------------------------------------------------------------------------- */
// 
/* -------------------------------------------------------------------------- */



int main(int argc, char *argv[])
{
  usage(argc, argv);
  
  char* AIR_FILE = argv[1];

  
  // nombre de threads employés sous fftw en TF3D
  int g_fftw_threads;


  // ***************************************************************************
  //
  // Commencement des calculs
  //
  // ***************************************************************************

  
  size_t i, j, k;


  
  // initialize multithreaded fftw 
  if (g_fftw_threads > 1)
    assert(fftw_init_threads());
  


  AIR_Volume<unsigned short int> V_sup_redon(AIR_FILE);
  V_sup_redon.set_data_linear_mode(true);
  V_sup_redon.update_definitions(); 
  V_sup_redon.allocate();
  V_sup_redon.update_data();
  V_sup_redon.binarize_data(0);


  // AIR_Volume<unsigned short int> V_supred_cubic(); // BUG! MARCHE PAS!
  AIR_Volume<unsigned short int> V_supred_cubic(0, 0, 0);
  V_supred_cubic.set_data_linear_mode(true);
  V_supred_cubic.get_definitions_from(V_sup_redon);
  V_supred_cubic.change_files("/ramdisk/cubic_bin_redon");
  V_supred_cubic.set_x_dim(512);
  V_supred_cubic.set_y_dim(512);
  V_supred_cubic.set_z_dim(512);
  V_supred_cubic.allocate();
  V_supred_cubic.fill(1);
  V_supred_cubic.copy_from_smaller(V_sup_redon);
  V_supred_cubic.write_hdr();
  V_supred_cubic.write_data();

  
  AIR_Volume<unsigned short int> V_sr_cubeshift(0, 0, 0);
  V_sr_cubeshift.set_data_linear_mode(true);
  V_sr_cubeshift.get_definitions_from(V_supred_cubic);
  V_sr_cubeshift.change_files("/ramdisk/cubic_shifted");
  V_sr_cubeshift.allocate();
  //V_sr_cubeshift.fill(30);
  
  V_supred_cubic.circshift_to(V_sr_cubeshift);
  V_sr_cubeshift.write_hdr();
  V_sr_cubeshift.write_data();


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
  //IntVol.set_bits(16);  // traitre si on oublie (assert dans write_hdr) // géré par templates
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
  V_sup_redon.change_files("/ramdisk/toto22"); //non, 22!
  V_sup_redon.write_hdr();
  V_sup_redon.write_data();

  
  
  //--------------------------------------------------
  // Test copies volume
    
  AIR_Volume<unsigned short int> VCube12(12, 12, 12, 1);
  VCube12.set_data_linear_mode(true);
  VCube12.allocate();
  //
  AIR_Volume<unsigned short int> VCube8(8, 8, 8, 1);
  VCube8.set_data_linear_mode(true);
  VCube8.allocate();
    
  for (i = 0; i < (8 * 8 * 8); i++)
    VCube8.data(i) = i;

  cout << endl;
  //  VCube8.show_values(cout);

  cout << endl << endl;
	

  VCube12.copy_from_smaller(VCube8);
  //VCube12.show_values(cout);
  
  //--------------------------------------------------


  
    


  //unsigned short int *sup_redon;
  //ARRAY_ALLOC(sup_redon, N_tab, unsigned short int);
  

    





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



