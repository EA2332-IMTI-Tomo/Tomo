#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>


using namespace std;



#include "main.h"
#include "IO.h"



void 
write_butterfly(RECON_TYPE* reel_arc, int tab_size, char* filename)
{
  BYTE* papillon_masque;
  ARRAY_NEW(papillon_masque, tab_size, BYTE);

  for (int i = 0; i < tab_size; i++)
    papillon_masque[i] = ((reel_arc[i] == 0) ? 0 : 1);
  
  FILE* fichier_papillon_masque = NULL;
  fichier_papillon_masque = fopen(filename, "wb");
  ASSERT(fichier_papillon_masque);
 
  fwrite(papillon_masque, sizeof(BYTE), tab_size, fichier_papillon_masque);

  fclose(fichier_papillon_masque);
  delete(papillon_masque);
}



void 
write_butterfly_orig(RECON_TYPE* reel_arc, int tab_size, char* filename)
{
  RECON_TYPE *papillon_masque=new RECON_TYPE[tab_size];
  float precision;
  

  for (int i = 0; i < tab_size; i++)
    papillon_masque[i] = ((reel_arc[i] == 0) ? 0 : 1);


  FILE* fichier_papillon_masque = NULL;
  fichier_papillon_masque = fopen(filename, "wb");
  ASSERT(fichier_papillon_masque);


  for(int cptb=0;cptb<tab_size;cptb++)
    {
      //Ecriture sur le disque
      precision=papillon_masque[cptb];
      fwrite(&precision,sizeof(precision),1,fichier_papillon_masque);
      //cout<<"\npapillon_masque(cpt)="<<papillon_masque[cptb];
    }

  
  fclose(fichier_papillon_masque);
  delete[] papillon_masque;
}


void
convert_write_results(int tab_size, char* reel_filename, char* imag_filename, char* module_filename, RECON_TYPE* final_reel_shift, RECON_TYPE* final_imag_shift, double* valeur_module_shift, bool skip_module)
{
  
  FILE* fichier_final_reel_shift = NULL;
  FILE* fichier_final_imag_shift = NULL;
  FILE* fichier_final_modul_shift = NULL;
  
  fichier_final_reel_shift = fopen(reel_filename, "wb");
  fichier_final_imag_shift = fopen(imag_filename, "wb");
  fichier_final_modul_shift = fopen(module_filename, "wb");

  ASSERT(fichier_final_reel_shift && fichier_final_imag_shift && fichier_final_modul_shift);
	  
  int elem_size = sizeof(OUTPUT_TYPE);

  CHRONO_SET(transfert);

  CHRONO_START(transfert);
  OUTPUT_TYPE* reel;
  ARRAY_ALLOC(reel, tab_size, OUTPUT_TYPE);
#pragma omp parallel for
  for (int i = 0; i < tab_size; i++)
    reel[i] = OUTPUT_TYPE(final_reel_shift[i]);
  CHRONO_STOP(transfert);
  //free(final_reel_shift);
  fwrite(reel, elem_size, tab_size, fichier_final_reel_shift);
  //free(reel);
	  
	  
  OUTPUT_TYPE* imag;
  CHRONO_START(transfert);
  imag = reel; //ARRAY_ALLOC(imag, tab_size, OUTPUT_TYPE);
#pragma omp parallel for
  for (int i = 0; i < tab_size; i++)
    imag[i] = OUTPUT_TYPE(final_imag_shift[i]);
  CHRONO_STOP(transfert);
  //free(final_imag_shift);
  fwrite(imag, elem_size, tab_size, fichier_final_imag_shift);
  //free(imag);


  if (! skip_module)
    {

      OUTPUT_TYPE* module;
      CHRONO_START(transfert);
      module = imag; // ARRAY_ALLOC(module, tab_size, OUTPUT_TYPE);
#pragma omp parallel for
      for (int i = 0; i < tab_size; i++)
	module[i] = OUTPUT_TYPE(valeur_module_shift[i]);
      CHRONO_STOP(transfert);
      // free(valeur_module_shift);
      fwrite(module, elem_size, tab_size, fichier_final_modul_shift);
      free(module);
      fclose(fichier_final_modul_shift);
    }


  CHRONO_END(transfert, "temps passé à allouer et à copier en mémoire");

	  
  printf("ecriture de final_reel_shift.bin, final_imag_shift.bin, final_modul_shift : OK \n");
  fclose(fichier_final_reel_shift);
  fclose(fichier_final_imag_shift);

}



// méthode qui était là avant et qui est un petit peu plus lente quand même
void
write_results2(int elem_size, int tab_size, char* reel_filename, char* imag_filename, char* module_filename, RECON_TYPE* final_reel_shift, RECON_TYPE* final_imag_shift, RECON_TYPE* valeur_module_shift)
{
  
  FILE* fichier_final_reel_shift = NULL;
  FILE* fichier_final_imag_shift = NULL;
  FILE* fichier_final_modul_shift = NULL;
  
  fichier_final_reel_shift = fopen(reel_filename, "wb");
  fichier_final_imag_shift = fopen(imag_filename, "wb");
  fichier_final_modul_shift = fopen(module_filename, "wb");

  OUTPUT_TYPE precision;
  
  for(int cpt = 0; cpt < (tab_size); cpt++)
    {
      precision = final_reel_shift[cpt];
      fwrite(&precision, sizeof(precision), 1, fichier_final_reel_shift);

      precision = final_imag_shift[cpt];
      fwrite(&precision, sizeof(precision), 1, fichier_final_imag_shift);

      precision = valeur_module_shift[cpt];
      //ecriture du module
      fwrite(&precision, sizeof(precision), 1, fichier_final_modul_shift);
    }

	  
  printf("ecriture de final_reel_shift.bin, final_imag_shift.bin, final_modul_shift : OK \n");
  fclose(fichier_final_reel_shift);
  fclose(fichier_final_imag_shift);
  fclose(fichier_final_modul_shift);
}



/*
void 
write_centers(int tab_size, int* tab, char* filename)
{
  FILE* fichier = fopen(filename, "wb");
  ASSERT(fichier != NULL);
	  
  int elem_size = sizeof(int);
  fwrite(tab, elem_size, tab_size, fichier);

  fclose(fichier);
}



void 
write_supRedon(int tab_size, SUPREDON_OUTPUT_TYPE* tab, char* filename)
{
  FILE* fichier = fopen(filename, "wb");
  ASSERT(fichier != NULL);
  	  
  int elem_size = sizeof(SUPREDON_OUTPUT_TYPE);
  fwrite(tab, elem_size, tab_size, fichier);


  fclose(fichier);
}
*/


/* déplacé en .h */
/*

  template <typename data_type, typename output_type>
  void
  write_convert_array(size_t tab_size, data_type* tab, char* filename)
  {
  FILE* fichier = fopen(filename, "wb");
  assert(fichier != NULL);
  	  
  int elem_size = sizeof(output_type);

  output_type* converted;
  ARRAY_ALLOC(converted, tab_size, output_type);

  for (int i = 0; i < tab_size; i++)
  converted[i] = (output_type) tab[i];
  
  fwrite(converted, elem_size, tab_size, fichier);
  free(converted);
  fclose(fichier);
  }



  template <typename data_type>
  void
  write_array(size_t tab_size, data_type* tab, char* filename)
  {
  FILE* fichier = fopen(filename, "wb");
  assert(fichier != NULL);
  	  
  int elem_size = sizeof(data_type);  
  fwrite(tab, elem_size, tab_size, fichier);

  fclose(fichier);
  }


*/


//if (g_save_sup_redon)
//  write_array<SUPREDON_OUTPUT_TYPE>(N_tab, sup_redon, str_concat(g_OUTPUT_DIR, OUTPUT_REDON_FILENAME));

	  
	// sizeof est compile-time, et precision est de type float en ce moment
	//int data_size = sizeof(OUTPUT_TYPE);
	  
	// conversion et sauvegarde sur disque
	// (passage de 64 à 32 bits pour gain de place)
	// 	convert_write_results(N_tab, str_concat(g_OUTPUT_DIR, OUTPUT_REAL_FILENAME), str_concat(g_OUTPUT_DIR, OUTPUT_IMAG_FILENAME), str_concat(g_OUTPUT_DIR, OUTPUT_MOD_FILENAME), final_reel_shift, final_imag_shift, valeur_module_shift, !(g_save_module && g_compute_module) );
