#ifndef __IO__
#define __IO__


// convertit en mémoire les résultats finaux en OUTPUT_TYPE (essentiellement de double en float, pour ne pas prendre trop de place sur le disque)
// alloue donc une instance supplémentaire en mémoire
// écrit les volumes sur disque
void
convert_write_results(int tab_size, char* reel_filename, char* imag_filename, char* module_filename, RECON_TYPE* final_reel_shift, RECON_TYPE* final_imag_shift, double* valeur_module_shift, bool skip_module);


// export de la fonction originale
void
write_results2(int elem_size, int N3D, char* reel_filename, char* imag_filename, char* module_filename, RECON_TYPE* final_reel_shift, RECON_TYPE* final_imag_shift, double* valeur_module_shift);

// création d'un papillon binaire à partir de réel arc et sauvegarde sur disque
void 
write_butterfly(RECON_TYPE* reel_arc, int tab_size, char* filename);


void 
write_butterfly_orig(RECON_TYPE* reel_arc, int tab_size, char* filename);


void 
write_centers(int tab_size, int* tab, char* filename);


void 
write_supRedon(int tab_size, unsigned short int* tab, char* filename);



template <typename data_type>
void 
write_array(size_t tab_size, data_type* tab, char* filename);


template <typename data_type, typename output_type>
void 
write_convert_array(size_t tab_size, data_type* tab, char* filename);


// ---------------------------------------------------------------------------


template <typename data_type, typename output_type>
void
write_convert_array(size_t tab_size, data_type* tab, char* filename)
{
  FILE* fichier = fopen(filename, "wb");
  ASSERT(fichier != NULL);
  	  
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
  ASSERT(fichier != NULL);
  	  
  int elem_size = sizeof(data_type);  
  fwrite(tab, elem_size, tab_size, fichier);

  fclose(fichier);
}



#endif // __IO__
