// =============================================================================
// =============================================================================
// fonctions-outil apellées à chaque angle
// =============================================================================
// =============================================================================



// =============================================================================
// =============================================================================

// fonction outil centrale de batch_launch, destinée à traiter un hologramme et à sauver le résultat dans fft_reel_shift / fft_imag_shift


void
Holo_Process::hologram_compute(size_t current_angle)
{
  
      //----------------------------------------------------------------------
      ///////////////////////// phase shifting ou pas ////////////////////////
      //----------------------------------------------------------------------

      // chargement selon fichier ou mémoire selon ce qui est sélectionné

      // cette opération, avec ou sans HA, transforme les (1/4)
      // images réelles en 1 hologramme complexe

      if (p_psh_mode) // PSH pur ou bien PSH+HA: même calcul
	phase_shifting(current_angle, plan_reel, plan_imag);
      else // HA pur
	no_phase_shifting(current_angle, plan_reel, plan_imag); 

      //show_hdp(); //cvWaitKey(500);


      //----------------------------------------------------------------------
      /////////////// TF et crop dans les fréquences instrumentales //////////
      //----------------------------------------------------------------------


      // --------------------------------------------------
      // A) TF

      // cas particuler HA pur
      if (p_off_axis_mode && ! p_psh_mode) 
	VF_Fourier2D -> import_R_from(plan_reel, true);
      else
	VF_Fourier2D -> import_from(plan_reel, plan_imag, true);
      

      // calcul
      VF_Fourier2D -> set_fourier_forward();
      // à rajouter selon Mat (5/03/15)
      VF_Fourier2D -> set_fourier_normalize();

      // de-tangling et re-circshift
      VF_Fourier2D -> export_to(plan_reel, plan_imag, true); 



      // --------------------------------------------------
      // B) Crop dans les fréquences instrumentales
      
      
      // HDP: les données sont alors dans fft_reel/imag_shift, qui ne font plus que  
      // 212*212 au lieu de 512^2
      
      // A) HDP pur: découpage carré selon Nxmax calculé
      if ( ! p_off_axis_mode ) 
	{
	  // 0.06s les deux en cumulé
	  square_crop(plan_reel, fft_reel_shift, p_window_size, 2 * p_Nxmax);
	  square_crop(plan_imag, fft_imag_shift, p_window_size, 2 * p_Nxmax);
	} 
      else // B) HA avec ou sans HDP: découpage selon cercle mesuré en show_fourier (x,y,r)
	{
	  // découpage d'un cercle (un carré, en fait) défini par un centre et un rayon
	  disk_crop(plan_reel, fft_reel_shift, p_window_size, p_off_axis_circle_cx, p_off_axis_circle_cy, p_off_axis_circle_r);
	  disk_crop(plan_imag, fft_imag_shift, p_window_size, p_off_axis_circle_cx, p_off_axis_circle_cy, p_off_axis_circle_r);
	}

}



// =============================================================================
// =============================================================================
// Mapping 3D


// angle courant | peak_x | peak_y 
void
Holo_Process::angle_mapping(size_t pl_current_angle, int pl_xmi, int pl_ymi, double* in_fft_r, double* in_fft_i)
{
  if (! p_psh_mode) // hack pourri et redondant
    {
      p_rayon = p_Nxmax = p_Nymax = p_off_axis_circle_r;
    }

  // coordonnée dans l'image2D centrée (xm0,ym0)=(0,0)=au centre de l'image
  const int l_xm0 = pl_xmi - p_Nxmax;
  // printf("xm0 %i \n", xm0);
  const int l_ym0 = pl_ymi - p_Nymax;

  if(l_xm0 == 0 && l_ym0 == 0)
    printf("(xm0,ym0)=(0,0) dans le plan %d\n", pl_current_angle);
      

  if((sqr(l_xm0) + sqr(l_ym0)) > c_rayon_inf)
    {
      i_centres_exclus++; 
      // MSG_ASSERT(false, "lol, encore jamais arrivé!!"); // sur le hors-axe, si
    }
  else
    {
      // sauvegarde des centres en coordonnées non centrée; on met le numero d'angle
      i_centre[pl_xmi * 2 * p_Nxmax + pl_ymi] = pl_current_angle;
      // pour vérifier les pb. A ouvrir en sizeof(int)=32 bits mais non signé

      // indice 3D
      // variables locales uniquement
      int k = 0;
      double arg_z_arc = 0;
      double z_arc = 0;
      double zm0;
      double zm0_carre = sqr(p_rayon) - sqr(l_xm0) - sqr(l_ym0);

      if(round(zm0_carre) > -1)
	{
	  zm0 = sqrt( zm0_carre );
	  // printf("zm0:%i\n",zm0);
	  i_nb_proj++;

	  // 4ms de base par angle => <1 avec omp
	  cpu_kernel_mapping(l_xm0, l_ym0, p_Nxmax, p_Nxmax, \
			     p_rayon, \  
			     p_delta_zmax, zm0,				\
			     c_dv0s2rf, c_dv1s2rf, c_dv0rf, c_dv2s2rf, c_dv1xdv0rf, \
			     i_sup_redon, i_reel_arc, i_imag_arc,	\
			     in_fft_r, in_fft_i);

	} //fin if zm0>-1
    } //fin else xm0_limite
}


// =============================================================================
// =============================================================================


void 
cpu_kernel_mapping(int l_xm0, int l_ym0, size_t p_Nxmax, size_t p_Nymax, \
		   float p_rayon, float p_delta_zmax, double zm0,  \
		   size_t c_dv0s2rf, size_t c_dv1s2rf, size_t c_dv0rf, size_t c_dv2s2rf, size_t c_dv1xdv0rf, \
		   unsigned short int *i_sup_redon, RECON_TYPE *i_reel_arc, RECON_TYPE *i_imag_arc, \
		   double *fft_reel_shift_norm, double *fft_imag_shift_norm)
{
  const size_t c_Nxmax2 = 2 * p_Nxmax, c_Nymax2 = 2 * p_Nymax;
  const float c_rayon_sqr = p_rayon * p_rayon;
  const float c_delta_zmax_sqr = p_delta_zmax * p_delta_zmax;
  const int c_scan_size = 2 * p_Nxmax * 2 * p_Nymax;

  // size_t fait exploser la formule
  int dv0s2rf = c_dv0s2rf;
  int dv1s2rf = c_dv1s2rf;
  int dv0rf = c_dv0rf;
  int dv2s2rf = c_dv2s2rf;
  int dv1xdv0rf = c_dv1xdv0rf;

  // on balaie une image de -106 à 106 en x et en y (p_Nymax = 106)
  // s: parcours tableau 1D
  // i, j: position ramenée à une matrice [212,212]
  // x, y: position dans une image dont le centre est en pos° [0, 0]

#pragma omp parallel for num_threads( OMP_THREADS ) 
  for (int s = 0; s < c_scan_size; s++)
    {
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
	  k = round(( - l_xm0 + _x + dv0s2rf) + (- l_ym0 + _y + dv1s2rf) * dv0rf + \
		    (z_arc + dv2s2rf) * dv1xdv0rf);

	  ASSERT( k >= 0 );

	  i_sup_redon[k] += 1;//pour calculer le support
	  
	  i_reel_arc[k] += fft_reel_shift_norm[s];//pour calculer l'image
	  i_imag_arc[k] += fft_imag_shift_norm[s];//pour calculer l'image
	}
    }
  	  
}


// =============================================================================
// =============================================================================
// Elimination du jumeau


void 
Holo_Process::twin_eliminate(int pl_xmi, int pl_ymi, double* io_fft_r, double* io_fft_i)
{
  const int l_xc = pl_xmi - p_Nxmax;
  const int l_yc = pl_ymi - p_Nymax;

  
      
  ////////////////////////////////////Virer l'objet jumeaux (0.001 à 5 s cumulé sur CPU)
  if(sqr(l_xc) + sqr(l_yc) > 900)//35*35=1225;objet jumeau pas trop près de l'objet
    { 
      i_jumeau_elimine++;

      const size_t c_Nxmax2 = 2 * p_Nxmax;
      const size_t c_Nymax2 = 2 * p_Nymax;

      //coordonnée objet jumeaux
      const int l_xmij = c_Nxmax2 - pl_xmi;
      const int l_ymij = c_Nymax2 - pl_ymi;


      const size_t c_mask_halfsize = c_mask_size / 2;
      double c_inv_255 = 1.0f / 255.0f;
	

      int borne_x = l_xmij + c_mask_halfsize;
      int borne_y = l_ymij + c_mask_halfsize;
      for(int x = l_xmij - c_mask_halfsize; x < borne_x; x++)
	{
	  for(int y = l_ymij - c_mask_halfsize; y < borne_y; y++)
	    {

	      //attention le masque repasse du coté gauche lorsque le jumeau touche le bord droit! A corriger
	      int x_2 = x;
	      int y_2 = y;
	      if(x > c_Nxmax2)
		x_2 = x - c_Nxmax2;
	      if(y > c_Nymax2)
		y_2 = y - c_Nymax2;
	      if(x < 0)
		x_2 = c_Nxmax2 + x;

	      if(y < 0)
		y_2 = c_Nymax2 + y;

	      int cptj = c_Nymax2 * y_2 + x_2;
	      double valeur_jumeau_sur_255 = cache_jumeau[ c_mask_size * (y - l_ymij + c_mask_halfsize) + (x - l_xmij+ c_mask_halfsize) ] * c_inv_255;
	      io_fft_r[cptj] *= valeur_jumeau_sur_255;
	      io_fft_i[cptj] *= valeur_jumeau_sur_255;
	    }//fin for y
	}//fin for x
    }//fin if objet jumeau pas trop pres

}


// =============================================================================
// =============================================================================
// Elimination du jumeau redux


void 
Holo_Process::twin_eliminate_redux(int pl_xmi, int pl_ymi, double* io_fft_r, double* io_fft_i, unsigned char* cache_jumeau)
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


      const size_t c_mask_halfsize = c_mask_size / 2; // 30/2 = 15
      double c_inv_255 = 1.0f / 255.0f;
	
      size_t borneinf_x = l_xmij - c_mask_halfsize;
      size_t bornesup_x = l_xmij + c_mask_halfsize;
      size_t borneinf_y = l_ymij - c_mask_halfsize;
      size_t bornesup_y = l_ymij + c_mask_halfsize;


#pragma omp parallel for num_threads( OMP_THREADS ) 
      for (size_t i = 0; i < c_NxNy4; i++)
	{
	  size_t y = i / c_Nxmax2;
	  size_t x = i % c_Nxmax2;
	    
	  if ( ( x > borneinf_x ) && ( x < bornesup_x ) &&  \	
	       ( y > borneinf_y ) && ( y < bornesup_y ) )
	    {
	      size_t index_dans_cache_jumeau = c_mask_size * ( y - l_ymij + c_mask_halfsize) + \
		x - l_xmij + c_mask_halfsize;
	      double valeur_jumeau_sur_255 = cache_jumeau[ index_dans_cache_jumeau ] * c_inv_255;
	      io_fft_r[i] *= valeur_jumeau_sur_255;
	      io_fft_i[i] *= valeur_jumeau_sur_255;
		
	    }
	}
		  
    }//fin if objet jumeau pas trop pres

}
    

// =============================================================================
// =============================================================================
// Recherche maximum et normalisation par pic


size_t
Holo_Process::peak_normalize(size_t crop_size, double* in_fft_r, double* in_fft_i, \
			     double* out_fft_r, double* out_fft_i)
{
  double max_part_reel, max_part_imag, max_module, max_mod_divisor;

  size_t cpt_max = 0;
  double max_val = 0; // valeur maximale du module des TFs, >0
      
  double tab[1];
                 
  // t(s) cpu par iter: 3e-5 s 
  // cumulé: 0.014s

  //shared(cpt_max, max_val) num_threads(4)
#pragma omp parallel for num_threads(4)
  for(size_t cpt = 1; cpt < crop_size; cpt++) 
    {
      double res = tab[0] = sqr(in_fft_r[cpt]) + sqr(in_fft_i[cpt]); 
	  
      if (res > max_val)
	{
	  cpt_max = cpt;
	  max_val = res; 
	}
    }
#pragma omp barrier     
  max_part_reel = in_fft_r[cpt_max];
  max_part_imag = in_fft_i[cpt_max];
  max_module = max_val; 
  max_mod_divisor = 1 / max_module;

      
  // t(s) cumulé cpu: 0.17s -> 0.04 à 0.07 avec *divisor // 0.04 omp
#pragma omp parallel for num_threads( 4 )
  for(size_t cpt = 0; cpt < crop_size; cpt++)
    {
      double _FFTRS = in_fft_r[cpt];
      double _FFTIS = in_fft_i[cpt];
      out_fft_r[cpt] = (_FFTRS * max_part_reel + _FFTIS * max_part_imag) * max_mod_divisor;
      out_fft_i[cpt] = (_FFTIS * max_part_reel - _FFTRS * max_part_imag) * max_mod_divisor;
    }
#pragma omp barrier     

  return cpt_max;
}


// =============================================================================
// =============================================================================
// Phase-shifting (appliqué en mode phase-shift seul ET combiné avec hors-axe)


void
Holo_Process::phase_shifting(size_t current_angle, double* plan_reel_out, double* plan_imag_out)
{
  // phase1 = holo1

  // Lecture des interférogrammes via fichiers
  if (! p_read_from_mem_mode )
    {
      if (i_holo_resize_needed_p) // il faut cropper les images pour découper une fenêtre carrée ^2
	{
	  compute_hologram_filename(p_images_radix, holo_filename, current_angle, 1);
	  i_vImgHoloReader -> reload_pgm(holo_filename);
	  vImgPhase1 -> copy_from_bigger(*i_vImgHoloReader, p_window_edge_x, p_window_edge_y); // vImgPhase1 et phase1 confondus
	  
	  compute_hologram_filename(p_images_radix, holo_filename, current_angle, 2);
	  i_vImgHoloReader -> reload_pgm(holo_filename);
	  vImgPhase2 -> copy_from_bigger(*i_vImgHoloReader, p_window_edge_x, p_window_edge_y); 

	  compute_hologram_filename(p_images_radix, holo_filename, current_angle, 3);
	  i_vImgHoloReader -> reload_pgm(holo_filename);
	  vImgPhase3 -> copy_from_bigger(*i_vImgHoloReader, p_window_edge_x, p_window_edge_y); 

	  compute_hologram_filename(p_images_radix, holo_filename, current_angle, 4);
	  i_vImgHoloReader -> reload_pgm(holo_filename);
	  vImgPhase4 -> copy_from_bigger(*i_vImgHoloReader, p_window_edge_x, p_window_edge_y); 
	}
      else // les images sont déjà acquises en carré
	{
	  compute_hologram_filename(p_images_radix, holo_filename, current_angle, 1); 
	  vImgPhase1 -> reload_pgm(holo_filename);
	  compute_hologram_filename(p_images_radix, holo_filename, current_angle, 2);
	  vImgPhase2 -> reload_pgm(holo_filename);
	  compute_hologram_filename(p_images_radix, holo_filename, current_angle, 3);  
	  vImgPhase3 -> reload_pgm(holo_filename);
	  compute_hologram_filename(p_images_radix, holo_filename, current_angle, 4);  
	  vImgPhase4 -> reload_pgm(holo_filename);
	}
    }
  else // lecture des interférogrammes via mémoire 
    {
      MSG_ASSERT( current_angle > 0, "erreur n° angle");
      
      size_t slice = (current_angle -1) * 4; 
      phase1 = (*i_V_hologramStack)[slice + 0];
      phase2 = (*i_V_hologramStack)[slice + 1];
      phase3 = (*i_V_hologramStack)[slice + 2];
      phase4 = (*i_V_hologramStack)[slice + 3];

      // vImgPhase1 -> update_pointer_single( (*i_V_hologramStack)[slice + 0] );        // ne met pas à jour phase1!!
      // extract_subImage(i_images_array[wrt_ptr + 0], phase1, p_image_dim_x, p_image_dim_y, p_window_edge_x, p_window_edge_y, p_window_size, p_window_size);
    }

  // important: current_angle est donné par un vecteur d'angles vérifiés corrects, donc pas nécessairement consécutifs
  // en version disque, c'est robuste car on lit à partir des fichiers testés
  // en version mémoire, ça fonctionne mais fragile. En effet, current_angle a été remplir pour décrire les numéros d'angles de 0 à n-1 sans interruption.
  
  //charge_decoupe_image_cv(phase1, holo_filename, tampon_image_cv, p_window_edge_x, p_window_edge_y, p_window_size, p_window_size);

  int fuse = 0;

  // ----------------------------------------
  // cas standard et rapide
  if (! p_modref_correction_mode)
    {

      // 0.29s cumulé. cast?
#pragma omp parallel for num_threads( 4 ) 
      for(size_t pixel = 0; pixel < c_N; pixel++)
	{
	  double mask = masque_hamming[pixel];
	  plan_reel_out[pixel] = (double)((phase1[pixel] - phase3[pixel]) * mask); 
	  plan_imag_out[pixel] = (double)((phase4[pixel] - phase2[pixel]) * mask); 
	}
#pragma omp barrier
      fuse++;
    }



  // ----------------------------------------
  // correction mod_ref et normalisation préalable
  if (p_modref_correction_mode)
    {    
      normalize_vectra1(phase1, phase1n, c_N);
      normalize_vectra1(phase2, phase2n, c_N);
      normalize_vectra1(phase3, phase3n, c_N);
      normalize_vectra1(phase4, phase4n, c_N);
  

#pragma omp parallel for num_threads( 4 ) 
      for(size_t pixel = 0; pixel < c_N; pixel++)
	{
	  
	  double p1 = phase1n[pixel];
	  double p2 = phase2n[pixel];
	  double p3 = phase3n[pixel];
	  double p4 = phase4n[pixel];

	  double p1mp3 = p1 - p3;
	  double p4mp2 = p4 - p2;

	  
	  double ampli_ref =  4.0f * sqrt(mod_refn[pixel]); //selon la formule
	  double mask = masque_hamming[pixel] / ampli_ref;
	  plan_reel_out[pixel] = (p1mp3 * mask); 
	  plan_imag_out[pixel] = (p4mp2 * mask); 
	  
      	}
#pragma omp barrier
      fuse++; 
    }

  ASSERT(fuse == 1);
}


// =============================================================================
// =============================================================================
// pas de phase shifting du tout (hors-axe pur)


void
Holo_Process::no_phase_shifting(size_t current_angle, double* plan_reel_out, double* plan_imag_out)
{

  // Lecture des interférogrammes via fichiers
  if (! p_read_from_mem_mode )
    { 
      compute_hologram_filename(p_images_radix, holo_filename, current_angle, 1);
      if (i_holo_resize_needed_p)
	{
	  i_vImgHoloReader -> reload_pgm(holo_filename);
	  vImgPhase1 -> copy_from_bigger(*i_vImgHoloReader, p_window_edge_x, p_window_edge_y);
	}
      else
	vImgPhase1 -> reload_pgm(holo_filename);
    }
  else // lecture des interférogrammes via mémoire 
    {
      phase1 = (*i_V_hologramStack)[ current_angle - 1 ];
    }
  
  int fuse = 0;


  // ----------------------------------------
  // cas standard et rapide
  if (! p_modref_correction_mode)
    {

      // 0.29s cumulé. cast?
#pragma omp parallel for num_threads( 4 ) 
      for(size_t pixel = 0; pixel < c_N; pixel++)
	{
	  plan_reel_out[pixel] = (double)(phase1[pixel] * masque_hamming[pixel]);
	  //plan_imag_out[pixel] = 0.0; // déjà fait dans import_R_from
	}
#pragma omp barrier
      fuse++;
    }



  // ----------------------------------------
  // cas alternatif: correction mod_ref et normalisation préalable
  if (p_modref_correction_mode)
    {    
      normalize_vectra1(phase1, phase1n, c_N);  

#pragma omp parallel for num_threads( 4 ) 
      for(size_t pixel = 0; pixel < c_N; pixel++)
	{	  
	  double p1 = phase1n[pixel];

	  double ampli_ref =  4.0f * sqrt(mod_refn[pixel]); //selon la formule
	  double mask = masque_hamming[pixel] / ampli_ref;
	  plan_reel_out[pixel] = (p1 * mask); 
	  plan_imag_out[pixel] = 0.0;
	  
      	}
#pragma omp barrier
      fuse++; 
    }

  ASSERT(fuse == 1);
}


// =============================================================================
// =============================================================================
// crop


// fonction permettant l'extraction de sous-images centrées ou non
// edge_x/y: désigne la position du coin supérieur gauche de la fenêtre de découpe 
// (dst_dim_x/y) dans l'image source (src_dimx/y)
template <typename T>
void
Holo_Process::extract_subImage(T* src, T* dst, size_t src_dimx, size_t src_dimy, size_t edge_x, size_t edge_y, size_t dst_dimx, size_t dst_dimy)
{
  MSG_ASSERT( (dst_dimx + edge_x < src_dimx) && \
	      (dst_dimy + edge_y < src_dimy),	\
	      "cropped window too far (x)");

  // on n'a besoin que du début de la ligne, on triera ensuite
  T* readPos = src + edge_y * src_dimx;
  T* writePos = dst;

  // taille d'un transfert de ligne
  size_t nbytes = dst_dimx * sizeof(T);

  // on effectue une série de transferts ligne par ligne
  // on copie à chaque itération la partie utile de la ligne
  for (size_t l = 0; l < dst_dimy; l++)
    {
      memcpy(writePos, readPos + edge_x, nbytes);
      readPos += src_dimx;
      writePos += dst_dimx;
    }
}


// découpage centré d'une image dst dans src
void
Holo_Process::square_crop(double *src, double *dst, size_t src_size, size_t dst_size)
{
  size_t edge = (src_size - dst_size) / 2;
  extract_subImage<double>(src, dst, src_size, src_size, edge, edge, dst_size, dst_size);
}


// découpage non centré du carré englobant un disque dont on précise la position du centre et le rayon
void
Holo_Process::disk_crop(double *src, double *dst, size_t src_size, size_t c_x, size_t c_y, size_t rayon) 
{
  size_t edge_x = c_x - rayon;
  size_t edge_y = c_y - rayon;
  size_t dst_size = 2 * rayon;
  extract_subImage<double>(src, dst, src_size, src_size, edge_x, edge_y, dst_size, dst_size);
}


void 
Holo_Process::compute_hologram_filename(const char* images_radix, char* filename, size_t cpt_angle, size_t cpt_dp)
{
  ASSERT((cpt_dp >= 1) && (cpt_dp <= 4));
  if (! i_proper_names_p)
    sprintf(filename, "%s%d-%03d." INPUT_FORMAT, images_radix, cpt_angle, cpt_dp);
  else
    sprintf(filename, "%s%03d-%03d." INPUT_FORMAT, images_radix, cpt_angle, cpt_dp);
}











/* utile en réflexion

// ----------------------------------------
// cas standard mais correction module
if ((! p_modref_correction_mode) && p_module_correction_mode)
{

#pragma omp parallel for num_threads( 4 ) 
for(size_t pixel = 0; pixel < c_N; pixel++)
{
double p1 = phase1[pixel];
double p2 = phase2[pixel];
double p3 = phase3[pixel];
double p4 = phase4[pixel];

double p1mp3 = p1 - p3;
double p4mp2 = p4 - p2;

double Module = 2 * sqrt( sqr(p1mp3) + sqr(p4mp2) ) / (p1 + p2 + p3 + p4);

if (Module > 0.1)
{
// attention, ajouter mod_ref
double mask = masque_hamming[pixel];
plan_reel_out[pixel] = (p1mp3 * mask); 
plan_imag_out[pixel] = (p4mp2 * mask); 
}
else
{ plan_reel_out[pixel] = plan_imag_out[pixel] = 0; }

}
#pragma omp barrier
fuse++;
}
*/
