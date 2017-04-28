    
  // ==========================================================
  // Exécution de FFT 3D précédée et suivie de shifts circulaires
  // ==========================================================

  
  // à ce stade, le volume ne fait plus que 508^3, soit 4*127 de large. 
  // Il y a une ligne de rognée en 128*128 qui se répercute


  // *******************************************************
  // circshift avant
  

  double *reel_arc_shift, *imag_arc_shift;
  TIMES("circshift avant TF3D\n", {
      
      /////////////////////creation des tableaux qui recoivent le circshift
      //partie reelle du papillon dans l'espace reciproque  //
      ARRAY_NEW(reel_arc_shift, N_tab, double);
            
      circshift3D_memcpy(reel_arc, reel_arc_shift, dv0rf, dv0rf, dv0rf);
      delete reel_arc;
  
      // partie imaginaire du papillon dans l'espace reciproque //
      ARRAY_NEW(imag_arc_shift, N_tab, double);

      circshift3D_memcpy(imag_arc, imag_arc_shift, dv0rf, dv0rf, dv0rf);
      delete imag_arc;

    });

 
  
  // *******************************************************
  // TF3D
  
  printf("*******************************************\n");
  printf("TF3D\n");
  printf("*******************************************\n");
  fingers();


  // à ce stade, le volume ne fait plus que 508^3, soit 4*127 de large. 
  // Il y a une ligne de rognée en 128*128 qui se répercute
  
  int N3D = N_tab;
  
  TIMES("TF3D", {    
      
#ifdef CUDA
      if (g_fft3d_cuda)
	TF3D_cuda(4 * Nxmax_Rf, 4 * Nxmax_Rf, 4 * Nxmax_Rf, reel_arc_shift, imag_arc_shift);
      else
#endif
	TF3D(N3D, 4 * Nxmax_Rf, reel_arc_shift, imag_arc_shift);      
    });

  double *final_reel = reel_arc_shift;
  double *final_imag = imag_arc_shift;


  // *******************************************************
  // circshift apres TF3D

  
  double *final_reel_shift, *final_imag_shift;
  //partie reelle du resultat final shifté (centré)
  ARRAY_NEW(final_reel_shift, N_tab, double);
  //partie imaginaire du resultat final shifté (centré)
  ARRAY_NEW(final_imag_shift, N_tab, double);
  
  
  TIMES("Circshift retour",
	{
	  circshift3D_memcpy(final_reel, final_reel_shift, dv0rf, dv0rf, dv0rf);
	  circshift3D_memcpy(final_imag, final_imag_shift, dv0rf, dv0rf, dv0rf);
	  delete final_reel;
	  delete final_imag;
	});





  // ===========================================================================
  // ===========================================================================
  // ===========================================================================


  
