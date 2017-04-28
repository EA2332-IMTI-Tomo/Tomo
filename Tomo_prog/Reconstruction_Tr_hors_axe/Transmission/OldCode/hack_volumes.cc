
cout << endl << "sauvegarde des volumes normalisés avant TF3D inverse et circshifts";
size_t l_matSize = cube_edge * cube_edge * cube_edge;

write_array<RECON_TYPE>(cube_edge * cube_edge * cube_edge, reel_arc, str_concat(g_OUTPUT_DIR, "/avantFFT_64_reel"));
write_array<RECON_TYPE>(cube_edge * cube_edge * cube_edge, imag_arc, str_concat(g_OUTPUT_DIR, "/avantFFT_64_imag"));

float * tab_arc_32;
ARRAY_NEW(tab_arc_32, l_matSize, float);
for (size_t i = 0; i < l_matSize; i++)
  tab_arc_32[i] = (float)reel_arc[i];
write_array<float>(l_matSize, tab_arc_32, str_concat(g_OUTPUT_DIR, "/avantFFT_32_reel"));
for (size_t i = 0; i < l_matSize; i++)
  tab_arc_32[i] = (float)imag_arc[i];
write_array<float>(l_matSize, tab_arc_32, str_concat(g_OUTPUT_DIR, "/avantFFT_32_imag"));
delete(tab_arc_32);
