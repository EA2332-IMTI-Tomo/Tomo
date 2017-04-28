// Ce fichier sert de patron pour un appel externe du programme
#include "cv.h"
#include "cpu_launcher.h"
#include "util_Image.h"

using namespace std;

bool PROPER_NAMES = true;
void 
dat_compute_hologram_filename(const char* images_radix, char* filename, size_t cpt_angle, size_t cpt_dp);


// ============================================================================
// ============================================================================

int 
main(int argc, char** argv)
{
  MSG_ASSERT((argc - 1 == 2), "fichier_cfg, rep_images_hdp");
  cout << endl << argc;
  cout << endl << argv[1];
  cout << endl << argv[2];

  unsigned char** interferograms;


  // ======================================================================
  // ce que manipTomo va produire

  // obligé de faire ça!
  size_t IMAGE_DIMX = 740;
  size_t IMAGE_DIMY = 572;
  size_t FINAL_ANGLE = 511;
  string images_radix = argv[2]; images_radix += "/session10032302-record";

  char* holo_filename;
  ARRAY_ALLOC(holo_filename, 1000, char);
  ARRAY_ALLOC(interferograms, FINAL_ANGLE * 4, unsigned char*);
  for (size_t j = 0; j < FINAL_ANGLE * 4; j++)
    ARRAY_ALLOC(interferograms[j], IMAGE_DIMX * IMAGE_DIMY, unsigned char);

  cout << endl << endl << "nombre d'angles considérés corrects" << FINAL_ANGLE << endl;
  for (size_t cursor = 0; cursor < FINAL_ANGLE; cursor++)
    {
      size_t current_angle = cursor + 1;
      cout << " | " << current_angle;
      
      for (size_t k = 0; k < 4; k++)
	{
	  dat_compute_hologram_filename(images_radix.c_str(), holo_filename, current_angle, k+1);
	  cout << endl << "nb angles eff:" << FINAL_ANGLE << "holo_filename:" << holo_filename;
	  remplit_tableau_cv(interferograms[cursor * 4 + k], holo_filename, IMAGE_DIMX, IMAGE_DIMY);
	}
    }

  // ======================================================================

  compute_cpu(argv[1], interferograms);
}











void 
dat_compute_hologram_filename(const char* images_radix, char* filename, size_t cpt_angle, size_t cpt_dp)
{
  ASSERT((cpt_dp >= 1) && (cpt_dp <= 4));
  if (! PROPER_NAMES)
    sprintf(filename, "%s%d-%03d." INPUT_FORMAT, images_radix, cpt_angle, cpt_dp);
  else
    sprintf(filename, "%s%03d-%03d." INPUT_FORMAT, images_radix, cpt_angle, cpt_dp);
}



