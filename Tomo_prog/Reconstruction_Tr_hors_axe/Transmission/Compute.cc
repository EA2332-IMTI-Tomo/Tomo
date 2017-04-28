// ultra-déconne à mort si mis après les classes templates à moi (thread surtout)
#include "vChrono.h" // requiert boost
#include "vChronos.h" 

#include <iostream>
#include <boost/thread.hpp>  
#include <boost/date_time.hpp>


using namespace std;

#include "cv.h"
#include "highgui.h"

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
using namespace cv;



#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>
#include <exception>


#include "recon_includes.h"
#include "util.h"



#include <AIR_Volume.h>
#include <COMPLEX_Volume.h>
#include <FFTW_Volume.h>


#include "FFTW_Image.hpp"

#include "Holo_Process.h"
#include "Compute.h"
#include "cvDisplayVolume.hpp"






//******************************************************************************
// 
//******************************************************************************


extern char *g_OUTPUT_DIR;
extern char *g_OUTPUT_RADIX;
// extern size_t g_fftw_threads;


#define IMG_TYPE_FILE unsigned char

//==============================================================================
// Thread de visualisation de volume
//==============================================================================


void 
BoostThreadFunc(cvDisplayVolume<RECON_TYPE> *slicer)
{
  do
    {
      slicer -> updateImage();
    }
  while( slicer -> showImage(30) );
    

  boost::posix_time::seconds workTime(3);
  // Pretend to do something useful...
  boost::this_thread::sleep(workTime);
         
  cout << "*****************************************************************************" << endl;
  cout << "cvDisplayVolume Worker: finished" << endl;
}


//==============================================================================
// version offline 
//==============================================================================


#include "compute1.cc"


//******************************************************************************
// version offline progressive 
//******************************************************************************


#include "compute2.cc"



//******************************************************************************
// version online progressive
//******************************************************************************


#include "compute3.cc"



//******************************************************************************
// version film offline
//******************************************************************************


#include "compute4.cc"


//******************************************************************************
// version spécifique au test interne de la lecture mémoire
//******************************************************************************


#include "compute6.cc"






//******************************************************************************
// version zoom x2 offline
//******************************************************************************


//#include "compute5.cc"
void
compute_cpu_zoom_offline(size_t cube_edge, size_t Nxmax, size_t Nymax, size_t Nxmax_Rf, \
			 size_t window_edge_x, size_t window_edge_y, size_t image_dim_x, size_t image_dim_y, \
			 size_t xm0_limite, size_t ym0_limite, float rayon, float delta_zmax, \
			 size_t  angle_start, size_t angle_count, size_t angle_jump, \
			 bool off_axis_and_psh, size_t offaxis_cx, size_t offaxis_cy, \
			 bool autofocus_gradient,			\
			 const string &images_radix, bool proper_numbers, bool modref_asked, const char* modref_file)
{
  ASSERT(false);
}



