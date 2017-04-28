#ifndef __CPU_LAUNCHER__
#define __CPU_LAUNCHER__



#include "vChrono.h" // requiert boost
#include "vChronos.h"


#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <exception>
#include <map>



#include "recon_includes.h"
#include "util.h"
#include "Compute.h"
#include "recon_parser.h"

#include "FFTW_Image.hpp"



// ============================================================================
// ============================================================================


void 
compute_cpu(const char* config_filename, unsigned char** interferograms);


#endif // __CPU_LAUNCHER__
