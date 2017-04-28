#ifndef __PARSER__
#define __PARSER__


#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <exception>
#include <map>

#include "recon_includes.h"

using namespace std;



void
init_vals(std::map<string, size_t> &map_dims, \ 
	  std::map<string, bool> &map_bools, \
	  std::map<string, string> &map_paths, \
	  std::map<string, float> &map_physics);



void
parse_file(const char* in_file, \ 
	   std::map<string, size_t> &map_dims, \ 
	   std::map<string, bool> &map_bools, \
	   std::map<string, string> &map_paths,
	   std::map<string, float> &map_physics);


void
parse_argts(size_t l_argc, char** l_argv, \ 
	    std::map<string, size_t> &map_dims, \ 
	    std::map<string, bool> &map_bools, \
	    std::map<string, string> &map_paths,
	    std::map<string, float> &map_physics,\
	    movieParams &movie_settings,\
	    enum ReconType &ReconMode);


void
compute_vals(std::map<string, size_t> &map_dims, \ 
	     std::map<string, float> &map_physics, \
	     std::map<string, float> &map_experiment);


void
check_vals(std::map<string, size_t> &map_dims, \
	   std::map<string, bool> &map_bools, \
	   std::map<string, string> &map_paths, \
	   std::map<string, float> &map_experiment);


#endif // __PARSER__
