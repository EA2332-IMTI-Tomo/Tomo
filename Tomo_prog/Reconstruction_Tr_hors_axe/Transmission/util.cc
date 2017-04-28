#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fftw3.h>
#include <cstring>
#include <fstream>


using namespace std;


//#include "main.h"
#include "macros.h"
#include "util.h"



void
birdy()
{
  printf("  \\,`//\n");
  printf(" _).. `_\n");
  printf("( __  -\\ \n");
  printf("    '`.\n");
  printf("   ( \\>\n");
  printf("   _||_ \n");
}


void
fingers()
{
  printf("     .\"\".    .\"\".\n");
  printf("     |  |   /  /\n");
  printf("     |  |  /  /\n");
  printf("     |  | /  /\n");
  printf("     |  |/  ;-.__ \n");
  printf("     |  ` _/  /  /\n");
  printf("     |  /` ) /  /\n");
  printf("     | /  /_/\\_/ \n");
  printf("     |/  /      |\n");
  printf("     (  ' \\ '-  |\n");
  printf("      \\    `.  /\n");
  printf("       |      |\n");
  printf("       |      |\n");
}



char*
str_alloc_cpy(char* str_src)
{
  char* str_dest;
  ASSERT(str_dest = (char *) calloc(1 + strlen(str_src), sizeof(char)));
  strcpy(str_dest, str_src);
  return str_dest;
}


// à virer à terme
char*
str_concat(const char* str1, const char* str2)
{
  char* str;
  int l1 = strlen(str1), l2 = strlen(str2);
  int i, j;

  ASSERT(str = (char*)calloc(1 + l1 + l2, sizeof(char))); 

  for (i = 0; i < l1; i++)
    str[i] = str1[i];

  for (j = 0; j < l2; j++)
    str[i + j] = str2[j];

  str[i + j] = '\0';

  return str;
}
  


/*
void
compute_hologram_filename(const char* images_radix, char* filename, size_t cpt_angle, size_t cpt_dp)
{
  //sprintf(filename, "%s%i-001." INPUT_FORMAT, path, i);
  ASSERT((cpt_dp >= 1) && (cpt_dp <= 4));
  sprintf(filename, "%s%d-%03d." INPUT_FORMAT, images_radix, cpt_angle, cpt_dp);
}
*/


/*
BOOL
ftest(char* filename)
{
  FILE* test_existence = fopen(filename, "rb");
  if (test_existence != NULL) {
    fclose(test_existence);
    return TRUE;
  }
  else {
    return FALSE;
  }
}




*/
