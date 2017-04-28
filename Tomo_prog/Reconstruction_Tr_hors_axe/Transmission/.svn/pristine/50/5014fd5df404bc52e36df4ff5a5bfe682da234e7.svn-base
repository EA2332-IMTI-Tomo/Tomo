#include <stdio.h>
#include <stlib>
#include <string.h>
#include <cstdlib>
#include <assert.h>

char*
str_concat(char* str1, char* str2)
{
  char* str;
  int msize = strlen(str1) + strlen(str2);

  assert(str = (char*)calloc(msize, sizeof(char))); 
  

  return str;
}


int
main(int argc, char** argv)
{
  assert(argc > 2); 

  fprintf(stdout, str_concat(argv[1], argv[2]));

  return EXIT_SUCCESS;
}


