#include "vectra.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <sys/stat.h>
#include <termios.h>
#include <unistd.h>
#include <fcntl.h>

#include <errno.h>
#include <unistd.h>



// -----------------------------------------------------------------------------


bool
vectra::file_exists_p(const char* filename)
{
  FILE* fp;
  bool rval;
  
  fp = fopen(filename, "r");
  
  if (fp == NULL)
    return false;
  else
    { fclose(fp); return true; }
}


// -----------------------------------------------------------------------------


bool
vectra::dir_exists_p(char* filename)
{
  struct stat st;
  return ((stat(filename ,&st) == 0) );
}


// -----------------------------------------------------------------------------


bool
vectra::kbhit(void)
{
  struct termios oldt, newt;
  int ch;
  int oldf;

  tcgetattr(STDIN_FILENO, &oldt);
  newt = oldt;
  newt.c_lflag &= ~(ICANON | ECHO);
  tcsetattr(STDIN_FILENO, TCSANOW, &newt);
  oldf = fcntl(STDIN_FILENO, F_GETFL, 0);
  fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);

  ch = getchar();

  tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
  fcntl(STDIN_FILENO, F_SETFL, oldf);

  if(ch != EOF)
    {
      ungetc(ch, stdin);
      return true;
    }

  return false;
}

#define _getch getchar


// -----------------------------------------------------------------------------


// http://stackoverflow.com/questions/675039/how-can-i-create-directory-tree-in-c-linux


bool
vectra::makedir(const char *path, mode_t mode)
{
  typedef struct stat Stat;
  Stat            st;
  int             status = 0;

  if (stat(path, &st) != 0)
    {
      /* Directory does not exist. EEXIST for race condition */
      if (mkdir(path, mode) != 0 && errno != EEXIST)
	status = -1;
    }
  else if (!S_ISDIR(st.st_mode))
    {
      errno = ENOTDIR;
      status = -1;
    }

  return(status != -1);
}


// -----------------------------------------------------------------------------


unsigned int 
vectra::file_size(FILE* fp)
{
  unsigned int fd; // UNIX file descriptor
  struct stat l_filestat;
  
  if (fp == NULL) return 0;
  
  // retrieve size of the file
  fd = fileno(fp);
  fstat(fd, &l_filestat);
  return(l_filestat.st_size);
}

