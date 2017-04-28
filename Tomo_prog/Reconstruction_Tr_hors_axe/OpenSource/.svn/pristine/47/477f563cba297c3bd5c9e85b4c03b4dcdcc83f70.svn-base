#include "vectra.h"

#include "assert.h"
#include <iostream>


using namespace std;


int 
main(int argc, char** argv)
{
  cout << endl << "arg1: fichier: " << argv[1];
  cout << endl << "il existe et est lisible: " << vectra::file_exists_p(argv[1]);

  cout << endl << "arg2: répertoire: " << argv[2];
  bool dire = vectra::dir_exists_p(argv[2]);
  cout << endl << "il existe?: " << dire;

  if (! dire)
    {
      cout << endl << "puisque non, on le crée";
      assert(vectra::makedir(argv[2], 0777));
    }


  cout << endl << "press key to exit";
  while (! vectra::kbhit());

  return 0;
}



