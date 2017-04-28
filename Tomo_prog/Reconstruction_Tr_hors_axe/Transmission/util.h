#ifndef __RECON_UTIL__
#define __RECON_UTIL__


// en gros, tout ici est dispensable.

void
birdy();

void
fingers();


// returns the adress of an allocated copy of str_src
char*
str_alloc_cpy(char* str_src);


// deprecated, à virer à terme
char* 
str_concat(const char* str1, const char* str2); 



#endif // __RECON_UTIL__
