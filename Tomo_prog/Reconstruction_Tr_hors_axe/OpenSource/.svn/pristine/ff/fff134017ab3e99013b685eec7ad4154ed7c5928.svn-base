#ifndef _TOOLS_
#define _TOOLS_

/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2001 the GLOS development team (http://glos.sourceforge.net) */


/* set TOOLS_T */
#include "./Tools_types.h"


/* round given value to nearest integer (floor if < floor +0.5, ceil
   otherwise) */
int 
i_round(double val);

/** prefixe i_ rajouté pour cause de MAJ mathcalls.h */

int 
i_round_f(float val);


/* return position of maximum element of given array (size > 0) */
int 
max_pos(TOOLS_T *array, int n);


/* */
TOOLS_T 
array_sum(TOOLS_T *tab, int size); 


/* return position of value in given array (size > 0), -1 if not found */
int 
in_array(TOOLS_T *array, int n, TOOLS_T value);


#endif /* TOOLS */
