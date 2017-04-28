/* GLOS (Graphic Library in Open Source), an ANSI Common Lisp OpenGL subset.
   Copyright (C) 2000 the GLOS development team (http://glos.sourceforge.net) */


#ifndef _VECTRA_MACROS_
#define _VECTRA_MACROS_


#include <time.h>
#include <assert.h>
#include <stdlib.h>


/** warning: numerical expressions are subject to syntax errors or
    operators precedence problems, unless the parameters are quoted in
    the body inside parenthesis */


/*******************************************************************************
 Booleans
*******************************************************************************/


#define BOOL unsigned short int
#define FALSE 0
#define TRUE 1



/*******************************************************************************
 Iterators
*******************************************************************************/


/* sais plus a quoi, mais ca sert */
#define WITH_BREAK(CODE_EXEC) \
{ \
  int I; \
  for (I = 0; I < 1; I++) \
    { \
	CODE_EXEC; \
    } \
}


/*******************************************************************************
 Connexity (img.h)
*******************************************************************************/


/* calls given function of arity 2 on every pixel coordinates around I, J 
   in 2D 4 connexity */
#define AROUND_4_CONNEX(I, J, FUN_OF_IJ) \
FUN_OF_IJ(I + 1, J); \
FUN_OF_IJ(I, J + 1); \
FUN_OF_IJ(I - 1, J); \
FUN_OF_IJ(I, J - 1);

/* 8 connexity without 4 connexity cases */
#define AROUND_8_ONLY_CONNEX(I, J, FUN_OF_IJ) \
FUN_OF_IJ(I + 1, J + 1); \
FUN_OF_IJ(I - 1, J + 1); \
FUN_OF_IJ(I - 1, J - 1); \
FUN_OF_IJ(I + 1, J - 1); 

/* calls given function of arity 2 on every pixel coordinates around I, J 
   in 2D 8 connexity. Order comply with Freeman code run. */
#define AROUND_8_FREEMAN_CONNEX(I, J, FUN_OF_IJ) \
FUN_OF_IJ(I + 1, J); \
FUN_OF_IJ(I + 1, J + 1); \
FUN_OF_IJ(I, J + 1); \
FUN_OF_IJ(I - 1, J + 1); \
FUN_OF_IJ(I - 1, J); \
FUN_OF_IJ(I - 1, J - 1); \
FUN_OF_IJ(I, J - 1); \
FUN_OF_IJ(I + 1, J - 1); 

/* */
#define AROUND_8_CONNEX(I, J, FUN_OF_IJ) \
AROUND_4_CONNEX(I, J, FUN_OF_IJ) \
AROUND_8_ONLY_CONNEX(I, J, FUN_OF_IJ)


/* the same , but if local boolean WITH_HALT is set to
   TRUE, iteration stops. */
#define AROUND_COND_8_CONNEX(I, J, FUN_OF_IJ) \
WITH_BREAK({ \
{ \
  BOOL WITH_HALT = FALSE; \
  FUN_OF_IJ(I + 1, J);     if (WITH_HALT) break; \
  FUN_OF_IJ(I, J + 1);     if (WITH_HALT) break; \
  FUN_OF_IJ(I - 1, J);     if (WITH_HALT) break; \
  FUN_OF_IJ(I, J - 1);     if (WITH_HALT) break; \
  FUN_OF_IJ(I + 1, J + 1); if (WITH_HALT) break; \
  FUN_OF_IJ(I - 1, J + 1); if (WITH_HALT) break; \
  FUN_OF_IJ(I - 1, J - 1); if (WITH_HALT) break; \
  FUN_OF_IJ(I + 1, J - 1); if (WITH_HALT) break; \
}});


/* the same in Freeman order, but if local boolean WITH_HALT is set to
   TRUE, iteration stops. */
#define AROUND_COND_8_FREEMAN_CONNEX(I, J, FUN_OF_IJ) \
WITH_BREAK({ \
{ \
  BOOL WITH_HALT = FALSE; \
  FUN_OF_IJ(I + 1, J);     if (WITH_HALT) break; \
  FUN_OF_IJ(I + 1, J + 1); if (WITH_HALT) break; \
  FUN_OF_IJ(I, J + 1);     if (WITH_HALT) break; \
  FUN_OF_IJ(I - 1, J + 1); if (WITH_HALT) break; \
  FUN_OF_IJ(I - 1, J);     if (WITH_HALT) break; \
  FUN_OF_IJ(I - 1, J - 1); if (WITH_HALT) break; \
  FUN_OF_IJ(I, J - 1);     if (WITH_HALT) break; \
  FUN_OF_IJ(I + 1, J - 1); if (WITH_HALT) break; \
}});


/* calls given function of arity 3 on every voxel coordinates around I, J, K 
   in 3D 6 connexity */
#define AROUND_6_CONNEX(I, J, K, FUN_OF_IJK) \
FUN_OF_IJK(I + 1, J, K); \
FUN_OF_IJK(I - 1, J, K); \
FUN_OF_IJK(I, J + 1, K); \
FUN_OF_IJK(I, J - 1, K); \
FUN_OF_IJK(I, J, K + 1); \
FUN_OF_IJK(I, J, K - 1);


/* calls given function of arity 3 on every voxel coordinates around I, J, K 
   in 3D 6 connexity */
#define AROUND_26_CONNEX(I, J, K, FUN_OF_IJK) \
/* first K slide around voxel */ \
FUN_OF_IJK(I + 1, J, K); \
FUN_OF_IJK(I - 1, J, K); \
FUN_OF_IJK(I + 1, J + 1, K); \
FUN_OF_IJK(I + 1, J - 1, K); \
FUN_OF_IJK(I - 1, J + 1, K); \
FUN_OF_IJK(I - 1, J - 1, K); \
FUN_OF_IJK(I, J + 1, K); \
FUN_OF_IJK(I, J - 1, K); \
/* the same at K+1 (one extra voxel) */ \
FUN_OF_IJK(I + 1, J, K + 1); \
FUN_OF_IJK(I - 1, J, K + 1); \
FUN_OF_IJK(I + 1, J + 1, K + 1); \
FUN_OF_IJK(I + 1, J - 1, K + 1); \
FUN_OF_IJK(I - 1, J + 1, K + 1); \
FUN_OF_IJK(I - 1, J - 1, K + 1); \
FUN_OF_IJK(I, J + 1, K + 1); \
FUN_OF_IJK(I, J - 1, K + 1); \
FUN_OF_IJK(I, J, K + 1); \
/* the same at K-1 */ \
FUN_OF_IJK(I + 1, J, K - 1); \
FUN_OF_IJK(I - 1, J, K - 1); \
FUN_OF_IJK(I + 1, J + 1, K - 1); \
FUN_OF_IJK(I + 1, J - 1, K - 1); \
FUN_OF_IJK(I - 1, J + 1, K - 1); \
FUN_OF_IJK(I - 1, J - 1, K - 1); \
FUN_OF_IJK(I, J + 1, K - 1); \
FUN_OF_IJK(I, J - 1, K - 1); \
FUN_OF_IJK(I, J, K - 1); 


/*******************************************************************************
 Generic 
*******************************************************************************/


/* execute given code unless condition is true */
#define UNLESS(COND, CODE_EXEC) \
if (! (COND)) \
  CODE_EXEC
/* usage:
UNLESS(t == -1),
{
  blah
})*/


/* reverse the values of A and B 
   requires extra variable of same type for value swapping */
#define REVERSE_VALUES(A, B, SWP) \
SWP = B; B = A; A = SWP;


/* display an integer in a lazy way */
#define DISP_INT(NUM) \
printf("\n%d", NUM);


/** times the execution of given code. Displays given message at start
    and end of execution, mentionning elapsed time in the former
    occurence */
#define TIME(MSG, CODE_EXEC) \
{ \
  time_t t_start; \
  time_t t_end; \
  double seconds_elapsed; \
  \
  t_start = time(NULL); \
  printf("\n" MSG ": START "); \
  \
  CODE_EXEC; \
  \
  t_end = time(NULL); \
  seconds_elapsed = difftime(t_end, t_start); \
  printf("\n" MSG ": COMPLETE in %5.1f seconds \n", seconds_elapsed); \
}
/* usage:
TIME("my message",
{
  code to time;
});
*/ 


#define DEFUN_MAX_FUNCTION(TYPE) \
TYPE \
max_##TYPE (TYPE* array, int start, int card) \
{ \
  int i, end = start + card; \
  TYPE max = array[start++]; \
  \
  for (i = start; i < end; i++) \
    { \
      if (array[i] > max) max = array[i]; \
    } \
  \
  return max; \
}
/* declares:
   TYPE max_TYPE(TYPE* array, int start, int card);
   
   which returns maximum value in array from a range of card elements
   starting from given index */


/*******************************************************************************
 Math
*******************************************************************************/


/* REM: passer en inline */
/* tests if DOW <= X <= UP, since this expression is not evaluated by
   egcs as we could expect */
#define BETWEEN(X, DOW, UP) \
(((DOW) <= (X)) && ((X) <= (UP)))


/* tests if X is between range [mean - sd, mean + sd] */
#define BETWEEN_SD(X, MEAN, SD) \
BETWEEN((X), ((MEAN) - (SD)), ((MEAN) + (SD)))


/* */
#define add_limited(NUM, ADD, LIM) \
  ((NUM) + (ADD) > (LIM)) ? (NUM) = (LIM) : (NUM) += (ADD); 


/* */
#define sub_limited(NUM, SUB, LIM) \
  ((NUM) - (SUB) <= (LIM)) ? (NUM) = (LIM) : (NUM) -= (SUB); 


/* */
#define sub_positive(NUM, SUB) \
  sub_limited(NUM, SUB, 0)


/* determines the square of a number (genericity needed) */
#define SQR(EXP) \
  ((EXP) * (EXP))
/* NB: parenthesis enable use of a numerical expression without
   precedence operator problems */


/** notice expressions are of form ( (EXP1) )*/


/* determines the absolute value of a number */
/* works for expressions but with multiple evals */
#define ABS(EXP) \
  (((EXP) > 0) ? (EXP) : -(EXP))


/* return the minimum value of both expressions values*/
#define MIN(EXP1, EXP2) \
  (((EXP1) > (EXP2)) ? (EXP2) : (EXP1))


/* return the maximum value of both expressions values*/
#define MAX(EXP1, EXP2) \
  (((EXP1) > (EXP2)) ? (EXP1) : (EXP2))

/*
#if defined(_MSC_VER) 
#  define MIN(a,b) std::_cpp_min(a,b)
#  define MAX(a,b) std::_cpp_max(a,b)
#else
#  define MIN(a,b) std::min(a,b)
#  define MAX(a,b) std::max(a,b)
#endif
*/

/* ensure VAR stays within DOWN =< VAR =< UP, set to boundary value
   otherwise */
#define FRAME_VALUE(VAR, DOWN, UP) \
if (VAR < DOWN) \
  VAR = DOWN; \
else if (VAR > UP) \
  VAR = UP;


/* */
#define BOUNDEDP(VAR, DOWN, UP) \
((VAR > DOWN) && (VAR < UP))

#define BOUNDED_EQ_P(VAR, DOWN, UP) \
((VAR >= DOWN) && (VAR <= UP))


/*******************************************************************************
 Memory-allocation related
*******************************************************************************/


#define ARRAY_ALLOC(TAB, DIM, DATATYPE) \
assert(TAB = (DATATYPE *) calloc(DIM, sizeof(DATATYPE))); 

#define ARRAY_SET(TAB, DIM, VAL) \
for (int i = 0; i < DIM; i++) \
{ \
   TAB[DIM] = VAL; \
}

#define ARRAY_CPY(TAB_SRC, TAB_DST, I, J) \
for (int n = I; n <= J; n++) \
{ \
   TAB_DST[n] = TAB_SRC[n]; \
}



#define ARRAY_REALLOC(TAB, NEWDIM, DATATYPE) \
assert(TAB = (DATATYPE *) realloc(TAB, NEWDIM * sizeof(DATATYPE))); 


/* allocates to given MAT pointer enough memory for a matrix of 
   X * Y elements of type DATATYPE. Initialization to 0.
   - Requires <assert.h> */
#define MATRIX_ALLOC(MAT, X, Y, DATATYPE) \
assert(MAT = (DATATYPE **) malloc(X * sizeof(DATATYPE*))); \
{ \
  int CTR; \
  for (CTR = 0; CTR < X; CTR++) \
    {assert(MAT[CTR] = (DATATYPE *) calloc(Y, sizeof(DATATYPE)));} \
} 
/* modified from img.h: 
   calloc used: elements initialized to 0
   CTR declared on heap and not inherited 
   X and Y reverted: now mat[X-1][Y-1] exists */


/* frees memory allocated by MATRIX_ALLOC */
#define MATRIX_FREE(MAT, X) \
{ \
  int I; \
  for (I = 0; I < X; I++) \
    free(MAT[I]); \
  free(MAT); \
} 
/* this could be a (inline?) function since type-independent*/


/* allocates to given VOL pointer enough memory for a volume of 
   X * Y * Z elements of type DATATYPE. Initialization to 0.
   - Requires <assert.h> */
#define VOLUME_ALLOC(VOL, X, Y, Z, DATATYPE) \
assert(VOL = (DATATYPE ***) malloc(X * sizeof(DATATYPE**))); \
{ \
  int I, J; \
  for (I = 0; I < X; I++) \
    { \
      assert(VOL[I] = (DATATYPE **) calloc(Y, sizeof(DATATYPE*))); \
      for (J = 0; J < Y; J++) \
        assert(VOL[I][J] = (DATATYPE *) calloc(Z, sizeof(DATATYPE))); \
    } \
} 
/* to declare a volume in SEG convention:
   VOLUME_ALLOC(vol, dimz, dimx, dimy);
   vol[z][x][y]; */


/** a verif sur poly michel, je crois que type facultatif */
#define MATRIX_ALLOC_2(MYMAT, X, Y, DATATYPE) \
{ \
  int I; \
  DATATYPE **MAT; \
  \
  assert(MAT = (int **) calloc(X, sizeof(*MAT))); \
  assert(*MAT=(int *) calloc (X * Y, sizeof(**MAT))); \
  for (I = 1; I < X; I++) \
    MAT[I] = MAT[I-1] + Y; \
  \
  MYMAT = **MAT; \
}


/* frees memory allocated by VOLUME_ALLOC */
#define VOLUME_FREE(VOL, X, Y) \
{ \
  int I, J; \
  for (I = 0; I < X; I++) \
    { \
      for (J = 0; J < Y; J++) \
        free(VOL[I][J]); \
    } \
} 
/* this could be a (inline?) function since type-independent*/


/* allocates to given MAT pointer enough memory for a matrix of 
   y * x elements of type DATATYPE.
   - Requires extra integer variable INT_CTR.
   - Requires <assert.h> */
#define REV_MATRIX_ALLOC(MAT, X, Y, DATATYPE, INT_CTR) \
assert(MAT = (DATATYPE **) malloc(Y * sizeof(DATATYPE*))); \
for (INT_CTR = 0; INT_CTR < Y; INT_CTR++) \
    {assert(MAT[INT_CTR] = (DATATYPE *) malloc(X * sizeof(DATATYPE)));} 
/* Attention: colle probablement pas avec SEG */


/*******************************************************************************
 n-arrays operations										
*******************************************************************************/


/* Reverses (destructively) the N first elements of given array */
#define ARRAY_REVERSE(TAB, N, SWP) \
{ \
  int I; \
  int STOP = (int)N / 2; \
  for (I = 0; I < STOP; I++) \
    {REVERSE_VALUES(TAB[I], TAB[N - 1 - I], SWP);} \
}


/* initializes matrix to given value */
#define MATRIX_INIT(MAT, X, Y, VAL) \
{ \
  int T_I, T_J; \
  for (T_I = 0; T_I < X; T_I++) \
    for (T_J = 0; T_J < Y; T_J++) \
      MAT[T_I][T_J] = VAL; \
}

/* copy matrix MAT1 to MAT2 */
#define MATRIX_COPY(MAT1, MAT2, X, Y) \
{ \
  int T_I, T_J; \
  for (T_I = 0; T_I < X; T_I++) \
    for (T_J = 0; T_J < Y; T_J++) \
      MAT2[T_I][T_J] = MAT1[T_I][T_J]; \
}


/* initializes matrix to given file content. TYPECODE = "f" for float,
   "d" for integer, "lf" for double, etc */
#define MATRIX_FILE_TXT(MAT, X, Y, FILENAME, TYPECODE) \
{ \
  FILE* fp; \
  int T_I, T_J; \
 \
  assert(fp = fopen(FILENAME, "r")); \
 \
  for (T_I = 0; T_I < X; T_I++) \
    for (T_J = 0; T_J < Y; T_J++) \
      assert(fscanf(fp, "%" TYPECODE, &MAT[T_I][T_J]) != EOF); \
 \
 fclose(fp); \
}


/* initializes matrix to given file content. TYPECODE = "f" for float,
   "d" for integer, etc */
#define MATRIX_FILE_RAW(MAT, X, Y, FILENAME, D_TYPE) \
{ \
  FILE* fp; \
  int T_I, T_J; \
 \
  assert(fp = fopen(FILENAME, "rb")); \
 \
  for (T_I = 0; T_I < X; T_I++) \
    for (T_J = 0; T_J < Y; T_J++) \
      assert(fread(&MAT[T_I][T_J], sizeof(D_TYPE), 1, fp) != EOF); \
 \
 fclose(fp); \
}


#define FILE_MATRIX_TXT(MAT, X, Y, FILENAME, TYPECODE) \
{ \
  FILE* fp; \
  int T_I, T_J; \
 \
  assert(fp = fopen(FILENAME, "w")); \
 \
  for (T_I = 0; T_I < X; T_I++) \
  { \
    for (T_J = 0; T_J < Y; T_J++) \
      fprintf(fp, "%" TYPECODE " ", MAT[T_I][T_J]); \
    fprintf(fp, "\n"); \
  } \
 \
 fclose(fp); \
}



#define MATRIX_DISPLAY(MAT, X, Y, TYPECODE) \
{ \
  int T_I, T_J; \
  \
  for (T_I = 0; T_I < X; T_I++) \
    { \
      printf("\n"); \
      for (T_J = 0; T_J < Y; T_J++) \
        printf("%" TYPECODE " ", MAT[T_I][T_J]); \
    } \
}


/* initializes volume to given value */
#define VOLUME_INIT(MAT, X, Y, Z, VAL) \
{ \
  int T_I, T_J, T_K; \
  for (T_I = 0; T_I < X; T_I++) \
    for (T_J = 0; T_J < Y; T_J++) \
      for (T_K = 0; T_K < Z; T_K++) \
        MAT[T_I][T_J][T_K] = VAL; \
}


/* transfer a whole slice (possibly of a volume, ie calling vol[k]) */
#define TRANSFER_SLICE(slice_in, slice_out, x, y) \
{ \
  int a, b; \
  for (a = 0; a < x; a++) \
    for (b = 0; b < y; b++) \
      slice_out[a][b] = slice_in[a][b]; \
}


/*******************************************************************************
 Class Helper Macros
*******************************************************************************/


/* automatically define inline accessors in class definition */

#define READER(TYPE, VARNAME) \
TYPE get_##VARNAME() const \
{ return this -> VARNAME; };

/* when the accessed variable is const */
#define CREADER(TYPE, VARNAME) \
const TYPE get_##VARNAME() const \
{ return this -> VARNAME; };

#define WRITER(TYPE, VARNAME) \
void set_##VARNAME(const TYPE &value) \
{ this -> VARNAME = value; };

#define READWRITER(TYPE, VARNAME) \
READER(TYPE, VARNAME) \
WRITER(TYPE, VARNAME)


#endif /* _VECTRA_MACROS_ */
