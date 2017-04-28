#ifndef _ROI_
#define _ROI_


/* #include "ROI_type.h" */
#include "macros.h"
#include "Point.h"


/*
  This module regroups every function related to regions of interest
  definition and handling.
  
  The ROI of a slice can be seen as a rectangle. Each slice has its
  own ROI. VOI conceptually consists in the number of the first and
  last slices of interests and the list of ROIs for each slice.

 */


/* -----------------------------------------------------------------------------
   
----------------------------------------------------------------------------- */


#define ROI_COLOR 128


/* -----------------------------------------------------------------------------
   Helper Macros
----------------------------------------------------------------------------- */


/* iterates execution of given code for all the values of local
   variables I and J, featuring points coordinates of ROI[][] */
#define RUN_THROUGH_ROI(K, CODE_EXEC) \
{ \
  int I, J; \
  int L1 = ROI[2][K], L2 =  ROI[3][K]; \
  \
  for (I = ROI[0][K]; I <= L1; I++) \
    for (J = ROI[1][K]; J <= L2; J++) \
      CODE_EXEC \
}
/* K must be positionned at current layer number, ROI[][] must be
   defined and set */


/* -----------------------------------------------------------------------------
   Generic ROI handling
----------------------------------------------------------------------------- */


/* checks wether given line is null, and if not, sets i_start to the
   index of first non-null point */
inline BOOL 
line_Occupied(ROI_T* line, int dim, int &start);


/* checks wether given layer column is null, and if not, sets start to
   the index of first non-null point */
inline BOOL 
column_Occupied(ROI_T** layer, int dim_x, int dim_y, 
		int col_number, int &start);


/* Determines the ROI on current layer in setting upper_left and
   lower_right */
void 
find_ROI_layer(ROI_T** layer, int dim_X, int dim_Y, 
	       Point<int> &upper_left, Point<int> &lower_right);


/* merges both ROIs represented by their upper_left and lower_right
   corners into resulting ROI in setting upper_left_1 and
   lower_right_1 */
void 
merge_ROIs(Point<int> &upper_left_1, Point<int> &lower_right_1, 
	   const Point<int> &upper_left_2, const Point<int> &lower_right_2);


#endif /* _ROI_ */
