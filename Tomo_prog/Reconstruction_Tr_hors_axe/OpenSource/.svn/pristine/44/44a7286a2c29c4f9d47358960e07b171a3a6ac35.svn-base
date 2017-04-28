#include "ROI.h"


inline BOOL 
line_Occupied(ROI_T* line, int dim, int &start)
{
  int i;

  for (i = 0; i < dim; i++)
    if (line[i]) 
      {
	start = i;
	return TRUE;
      }

  return FALSE;
}


/* checks wether given layer column is null, and if not, sets start to
   the index of first non-null point */
inline BOOL 
column_Occupied(ROI_T** layer, int dim_x, int dim_y, 
		int col_number, int &start)
{
  int i;
  assert(col_number < dim_x);

  for (i = 0; i < dim_y; i++)
    if (layer[i][col_number]) 
      {
	start = i;
	return TRUE;
      }
  
  return FALSE;
}


void 
find_ROI_layer(ROI_T** layer, int dim_X, int dim_Y, 
	       Point<int> &upper_left, Point<int> &lower_right)
{
  int i, j;
  int start;
  Point<int> upper, lower, left, right;

  for (i = 0; i < dim_X; i++)
    {
      if (line_Occupied(layer[i], dim_Y, start))
	{
	  upper.set(i, start, 0);
	  break;
	}
    }

  for (i = dim_X - 1; i >= 0; i--)
    {
      if (line_Occupied(layer[i], dim_Y, start))
	{
	  lower.set(i, start, 0);
	  break;
	}
    }

  for (j = 0; j < dim_Y; j++)
    {
      if (column_Occupied(layer, dim_X, dim_Y, j, start))
	{
	  left.set(start, j, 0);
	  break;
	}
    }

  for (j = dim_Y - 1; j >= 0; j--)
    {
      if (column_Occupied(layer, dim_X, dim_Y, j, start))
	{
	  right.set(start, j, 0);
	  break;
	}
    }

  upper_left.set(upper.get_X(), left.get_Y(), 0);
  lower_right.set(lower.get_X(), right.get_Y(), 0);
}


void 
merge_ROIs(Point<int> &upper_left_1, Point<int> &lower_right_1, 
	   const Point<int> &upper_left_2, const Point<int> &lower_right_2)
{
  Point<int> ul_t, lr_t;
  
  ul_t.set(MIN(upper_left_1.get_X(),upper_left_2.get_X()), 
	   MIN(upper_left_1.get_Y(),upper_left_2.get_Y()), 
	   0);
  lr_t.set(MAX(lower_right_1.get_X(),lower_right_2.get_X()), 
	   MAX(lower_right_1.get_Y(),lower_right_2.get_Y()), 
	   0);
  
  upper_left_1 = ul_t;
  lower_right_1 = lr_t;
}
