/* EXISTS AN INTERFACCE PROGRAM WITH MATLAB : GRADMEX.C               *
 *====================================================================*
 * Name of the function : gradient (void)                             *
 * Author               : Manuel DAVY - IRCYN                         *
 * Date of creation     : 02 - 02 - 1999                              *
 *--------------------------------------------------------------------*
 * THE ALGORITHM                                                      *
 * Computes the 2D gradient of a 2D matrix of potential               *
 *          d matrix(x,y)                                             *
 * grad_x = ------------                                              *
 *             dx                                                     *
 *                                                                    *
 *         d matrix(x,y)                                              *
 * grad_y = ------------                                              *
 *             dy                                                     *
 * The algorithm consists in computing, for example in x , the central* 
 * slope :                                                            *
 *                                                                    *
 *                    matrix(i+1,j) - matrix (i-1,j)                  *
 * grad_x(i,j) = ----------------------------------------             *
 *                             2 * step_x                             *
 *                                                                    *
 * On the edge, the forward differnece is computed                    *
 *====================================================================*
 * INPUT VARIABLES                                                    *
 * Name           |                   role                            *
 * matrix         | Matrix of potential for which the gradient is     *
 *                | computed. Stored in a vector (in fact, a matrix   *
 *                | x-y)                                              *
 * size_x         | size of the matrix matrix along x (nb of lines)   *
 * size_y         | size of the matrix matrix along y (nb of columns) *
 * step_x         | step between two lines in Matrix                  *
 * step_y         | step between two columns in Matrix                *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES                                                   *
 * Name           |                   role                            *
 * grad_x         | Matrix (stored as a vector) for the gradient      *
 *                | along x                                           *
 * grad_y         | idem for y                                        *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES                                                 *
 * Name           |                   role                            *
 * x              | variable along the lines of matrix                *
 * y              | variable along the columns of matrix              *
 *====================================================================*/
void
gradient (double *matrix,
	  int size_x, int size_y,
	  double step_x, double step_y,
	  double *grad_x, double *grad_y)
{
  int            x, y;


  /*-----------------------------------------------*/
  /* for the edges, compute the forward difference */
  /*-----------------------------------------------*/

  /* make sure that the matrix is not a vector */
  /* in that case, the gradient over the unitary
     direction is null */
  if ((size_x <= 1) && (size_y <= 1))
    {
      grad_x[0] = 0;
      grad_y[0] = 0;
    }
  else
    {
      if (size_x > 1)		/* more than one line in matrix */
	{
	  for (y = 0; y < size_y; y++)
	    {
	      grad_x[idx (0, y, size_x)] = (matrix[idx (1, y, size_x)]
					    - matrix[idx (0, y, size_x)])
		/ step_x;
	      if (size_x >= 2)
		{
		  grad_x[idx (size_x - 1, y, size_x)] =
		    (matrix[idx (size_x - 1, y, size_x)]
		     - matrix[idx (size_x - 2, y, size_x)])
		    / step_x;
		}
	    }
	}
      else
	/* only one line in matrix */
	{
	  for (y = 0; y < size_y; y++)
	    {
	      grad_x[y] = 0;
	    }
	}
      if (size_y > 1)		/* more than one column in matrix */
	{
	  for (x = 0; x < size_x; x++)
	    {
	      grad_y[idx (x, 0, size_x)] = (matrix[idx (x, 1, size_x)]
					    - matrix[idx (x, 0, size_x)])
		/ step_y;
	      if (size_y >= 2)
		{
		  grad_y[idx (x, size_y - 1, size_x)] =
		    (matrix[idx (x, size_y - 1, size_x)]
		     - matrix[idx (x, size_y - 2, size_x)])
		    / step_y;
		}
	    }
	}
      else
	/* only one column in matrix */
	{
	  for (x = 0; x < size_x; x++)
	    {
	      grad_y[x] = 0;
	    }
	}
      /* for the other points, computes the central difference */
      if (size_x >= 2)
	{
	  for (x = 1; x < size_x - 1; x++)
	    for (y = 0; y < size_y; y++)
	      {
		grad_x[idx (x, y, size_x)] = (matrix[idx (x + 1, y, size_x)]
					      - matrix[idx (x - 1, y, size_x)])
		  / (2.0 * step_x);

	      }
	}
      if (size_y >= 2)
	{
	  for (y = 1; y < size_y - 1; y++)
	    for (x = 0; x < size_x; x++)
	      {
		grad_y[idx (x, y, size_x)] = (matrix[idx (x, y + 1, size_x)]
					      - matrix[idx (x, y - 1, size_x)])
		  / (2.0 * step_y);

	      }
	}
    }
}
