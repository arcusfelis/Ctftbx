/* EXISTS AN INTERFACE PROGRAM TO MATLAB : CTFRREAS.C                */
/*====================================================================*
 * Name of the function : reassign (double)                           *
 * Author               : Manuel DAVY                                 * 
 * Date of creation     : 26 - 01 - 1999                              *
 *--------------------------------------------------------------------*
 * THE ALGORITHM                                                      *
 * Given a TFR and a field of vectors  , reassigns the TFR pixels     *
 * according to the field along the directions 'time' (x) and 'freq'  *
 * (y). sometimes, the pixel should be reassigned outside the TFR     *
 * matrix. In time, such pixels are left along the edges (time=0  or  *
 * time = max_time). In frequency, a circular rotation is done in     *
 * order to reassign the pixels somewhere.                            *
 *====================================================================*
 * INPUT VARIABLES                                                    *
 * Name              |                role                            *
 * TFR_to_reassign : | the TFR to be reassigned                       *
 *   .N_time         | Number of time bins i.e number of columns      *
 *   .N_time         | Number of frequency bins i.e. number of rows   *
 *   .is_complex     | must be initialized to the same value as       *
 *                   | 'TFR_reassigned.is_complex'                    *
 *   .real_part      | real part of the tfr to reassign               *
 *   .imag_part      | imag. part of the tfr (only if the field       *
 *                   | 'TFR_to_reassign.is_complex==TRUE'             *
 *                   |                                                *
 * TFR_reassigned    | result of the reassignment. The previous fileds*
 *                   | '.N_time', '.N_freq' and '.is_complex' in      *
 *                   | TFR_to_reassign) must be set to the same values*
 *                   |                                                *
 * field_x           | component of the reassignment along the axis   *
 *                   | 'time'                                         *
 * field_y           | component of the reassignment along the axis   *
 *                   | 'frequency'                                    *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES                                                   *
 * Name              |                role                            *
 * TFR_reassigned    | Result of the reassignment                     *
 *     .real_part    | real part of the reassigned TFR                *
 *     .imag_part    | imag part of the reassigned TFR                *
 *--------------------------------------------------------------------*
 * LOCAL VARIABLES                                                    *
 * Name              |               role                             *
 * time, freq        |  row  and column index in the matrix TFR       *
 * index             |  index in the matrix seen as a vector          *
 * index_time        |new position of the TFR pixel located in 'time' *
 * index_freq        |new position of the TFR pixel located in 'freq' *
 *====================================================================*
 * SUBROUTINES USED HERE                                              *
 *--------------------------------------------------------------------*
 * Name   | int idx(int line, int row, int nb_row)                    *
 * Action | computes the vector index for an element in a matrix given* 
 *        | the line and column indices and the number of columns     *
 * Place  | divers.c                                                  *
 *====================================================================*/


void
reassign (type_TFR TFR_to_reassign, double *field_x, double *field_y,
	  type_TFR TFR_reassigned)
{
  int            index_time, index_freq;
  int            time, freq, index;

 /*--------------------------------------------------------------------*/
 /*                      Test the input variables                      */
 /*--------------------------------------------------------------------*/
  if ((TFR_to_reassign.is_complex != TFR_reassigned.is_complex) ||
      (TFR_to_reassign.N_time != TFR_reassigned.N_time) || 
      (TFR_to_reassign.N_freq != TFR_reassigned.N_freq))
    {
      printf ("reassign.c :  TFR_reassigned and TFR_to_reassign are not compatible\n");
      exit(0);
    }


 /*--------------------------------------------------------------------*/
 /*               initialization of the reassigned matrix              */
 /*--------------------------------------------------------------------*/
  for (index = 0; index < (TFR_to_reassign.N_time * TFR_to_reassign.N_freq);
       index++)
    {
      TFR_reassigned.real_part[index] = 0.0;
    }

  if (TFR_reassigned.is_complex == TRUE)
    {
      for (index = 0;
	   index < (TFR_to_reassign.N_time *TFR_to_reassign.N_freq);
	   index++)
	{
	  TFR_reassigned.imag_part[index] = 0.0;
	}
    }
 /*--------------------------------------------------------------------*/
 /*                    reassignement of the matrix                     */
 /*--------------------------------------------------------------------*/
  for (time = 0; time < TFR_to_reassign.N_time; time++)
    for (freq = 0; freq < TFR_to_reassign.N_freq; freq++)
      {
	index = idx (freq, time, TFR_to_reassign.N_freq);
	/* computation of the final position of the pixel located */
        /* in (time, freq) */
	index_time = ROUND (time + field_x[index]);
	index_freq = ROUND (freq + field_y[index]);

	/* case of the negative indices            */
	if (index_time < 0) /* in time, the pixels are reassigned along */
	  index_time = 0   ;/* the edges */

	
	while (index_freq < 0) /* in frequency, the pixels are reassigned */
	  {                 /* according to a rotation */
	    index_freq = index_freq  + TFR_to_reassign.N_freq;
	  }


	/* case of too large indices */
	if  (index_time >= TFR_to_reassign.N_time)
	  /* in time, the pixels are reassigned along the edges*/
	  {
	    index_time = TFR_to_reassign.N_time - 1;
	  }
	while (index_freq >= TFR_to_reassign.N_freq)
	  /* in frequency, the pixels are reassigned */
	  {     /* according to a rotation */
	    index_freq = index_freq  - TFR_to_reassign.N_freq;
	  }
 


	/* reassignement */
	index = idx (index_freq, index_time, TFR_to_reassign.N_freq);
	TFR_reassigned.real_part[index] = TFR_reassigned.real_part[index]
	  + TFR_to_reassign.real_part[idx (freq, time, TFR_reassigned.N_freq)];

	if (TFR_to_reassign.is_complex == TRUE)
	  {
	    TFR_reassigned.imag_part[index] = TFR_reassigned.imag_part[index]
	      + TFR_to_reassign.imag_part[idx (freq, time, TFR_to_reassign.N_freq)];

	  }
      }

}
