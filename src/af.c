/* EXISTS AN INTERFACE PROGRAM TO MATLAB : AMBIMEX.C		      *
 *====================================================================*
 * Name of the function :  af.c (void)   			      *
 * Author               :  Manuel DAVY                                *
 * Date of creation     :  10 - 01 - 1999                             *
 *--------------------------------------------------------------------*
 * THE ALGORITHM             				              *
 *								      *
 * computes the instanteneous product by column (for a given delay)   *
 * and the Fourier transforms it to the ambiguity plane.	      *
 *								      *
 *====================================================================*
 * INPUT VARIABLES   					              *
 * Name              |                role               	      *
 * Signal            | the signal to be considered here               *
 * Signal.length     | must be initialized                            *
 * Signal.is_complex | must be initialized                            *
 * Signal.real_part  | must be initialized                            *
 * Signal.imag_part  | initialized if .is_complex == TRUE             *
 *                   |                                                *
 * AF                |   The matrix containing the AF.                *
 * AF.N_doppler      | number of rows in the final af matrix          *
 * AF.N_delay        | number of columns in the af matrix             *
 * AF.is_complex     | must be set to TRUE                            *
 * AF.delay_bins     | must specify the delay bins where the af has to*
 *                   | be computed 				      *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES    						      *
 * Name              |                role                	      *
 * AF.real_part      | Real part of the AF                            *
 * AF.imag_part      | Imaginary part of the AF                       *
 * AF.doppler_bins   | vector of doppler bins where the Af is computed*
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES 						      *
 * Name              |                role                 	      *
 * puiN              |   Next power of two to AF.N_doppler	      *
 * col, line         |   indices to scan the columns and lines of the *
 *                   |   matrices and vectors			      *
 * index1            |   indices to localize an element in a matrix   *
 * index2            |   (stored as a vector)			      *
 *                   |             				      *
 * afr,afi           |   vectors containing the current real and imag *
 *                   |   parts of the instantaneous product, and,     *
 *                   |   after the fft, the ambiguity function	      *
 *                   |                                                *
 * rem               |   remaining after the euclidean div. of the    *
 *                   |   length of the signal with 2                  *
 *====================================================================*
 * SUBROUTINES USED HERE				      	      *
 *--------------------------------------------------------------------*
 * Name   | int idx(int line, int row, int nb_row)                    *
 * Action | computes the vector index for an element in a matrix given*
 *        | the line and column indices and the number of lines       *
 * Place  | divers.c                                                  *
 *--------------------------------------------------------------------*
 * Name   | void fft(int n, int m, double *x, double *y)              *
 * Action | Computes the fft                                          *
 * Place  | divers.c                                                  *
 *--------------------------------------------------------------------*
 * Name   | int po2(int x)                                            *
 * Action | Computes the next power of two of x                       *
 * Place  | divers.c                                                  *
 *====================================================================*/


/*====================================================================*
 * THIS FUNCTION         				      	      *
 *--------------------------------------------------------------------*
 * Name   |   void af (type_signal Signal, type_AF AF)                *
 * Action | Computes the ambiguity function of a signal for the delay *
 *        | bins passed in the field AF.delay_bins                    *
 * Place  | af.c                                                      *
 *====================================================================*/

void
af (type_signal Signal, type_AF AF)

{
  int            puiN, tau, col, line;
  int            index1, index2, rem;
  double        *afr, *afi;


  /* tests the initialization of the fields */
  if (AF.is_complex == FALSE)
    {
      printf ("af.c : the AF matrix must be complex\n");
      exit(0);
    }

  if (AF.N_delay <= 0)
    {
      printf ("af.c : the field AF.N_delay is not correctly set \n");
      exit(0);
    }

  if (AF.N_doppler <= 0)
    {
      printf ("af.c : the field AF.N_doppler is not correctly set \n");
      exit(0);
    }

  /* when the delay_bins are not spaced form one point to the following in the signal
     there is an implicit decimation of the signal, but without filtering :
     DANGER for the user! */
  if (AF.delay_bins[0] + 1 != AF.delay_bins[1])
    {
      printf ("af.c : Warning -> this delay vector implies\n");
      printf ("       a non controled decimation of the signal !\n");
    }

  /* tests whether Signal.length is even or odd */
 
  rem = ISODD(Signal.length);

  /* memory allocation for an intermediary vector */

  afr = (double *) ALLOC (AF.N_doppler, sizeof (double));
  afi = (double *) ALLOC (AF.N_doppler, sizeof (double));


  puiN = po2 (AF.N_doppler);
  for (line = 0; line < AF.N_doppler; line++)
    {
      afr[line] = 0;
      afi[line] = 0;
    }
  /*construction of the AF matrix : by columns */
  for (col = 0; col < AF.N_delay; col++)
    {
      /* current value of the delay : tau */
      tau = (int) AF.delay_bins[col];

      /* for this delay : tau, computation of the instantaneous product
         for all the time values */
      for (line = (ABS (tau)); line < (AF.N_doppler - ABS (tau)); line++)
	{
	  /* when the signal is complex-valued */
	  if (Signal.is_complex == TRUE)
	    {
	      afr[line] = Signal.real_part[line + tau]
		* Signal.real_part[line - tau]
		+ Signal.imag_part[line + tau]
		* Signal.imag_part[line - tau];
	      afi[line] = Signal.real_part[line - tau]
		* Signal.imag_part[line + tau]
		- Signal.real_part[line + tau]
		* Signal.imag_part[line - tau];
	    }
	  else
	    /* when the signal is real-valued */
	    {
	      afr[line] = Signal.real_part[line + tau]
		* Signal.real_part[line - tau];
	      afi[line] = 0.0;
	    }
	}
      /* fft on the columns of the instantaneous product matrix for the current delay */
      fft (AF.N_doppler, puiN, afr, afi);

      /* recopy of the doppler vector in the af matrix */
      /* first half : copied in the second half of the AF */

       fftshift(afr,afr,AF.N_doppler); 
       fftshift(afi,afi,AF.N_doppler); 
       for (line=0;line<AF.N_doppler; line++)
	 {
	   index1 = idx (line, col, AF.N_doppler);
	   AF.real_part[index1] = afr[line];
	   AF.imag_part[index1] = afi[line];
	 }
  
  

      /* reinitialization of the vectors */
      for (line = 0; line < AF.N_doppler; line++)
	{
	  afr[line] = 0;
	  afi[line] = 0;
	}
    }

  /* updating of the .doppler_bins field */
  for (line = 0; line < AF.N_doppler; line++)
    {
      AF.doppler_bins[line] = -0.5 + (double) ((line + (rem / 2.0)) / AF.N_doppler);
    }


  /* don't forget the memory !! */
  FREE (afr);
  FREE (afi);

}
