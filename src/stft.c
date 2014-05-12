/* EXISTS AN INTERFACE PROGRAM TO MATLAB : STFTMEX.C		             *
 *====================================================================*
 * Name of the function : stft.c (void)                  	          *
 * Author               : Manuel DAVY - IRCyN                         *
 * Date of creation     : 10 - 02 - 1999                              *
 *--------------------------------------------------------------------*
 * THE ALGORITHM             				                               *
 *								                                              *
 * Given a signal to analyze in time and frequency, computes the      *
 * Short Time Fourier Transform (STFT) :			                      *
 *								                                              *
 *		            /            -j 2 pi f s			                   *
 *    STFT(t,f) = | x(s)h(t-s)e      	    ds			                *
 *		            /						                                  *
 *								                                              *
 * This function is complex valued. Its computation requires a window,*
 * its displacement positions and the number of frequency bins to be  *
 * computed with discrete variables.				                      *
 *								                                              *
 *====================================================================*
 * INPUT VARIABLES   					                                  *
 * Name              |                role               	          *
 * Signal            |   The signal to analyze. No field modified     *
 *                   |          				                            *
 * Window            |   Vector containing the points of the window   *
 * Window_Length     |   Number of points of the window (ODD number !)*
 *                   |          				                            *
 * tfr               |  Matrix containing the resulting TFR (complex) *
 * tfr.time_instants |  positions of the smoothing window             *   
 * tfr.N_time        |  length of '.time_instants' = number of cols.  *
 *                   |  in the tfr matrix                             *
 * tfr.N_freq        |  number of frequency bins = number of rows     *
 *                   |  in the tfr matrix                             *
 * tfr.is_complex    |  must be set to TRUE (a stft is complex !)     *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES    						                               *
 * Name              |                role                	          *
 * tfr.real_part     |  the output tfr matrix  (real_part)            *
 * tfr.imag_part     |  the output tfr matrix  (imag part)            *
 * norm_vector       |  Value of the normalization factor applied at  *
 *                   |  the points of computation of the stft i.e.    *
 *                   |  tfr.time_instants 			                   *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES 						                               *
 * Name              |                 role                 	       *
 *                   |        				                            *
 *  Nfft             | Next power of two to tfr.N_freq                *
 *                   |               				                      *
 * column, line      |   variables of displacement in the matrices    *
 * tau               |   local time variable (the 's' in the equation *
 *                   |   above  				                            *
 * time              |   Current instant of computation of the spectro*
 * taumin            |   local time variable bounds. Used to take into*
 * taumax            |   accound the beginning and the end of the     *
 *                   |   signal, where the window is cut	             *
 * normh             |   current normalization value : depends on     *
 *                   |   wether the window is cut or not (near the    *
 *                   |   edges) 				                            *
 * wind_sig_real     |   Real and imaginary parts of the windowed     *
 * wind_sig_imag     |   signal at the current position of the window *
 * inter*            |   several intermediary variables		          *
 *====================================================================*
 * SUBROUTINES USED HERE				      	                         *
 *--------------------------------------------------------------------*
 * Name   | int idx(int line, int row, int nb_row)                    *
 * Action | computes the vector index for an element in a matrix given*
 *        | the line and column indices and the number of lines       *
 * Place  | divers.c                                                  *
 *--------------------------------------------------------------------*
 * Name   | int irem( double x, double y)                             *
 * Action | computes the remainder after Euclidean division of double *
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

void
stft (type_signal Signal,
      double *Window, int Window_Length,
      type_TFR tfr, double *norm_vector)

{
  int            Mfft, Nfft, column, row, time;
  int            taumin, taumax, tau;
  int            half_Window_Length;
  double        *wind_sig_real, *wind_sig_imag;		/* windowed signal */
  double         inter_real, inter_imag, normh;
  double         inter;

 /*--------------------------------------------------------------------*/
 /*                      Test the input variables                      */
 /*--------------------------------------------------------------------*/


   if (tfr.is_complex == FALSE)
    {
      printf ("stft.c : The tfr matrix must be complex valued\n");
      exit(0);
    }

  if (tfr.N_freq <= 0)
    {
      printf ("stft.c : The field tfr.N_freq is not correctly set\n");
      exit(0);
    }

  if (tfr.N_time <= 0)
    {
      printf ("stft.c : The field tfr.N_time is not correctly set\n");
      exit(0);
    }

 /*--------------------------------------------------------------------*/
 /*                   checks that the window length is odd             */
 /*--------------------------------------------------------------------*/

  if (ISODD(Window_Length) == 0)
    {
      printf ("stft.c : The window Length must be an ODD number\n");
      exit(0);
    }

  half_Window_Length = (Window_Length - 1) / 2;
  inter = 0.0;
  for (row = 0; row <Window_Length; row++)
	 {
	  inter = inter + sqr (Window[row]);
	 }
  normh = sqrt (inter);
 /*--------------------------------------------------------------------*/
 /*           creation of the vector of frequency bins  (output)       */
 /*--------------------------------------------------------------------*/
  Nfft = po2 (tfr.N_freq);

  for (row = 0; row < tfr.N_freq; row++)
    {
      tfr.freq_bins[row] = (double) row / tfr.N_freq;
    }
 /*--------------------------------------------------------------------*/
 /*                memory allocation for the windowed signal           */
 /*--------------------------------------------------------------------*/
  wind_sig_real = (double *) ALLOC (tfr.N_freq, sizeof (double));
  wind_sig_imag = (double *) ALLOC (tfr.N_freq, sizeof (double));

 /*--------------------------------------------------------------------*/
 /*      computation of the fft for the current windowed signal        */
 /*--------------------------------------------------------------------*/
  for (column = 0; column < tfr.N_time; column++)
    {

      /* initialization of the intermediary vectors */
      for (row = 0; row < tfr.N_freq; row++)
	{
	  wind_sig_real[row] = 0.0;
	  wind_sig_imag[row] = 0.0;
	}

      /* time instants of interest to compute the stft */
      time = ((int) tfr.time_instants[column]) - 1;

      /* the signal is multipied by the window between the instants
         time-taumin and time+taumax */
      /* when the window is wider than the number of desired frequencies (tfr.N_freq),
         the range is limited to tfr.N_freq */
      taumin = MIN (tfr.N_freq / 2, half_Window_Length);
      taumin = MIN (taumin, time);

      taumax = MIN ((tfr.N_freq / 2 - 1), half_Window_Length);
      taumax = MIN (taumax, (Signal.length - time - 1));

      /* Computation of a normalization factor, 
         equal to the quadratic norm of the window */


      norm_vector[column] = 1.0 / normh;
      /* The signal is windowed around the current time */
      for (tau = -taumin; tau <= taumax; tau++)
	{
	  row = irem( (tfr.N_freq+tau), tfr.N_freq ) ;
	  wind_sig_real[row] = Signal.real_part[time + tau]
	                     * Window[half_Window_Length + tau]
                             / normh;
	  if (Signal.is_complex == TRUE)
	    {
	      wind_sig_imag[row] = Signal.imag_part[time + tau]
                                 * Window[half_Window_Length + tau] 
                                 / normh;
	    }
	}

      /* fft of the windowed signal */
      fft (tfr.N_freq, Nfft, wind_sig_real, wind_sig_imag);


      /* the first half of the fft is put in the stft matrix  */
      for (row = 0; row < tfr.N_freq; row++)
	{
	  tfr.real_part[idx (row, column, tfr.N_freq)] = wind_sig_real[row];
	  tfr.imag_part[idx (row, column, tfr.N_freq)] = wind_sig_imag[row];

	}
    }
 /*--------------------------------------------------------------------*/
 /*                free the memory used in this program                */
 /*--------------------------------------------------------------------*/
  FREE (wind_sig_real);
  FREE (wind_sig_imag);

}
