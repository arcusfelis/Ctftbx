/* EXISTS AN INTERFACE PROGRAM TO MATLAB : CTFRMHS.C                          *
 *============================================================================*
 * Name of the function : mhs.c (void)                                        *
 * Authors              : Emmanuel Roy - Manuel DAVY                          *
 * Date of creation     : 10 - 02 - 1999                                      *
 *----------------------------------------------------------------------------*
 * THE ALGORITHM                                                              *
 *                                                                            *
 * Given a signal to analyze in time and frequency, computes the Margenau-Hill*
 * Spectrogram Distribution (MHS) :                                           *
 *                                                                            *
 *                                                                            *
 *                       -1                                                   *
 *     MHS(t,f) = Real{ Kgh  F(t,f,g) F*(t,f,h) }                             *
 *                                                                            *
 *                  /                                                         *
 *    where  Kgh =  | h(u) g*(u) du                                           *
 *                 /                                                          *
 *                                                                            *
 *    F(t,f,h) and F(t,f,h) are the short-time Fourier transform of the signal*
 *    with h and g as analysis windows                                        *
 *                                                                            *
 * This function is real valued. Its computation requires a real or complex   *
 * signal, a vector containing time instants, the number of frequency bins, 2 *
 * frequency smoothing (analysis) windows.                                    *
 *                                                                            *
 *============================================================================*
 * INPUT VARIABLES                                                            *
 * Name                |              role                                    *
 * Signal              | The signal to analyze. No field modified             *
 *                     |                                                      *
 * WindowG             | Vector containing the points of the frequency window *
 * WindowG_Length      | Number of points of the g window (ODD number !)      *
 *                     |                                                      *
 * WindowH             | Vector containing the points of the frequency window *
 * WindowH_Length      | Number of points of the h window (ODD number !)      *
 *                     |                                                      *
 * tfr                 | Matrix containing the resulting TFR (real)           *
 * tfr.time_instants   | positions of the smoothing window                    *
 * tfr.N_time          | length of '.time_instants' = number of cols.         *
 *                     | in the tfr matrix                                    *
 * tfr.N_freq          | number of frequency bins = number of rows in the tfr *
 *                     | matrix                                               *
 * tfr.is_complex      | must be set to FALSE (a MHS tfr is real-valued)      *
 *                     |                                                      *
 *----------------------------------------------------------------------------*
 * OUTPUT VARIABLES                                                           *
 * Name                |                role                                  *
 * tfr.real_part       | the output tfr matrix  (real_part)                   *
 * tfr.freq_bins       | vector of frequency bins (freqs where the tfr matrix *
 *                     | is computed)                                         *
 *----------------------------------------------------------------------------*
 * INTERNAL VARIABLES                                                         *
 * Name                |                 role                                 *
 *                     |                                                      *
 * Nfft                | Next power of two to tfr.N_freq                      *
 * column, row         | variables of displacement in the matrices            *
 * time                | local time-instant variable to compute the tfr       *
 *                     |                                                      *
 * half_WindowG_Length | half-length of the frequency smoothing window g      *
 *                     |                                                      *
 * half_WindowH_Length | half-length of the frequency smoothing window h      *
 * normH               | normalization factor for the frequency window        *
 *                     |                                                      *
 * Lgh                 | variable to compute the lower length between         *
 *                     | half_WindowG_Length and half_WindowH_Length          *
 * Kgh                 | normalization factor                                 *
 *                     |                                                      *
 * tau                 | time-lag variable                                    *
 * taumin              | local time-lag variable bounds. Used to take into    *
 * taumax              | accound the beginning and the end of the             *
 *                     | signal, where the window is cut                      *
 *                     |                                                      *
 * windG_sig_real      | real and imaginary parts of the windowed signal      *
 * windG_sig_imag      | (analysis window g)                                  *
 *                     |                                                      *
 * windH_sig_real      | real and imaginary parts of the windowed signal      *
 * windH_sig_imag      | (analysis window h)                                  *
 *                     |                                                      *
 *============================================================================*
 * SUBROUTINES USED HERE                                                      *
 *----------------------------------------------------------------------------*
 * Name   | int idx(int i_row, int j_col, int nb_row)                         *
 * Action | computes the vector index for an element in a matrix given the row*
 *        | and column indices (i,j) and the total number of row              *
 * Place  | divers.c                                                          *
 *----------------------------------------------------------------------------*
 * Name   | int irem( double x, double y)                                     *
 * Action | computes the remainder after Euclidean division of double         *
 * Place  | divers.c                                                          *
 *----------------------------------------------------------------------------*
 * Name   | void fft(int n, int m, double *x, double *y)                      *
 * Action | Computes the fft                                                  *
 * Place  | divers.c                                                          *
 *----------------------------------------------------------------------------*
 * Name   | int po2(int x)                                                    *
 * Action | Computes the next power of two of x                               *
 * Place  | divers.c                                                          *
 *============================================================================*/

void
mhs (type_signal Signal,
      double *WindowG, int WindowG_Length,
      double *WindowH, int WindowH_Length,
      type_TFR tfr )

{
  int            Nfft, column, row, time;
  int            taumin, taumax, tau;
  int            half_WindowG_Length, half_WindowH_Length;
  double        *windG_sig_real, *windG_sig_imag; /* windowed signal */
  double        *windH_sig_real, *windH_sig_imag; /* windowed signal */
  double         normH;
  int            Lgh, points;
  double         Kgh;

 /*--------------------------------------------------------------------*/
 /*                      Test the input variables                      */
 /*--------------------------------------------------------------------*/


   if (tfr.is_complex == TRUE)
    {
      printf ("mhs.c : The tfr matrix must be real valued\n");
      exit(0);
    }

  if (tfr.N_freq <= 0)
    {
      printf ("mhs.c : The field tfr.N_freq is not correctly set\n");
      exit(0);
    }

  if (tfr.N_time <= 0)
    {
      printf ("mhs.c : The field tfr.N_time is not correctly set\n");
      exit(0);
    }

 /*--------------------------------------------------------------------*/
 /*                   checks that the window length is odd             */
 /*--------------------------------------------------------------------*/
 if (ISODD(WindowG_Length) == 0)
    {
      printf ("mhs.c : The window G Length must be an ODD number\n");
      exit(0);
    }


  if (ISODD(WindowH_Length) == 0)
    {
      printf ("mhs.c : The window H Length must be an ODD number\n");
      exit(0);
    }

  half_WindowG_Length = (WindowG_Length - 1) / 2;
  half_WindowH_Length = (WindowH_Length - 1) / 2;

  normH=WindowH[half_WindowH_Length];

  for(row = 0; row < WindowH_Length; row++)
    {
      WindowH[row] = WindowH[row]/normH;
    }

  
  Lgh = MIN( half_WindowG_Length , half_WindowH_Length );
  Kgh = 0.0;
  for( points=-Lgh ; points<=Lgh ; points++)
     {
      Kgh = Kgh + WindowH[half_WindowH_Length+points]
                 *WindowG[half_WindowG_Length+points];
      }

  for(row = 0; row < WindowH_Length; row++)
    {
      WindowH[row] = WindowH[row]/Kgh;
    }
  
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
  windG_sig_real = (double *) ALLOC (tfr.N_freq, sizeof (double));
  windG_sig_imag = (double *) ALLOC (tfr.N_freq, sizeof (double));

  windH_sig_real = (double *) ALLOC (tfr.N_freq, sizeof (double));
  windH_sig_imag = (double *) ALLOC (tfr.N_freq, sizeof (double));

  for (row = 0; row < tfr.N_freq; row++)
	   {
	    windG_sig_real[row] = 0.0;
	    windG_sig_imag[row] = 0.0;
       windH_sig_real[row] = 0.0;
	    windH_sig_imag[row] = 0.0;
       }
 /*--------------------------------------------------------------------*/
 /*      computation of the fft for the current windowed signal        */
 /*--------------------------------------------------------------------*/
  for (column = 0; column < tfr.N_time; column++)
    {

      /* time instants of interest to compute the tfr */
      time = ((int) tfr.time_instants[column]) - 1;

      taumin = MIN (tfr.N_freq / 2, half_WindowG_Length);
      taumin = MIN (taumin, time);

      taumax = MIN ((tfr.N_freq / 2 - 1), half_WindowG_Length);
      taumax = MIN (taumax, (Signal.length - time - 1));

      /* The signal is windowed around the current time */
      for (tau = -taumin; tau <= taumax; tau++)
	     {
	      row = irem( (tfr.N_freq+tau), tfr.N_freq ) ;
	      windG_sig_real[row] = Signal.real_part[time + tau]
	                         * WindowG[half_WindowG_Length + tau];
	      if (Signal.is_complex == TRUE)
	       {
	        windG_sig_imag[row] = Signal.imag_part[time + tau]
                                   * WindowG[half_WindowG_Length + tau];
	       }
        }

      /* fft of the windowed signal */
      fft (tfr.N_freq, Nfft, windG_sig_real, windG_sig_imag);


      taumin = MIN (tfr.N_freq / 2, half_WindowH_Length);
      taumin = MIN (taumin, time);

      taumax = MIN ((tfr.N_freq / 2 - 1), half_WindowH_Length);
      taumax = MIN (taumax, (Signal.length - time - 1));

      /* The signal is windowed around the current time */
      for (tau = -taumin; tau <= taumax; tau++)
	     {
	      row = irem( (tfr.N_freq+tau), tfr.N_freq ) ;
	      windH_sig_real[row] = Signal.real_part[time + tau]
	                         * WindowH[half_WindowH_Length + tau];
	      if (Signal.is_complex == TRUE)
	       {
	        windH_sig_imag[row] = Signal.imag_part[time + tau]
                                   * WindowH[half_WindowH_Length + tau];
	       }
        }

      /* fft of the windowed signal */
      fft (tfr.N_freq, Nfft, windH_sig_real, windH_sig_imag);



      /* the first half of the fft is put in the tfr matrix  */
      for (row = 0; row < tfr.N_freq; row++)
	{
	  tfr.real_part[idx (row,column,tfr.N_freq)] = 
	             windG_sig_real[row]*windH_sig_real[row]
                   + windG_sig_imag[row]*windH_sig_imag[row];

     windG_sig_real[row] = 0.0;
	    windG_sig_imag[row] = 0.0;
       windH_sig_real[row] = 0.0;
	    windH_sig_imag[row] = 0.0;
        }
    }
 /*--------------------------------------------------------------------*/
 /*                free the memory used in this program                */
 /*--------------------------------------------------------------------*/
  FREE (windG_sig_real);
  FREE (windG_sig_imag);
  FREE (windH_sig_real);
  FREE (windH_sig_imag);

}
