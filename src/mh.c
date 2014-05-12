/* EXISTS AN INTERFACE PROGRAM TO MATLAB : CTFRMH.C                           *
 *============================================================================*
 * Name of the function : mh.c (void)                                         *
 * Authors              : Emmanuel Roy - Manuel DAVY                          *
 * Date of creation     : 10 - 02 - 1999                                      *
 *----------------------------------------------------------------------------*
 * THE ALGORITHM                                                              *
 *                                                                            *
 * Given a signal to analyze in time and frequency, computes the Margenau-Hill*
 * Time-Frequency Distribution (MH) :                                         *
 *                                                                            *
 *                   /                                                        *
 *                1  |                                   -j2pi f tau          *
 *     MH(t,f) =  -  | ( x(t+tau)x*(t) + x(t)x*(t-tau) )e            dtau     *
 *                2  |                                                        *
 *                  /                                                         *
 *                                                                            *
 * This function is real valued. Its computation requires a real or complex   *
 * signal, a vector containing time instants and the number of frequency bins.*
 *                                                                            *
 *============================================================================*
 * INPUT VARIABLES                                                            *
 * Name                |              role                                    *
 * Signal              | The signal to analyze. No field modified             *
 *                     |                                                      *
 * tfr                 | Matrix containing the resulting TFR (real)           *
 * tfr.time_instants   | positions of the smoothing window                    *
 * tfr.N_time          | length of '.time_instants' = number of cols.         *
 *                     | in the tfr matrix                                    *
 * tfr.N_freq          | number of frequency bins = number of rows in the tfr *
 *                     | matrix                                               *
 * tfr.is_complex      | must be set to FALSE (a MH tfr is real-valued)      *
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
 * tau                 | time-lag variable                                    *
 * taumin              | local time-lag variable bounds. Used to take into    *
 * taumax              | accound the beginning and the end of the             *
 *                     | signal, where the window is cut                      *
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
mh (type_signal Signal, type_TFR tfr)

{
  int            Nfft, column, row, time;
  int            taumin, taumax, tau;
  double        *lacf_real, *lacf_imag;		/* local autocorrelation function */

 /*--------------------------------------------------------------------------*/
 /*                      Test the input variables                            */
 /*--------------------------------------------------------------------------*/


   if (tfr.is_complex == TRUE)
    {
      printf ("mh.c : The tfr matrix must be complex valued\n");
      exit(0);
    }

  if (tfr.N_freq <= 0)
    {
      printf ("mh.c : The field tfr.N_freq is not correctly set\n");
      exit(0);
    }

  if (tfr.N_time <= 0)
    {
      printf ("mh.c : The field tfr.N_time is not correctly set\n");
      exit(0);
    }


 /*--------------------------------------------------------------------------*/
 /*             creation of the vector of frequency bins  (output)           */
 /*--------------------------------------------------------------------------*/
  Nfft = po2 (tfr.N_freq);
  
  for (row = 0; row < tfr.N_freq; row++)
    {
      tfr.freq_bins[row] = (double) row / tfr.N_freq;
    }
 /*--------------------------------------------------------------------------*/
 /*                  memory allocation for the windowed signal               */
 /*--------------------------------------------------------------------------*/
  lacf_real = (double *) ALLOC (tfr.N_freq , sizeof (double));
  lacf_imag = (double *) ALLOC (tfr.N_freq , sizeof (double));
  
  /* initialization of the intermediary vectors */
  for (row = 0; row < tfr.N_freq ; row++)
    {
      lacf_real[row] = 0.0;
      lacf_imag[row] = 0.0;
    }
  
 /*--------------------------------------------------------------------------*/
 /*          computation of the fft for the current windowed signal          */
 /*--------------------------------------------------------------------------*/
  for (column = 0; column < tfr.N_time; column++)
    {
      
      /* time instants of interest to compute the tfr */
      time = ((int) tfr.time_instants[column]) - 1;
      
      /* taumax and taumin enable the computation near the edges */
      taumin=MIN( (tfr.N_freq - time), (Signal.length - time - 1) );
      taumax = time;
      
      /* The signal is windowed around the current time */
      for (tau = -taumin; tau <= taumax; tau++)
	{
	  row = irem( (tfr.N_freq+tau),  tfr.N_freq ) ;
	  
	  if (Signal.is_complex == TRUE)
	  /* when the signal is complex valued */
	    {
	      lacf_real[row] =    Signal.real_part[time]
		                * Signal.real_part[time - tau]
		               +  Signal.imag_part[time]
                                * Signal.imag_part[time - tau];

	      lacf_imag[row] =    Signal.imag_part[time]
                                * Signal.real_part[time - tau]
                               -  Signal.real_part[time]
                                * Signal.imag_part[time - tau];
	    }
	  else
	  /* when the signal is real valued */
	    {
	      lacf_real[row] =    Signal.real_part[time]
                                * Signal.real_part[time - tau];

	      lacf_imag[row] = 0.0;
	    }
        }
      
      
      /* fft of the local autocorrelation function lacf */
      fft (tfr.N_freq, Nfft, lacf_real, lacf_imag);
      
      
      /* the fft is put in the tfr matrix  */
      for (row = 0; row < tfr.N_freq; row++)
	{
	  tfr.real_part[idx (row,column,tfr.N_freq)]= lacf_real[row];
	  lacf_real[row] = 0.0;
	  lacf_imag[row] = 0.0;
        }
    }
 /*--------------------------------------------------------------------------*/
 /*                  free the memory used in this program                    */
 /*--------------------------------------------------------------------------*/
  FREE (lacf_real);
  FREE (lacf_imag);
}
