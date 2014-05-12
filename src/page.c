/* EXISTS AN INTERFACE PROGRAM TO MATLAB : CTFRPAGE.C                         *
 *============================================================================*
 * Name of the function : page.c (void)                                       *
 * Authors              : Emmanuel Roy - Manuel DAVY                          *
 * Date of creation     : 10 - 02 - 1999                                      *
 *----------------------------------------------------------------------------*
 * THE ALGORITHM                                                              *
 *                                                                            *
 * Given a signal to analyze in time and frequency, computes the Page         *
 * Distribution (PAGE) :                                                      *
 *                                                                            *
 *                           /  / t                        \*                 *
 *                           |  |         -j2pi f tau      |     -j2pi f t    *
 *   PAGE(t,f) = 2 Real{ x(t)|  | x(tau) e            dtau |   e           }  *
 *                           |  |                          |                  *
 *                           \  /-infty                   /                   *
 *                                                                            *
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
 * tfr.is_complex      | must be set to FALSE (a PAGE tfr is real-valued)     *
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
 *                     |                                                      *
 * lacf_real           | real and imaginary parts of the local autocorrelation*
 * lacf_imag           | function of the signal                               *
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
page(type_signal Signal, type_TFR tfr)

{
  int            Nfft, column, row, time;
  int            taumin, taumax, tau;
  double        *lacf_real, *lacf_imag; /*local autocorrelation function */

 /*--------------------------------------------------------------------*/
 /*                      Test the input variables                      */
 /*--------------------------------------------------------------------*/


   if (tfr.is_complex == TRUE)
    {
      printf ("page.c : The tfr matrix must be complex valued\n");
      exit(0);
    }

  if (tfr.N_freq <= 0)
    {
      printf ("page.c : The field tfr.N_freq is not correctly set\n");
      exit(0);
    }

  if (tfr.N_time <= 0)
    {
      printf ("page.c : The field tfr.N_time is not correctly set\n");
      exit(0);
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
  lacf_real = (double *) ALLOC (tfr.N_freq , sizeof (double));
  lacf_imag = (double *) ALLOC (tfr.N_freq , sizeof (double));
  
  /* initialization of the intermediary vectors */
  for (row = 0; row < tfr.N_freq ; row++)
    {
      lacf_real[row] = 0.0;
      lacf_imag[row] = 0.0;
   }
  
  /*--------------------------------------------------------------------*/
  /*      computation of the fft for the current windowed signal        */
  /*--------------------------------------------------------------------*/
  for (column = 0; column < tfr.N_time; column++)
    {
      
      /* time instants of interest to compute the tfr */
      time = ((int) tfr.time_instants[column]) - 1;

      /* taumin and taumax enable the computation near the edges */
      taumin = MIN( (tfr.N_freq - time-1), (Signal.length - time-1) );
      taumax = time;
      
      /* for each delay value, the laf is computed and ffted */
      /* in order to pass from the (t,tau) domain to the (t,f) */
      /* domain */
      for (tau = -taumin; tau <= taumax; tau++)
	{
	  row = irem( (tfr.N_freq+tau), tfr.N_freq ) ;
	  
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
  /*--------------------------------------------------------------------*/
  /*                free the memory used in this program                */
  /*--------------------------------------------------------------------*/
  FREE (lacf_real);
  FREE (lacf_imag);
}
