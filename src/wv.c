/* EXISTS AN INTERFACE PROGRAM TO MATLAB : CTFRWV.C                           *
 *============================================================================*
 * Name of the function : wv.c (void)                                         *
 * Authors              : Emmanuel Roy - Manuel DAVY                          *
 * Date of creation     : 10 - 02 - 1999                                      *
 *----------------------------------------------------------------------------*
 * THE ALGORITHM                                                              *
 *                                                                            *
 * Given a signal to analyze in time and frequency, computes the Wigner-Ville *
 * Time-Frequency Distribution (WV) :                                         *
 *                                                                            *
 *               /                        -j2pi f tau                         *
 *     WV(t,f) = | x(t+tau/2x*(t-tau/2) e             dtau                    *
 *              /                                                             *
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
 * tfr.is_complex      | must be set to FALSE (a WV tfr is real-valued)      *
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
 * taumax              | local time-lag variable bounds. Used to take into    *
 *                     | accound the beginning and the end of the             *
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
wv (type_signal Signal, type_TFR tfr)

{
  int            Nfft, column, row, time;
  int            taumax, tau;
  double        *lacf_real, *lacf_imag;	/* local autocorrelation function */

 /*--------------------------------------------------------------------*/
 /*                      Test the input variables                      */
 /*--------------------------------------------------------------------*/


   if (tfr.is_complex == TRUE)
    {
      printf ("wv.c : The tfr matrix must be real valued\n");
      exit(0);
    }

  if (tfr.N_freq <= 0)
    {
      printf ("wv.c : The field tfr.N_freq is not correctly set\n");
      exit(0);
    }

  if (tfr.N_time <= 0)
    {
      printf ("wv.c : The field tfr.N_time is not correctly set\n");
      exit(0);
    }


  /*--------------------------------------------------------------------*/
  /*           creation of the vector of frequency bins  (output)       */
  /*--------------------------------------------------------------------*/
  Nfft = po2 (tfr.N_freq);
  
  for (row = 0; row < tfr.N_freq; row++)
    {
      tfr.freq_bins[row] = (double) (0.5 * row) / tfr.N_freq;
    }
 /*--------------------------------------------------------------------*/
 /*       memory allocation for the local autocorrelation function     */
 /*--------------------------------------------------------------------*/
  lacf_real = (double *) ALLOC (tfr.N_freq , sizeof (double));
  lacf_imag = (double *) ALLOC (tfr.N_freq , sizeof (double));

 /* initialization of these vectors */
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

      /* time instants of interest to compute the Wigner distrib. */
      time = ((int) tfr.time_instants[column]) - 1;
      
      /* taumax enables the computation near the edges */
      taumax = MIN (time, (Signal.length - time - 1));
      taumax = MIN (taumax, (tfr.N_freq / 2 - 1));
      
      /* for each delay value, the laf is computed and ffted */
      /* in order to pass from the (t,tau) domain to the (t,f) */
      /* domain */
      for (tau = -taumax; tau <= taumax; tau++)
	{
	  row = irem((tfr.N_freq+tau), tfr.N_freq ) ;
	  /* when the signal is complex valued */
	  if (Signal.is_complex == TRUE)
	    {
	      lacf_real[row] =   Signal.real_part[time + tau]
		               * Signal.real_part[time - tau]
		             +   Signal.imag_part[time + tau]
                               * Signal.imag_part[time - tau];

 	      lacf_imag[row] =   Signal.imag_part[time + tau]
                               * Signal.real_part[time - tau]
		             -   Signal.real_part[time + tau]
                               * Signal.imag_part[time - tau];
	    }
	  /* when the signal is real valued */
	  else
	    {
	      lacf_real[row] =   Signal.real_part[time + tau]
                               * Signal.real_part[time - tau];

	      lacf_imag[row] = 0.0;
	    }
        }


       tau=floor(tfr.N_freq/2);
       if ((time<=Signal.length-tau-1)&(time>=tau))
       {
        if (Signal.is_complex == TRUE)
        {
         lacf_real[tau] =  Signal.real_part[time+tau]*Signal.real_part[time-tau]
                        +Signal.imag_part[time+tau]*Signal.imag_part[time-tau];
         lacf_imag[tau] = 0;
        }
        else
        {
         lacf_real[tau] =  Signal.real_part[time+tau]*Signal.real_part[time-tau];
         lacf_imag[tau] = 0;
        }
       }

      /* fft of the local autocorrelation function lacf */
      fft (tfr.N_freq, Nfft, lacf_real, lacf_imag);

      
      /* the fft is put in the wv matrix and reinitialized */
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
