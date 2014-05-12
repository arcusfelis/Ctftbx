/* EXISTS AN INTERFACE PROGRAM TO MATLAB : CTFRZAM.C                          *
 *============================================================================*
 * Name of the function : zam.c (void)                                        *
 * Authors              : Emmanuel Roy - Manuel DAVY                          *
 * Date of creation     : 10 - 02 - 1999                                      *
 *----------------------------------------------------------------------------*
 * THE ALGORITHM                                                              *
 *                                                                            *
 * Given a signal to analyze in time and frequency, computes the              *
 * Zhao-Atlas-Marks Distribution (ZAM) :                                      *
 *                                                                            *
 *            /         / t+|tau|/2                                           *
 *            |        |                                 -j2pi f tau          *
 * ZAM(t,f) = | h(tau) |         x(mu+tau/2)x*(mu-tau/2)e            dmu dtau *
 *            |        |                                                      *
 *           /        / t-|tau|/2                                             *
 *                                                                            *
 * This function is real valued. Its computation requires a real or complex   *
 * signal, a vector containing time instants, the number of frequency bins, a *
 * time smoothing window and a frequency smoothing window.                    *
 *                                                                            *
 *============================================================================*
 * INPUT VARIABLES                                                            *
 * Name                |              role                                    *
 * Signal              | The signal to analyze. No field modified             *
 *                     |                                                      *
 * WindowT             | Vector containing the points of the time moothing    *
 *                     | window                                               *
 * WindowT_Length      | Number of points of the time window (ODD number !)   *
 *                     |                                                      *
 * WindowF             | Vector containing the points of the frequency window *
 * WindowF_Length      | Number of points of the window (ODD number !)        *
 *                     |                                                      *
 * tfr                 | Matrix containing the resulting TFR (real)           *
 * tfr.time_instants   | positions of the smoothing window                    *
 * tfr.N_time          | length of '.time_instants' = number of cols.         *
 *                     | in the tfr matrix                                    *
 * tfr.N_freq          | number of frequency bins = number of rows in the tfr *
 *                     | matrix                                               *
 * tfr.is_complex      | must be set to FALSE (a ZAM tfr is real-valued)      *
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
 * half_WindowT_Length | half-length of the time smoothing window             *
 * normT               | normalization factor for the time smoothing window   *
 *                     |                                                      *
 * half_WindowF_Length | half-length of the frequency smoothing window        *
 * normF               | normalization factor for the frequency window        *
 *                     |                                                      *
 * tau                 | time-lag variable                                    *
 * taumin              | local time-lag variable bounds. Used to take into    *
 * taumax              | accound the beginning and the end of the             *
 *                     | signal, where the window is cut                      *
 *                     |                                                      *
 * mu                  | time-smoothing variable                              *
 * mumin               | local time-smoothing variable bounds. Used to take   *
 * mumax               | into accound the beginning and the end of time       *
 *                     | smoothing procedure                                  *
 *                     |                                                      *
 * lacf_real           | real and imaginary parts of the local autocorrelation*
 * lacf_imag           | function of the signal                               *
 *                     |                                                      *
 * R1_real R1_imag     | used to compute real and imaginary parts of the time *
 * R2_real R2_imag     | smoothed-windowed local autocorrelation function     *
 *                     |                                                      *
 *============================================================================*
 * SUBROUTINES USED HERE                                                      *
 *----------------------------------------------------------------------------*
 * Name   | int idx(int i_row, int j_col, int nb_row)                         *
 * Action | computes the vector index for an element in a matrix given the row*
 *        | and column indices (i,j) and the total number of row              *
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
zam (type_signal Signal,
    double *WindowT, int WindowT_Length,
    double *WindowF, int WindowF_Length,
    type_TFR tfr)

{
  int            Nfft, column, row, time;
  int            half_WindowT_Length, half_WindowF_Length;
  int            taumin, taumax, tau;
  int            mumin, mumax, mu;
  double        *lacf_real, *lacf_imag;/* local autocorrelation function */
  double         normT, normF;
  double         R1_real, R1_imag, R2_real, R2_imag;

 /*--------------------------------------------------------------------*/
 /*                      Test the input variables                      */
 /*--------------------------------------------------------------------*/


   if (tfr.is_complex == TRUE)
    {
      printf ("zam.c : The tfr matrix must be real valued\n");
      exit(0);
    }

  if (tfr.N_freq <= 0)
    {
      printf ("zam.c : The field tfr.N_freq is not correctly set\n");
      exit(0);
    }

  if (tfr.N_time <= 0)
    {
      printf ("zam.c : The field tfr.N_time is not correctly set\n");
      exit(0);
    }

  if (ISODD(WindowT_Length) == 0)
    {
      printf ("zam.c : The time-window Length must be an ODD number\n");
      exit(0);
    }

  if (ISODD(WindowF_Length) == 0)
    {
      printf ("zam.c : The frequency-window Length must be an ODD number\n");
      exit(0);
    }
  /*--------------------------------------------------------------------*/
  /*                    Determines some internal constants              */
  /*--------------------------------------------------------------------*/

  half_WindowT_Length = (WindowT_Length - 1) / 2;
  normT=WindowT[half_WindowT_Length];
  /* normalization of the time smoothing window */
 for(row = 0; row < WindowT_Length; row++)
    {
      WindowT[row]=WindowT[row]/normT;
    }

  half_WindowF_Length = (WindowF_Length - 1) / 2;
  normF=WindowF[half_WindowF_Length];

  /* normalization of the frequency smoothing window */
  for(row = 0 ; row < WindowF_Length ; row++)
    {
      WindowF[row]=WindowF[row]/normF;
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

      taumax = MIN((time+half_WindowT_Length),(Signal.length-time-1+half_WindowT_Length));
      taumax = MIN(taumax,(tfr.N_freq / 2 - 1));
      taumax = MIN(taumax, half_WindowF_Length);

     /* initialization of the first local autocorrelation function */
      if (Signal.is_complex == TRUE)
	{
	  lacf_real[0] =   Signal.real_part[time]
	                 * Signal.real_part[time]
                       +   Signal.imag_part[time]
                         * Signal.imag_part[time];
	  /* the imag part is always zero because the imag part
	     of any complex 'x * conjugate(x)' is zero */

	  lacf_imag[0] =  0.0 ;
	}
      else /* the signal is real-valued */
	{
	  lacf_real[0] =   Signal.real_part[time]
                         * Signal.real_part[time];

	  lacf_imag[0] =   0.0;
	}

      /* The signal is windowed around the current time */
      for (tau = 1; tau <= taumax; tau++)
	{
	  R1_real = 0.0;
	  R2_real = 0.0;
	  R1_imag = 0.0;
	  R2_imag = 0.0;

	  /* bound of mu in order to take into account the edges */
	  mumin=MIN(tau,half_WindowT_Length);
	  mumin=MIN(mumin,(Signal.length-time-1-tau));

	  mumax=MIN(tau,half_WindowT_Length);
	  mumax=MIN(mumax,time-tau);

         for(mu=-mumin;mu<=mumax;mu++)
	   {
	     if (Signal.is_complex == TRUE)
	       {
		 R1_real = R1_real +   (Signal.real_part[time+tau-mu]
				      * Signal.real_part[time-tau-mu]
				    +   Signal.imag_part[time+tau-mu]
				      * Signal.imag_part[time-tau-mu])
		                   *    WindowT[half_WindowT_Length+mu];

		 R1_imag = R1_imag +    (Signal.imag_part[time+tau-mu]
				       * Signal.real_part[time-tau-mu]
				      -  Signal.real_part[time+tau-mu]
				       * Signal.imag_part[time-tau-mu])
                                    *    WindowT[half_WindowT_Length+mu];

		 R2_real = R2_real +    (Signal.real_part[time-tau-mu]
				       * Signal.real_part[time+tau-mu]
				      +  Signal.imag_part[time-tau-mu]
				       * Signal.imag_part[time+tau-mu])
		                    *    WindowT[half_WindowT_Length+mu];

		 R2_imag = R2_imag +   (Signal.imag_part[time-tau-mu]
				      * Signal.real_part[time+tau-mu]
                                     -  Signal.real_part[time-tau-mu]
				      * Signal.imag_part[time+tau-mu])
                                    *   WindowT[half_WindowT_Length+mu];
	       }
	     else
	       {
		 R1_real = R1_real +   (Signal.real_part[time+tau-mu]
				      * Signal.real_part[time-tau-mu])
                                     *  WindowT[half_WindowT_Length+mu];

		 R1_imag = 0.0;

		 R2_real = R2_real +   (Signal.real_part[time-tau-mu]
				      * Signal.real_part[time+tau-mu])
                                     *  WindowT[half_WindowT_Length+mu];

		 R2_imag = 0.0;
	       }
          }

         lacf_real[tau] = R1_real * WindowF[half_WindowF_Length+tau];
         lacf_imag[tau] = R1_imag * WindowF[half_WindowF_Length+tau];
         lacf_real[tfr.N_freq-tau] = R2_real
	                           * WindowF[half_WindowF_Length-tau];
         lacf_imag[tfr.N_freq-tau] = R2_imag
                                   * WindowF[half_WindowF_Length-tau];

        }


       tau=floor(tfr.N_freq/2);
       if ((time<=Signal.length-tau-1)&(time>=tau)&(tau<=half_WindowF_Length))
       {
        R1_real = 0.0;
	     R2_real = 0.0;
	     R1_imag = 0.0;
	     R2_imag = 0.0;

	  /* bound of mu in order to take into account the edges */
	  mumin=MIN(tau,half_WindowT_Length);
	  mumin=MIN(mumin,(Signal.length-time-1-tau));

	  mumax=MIN(tau,half_WindowT_Length);
	  mumax=MIN(mumax,time-tau);

         for(mu=-mumin;mu<=mumax;mu++)
	   {
	     if (Signal.is_complex == TRUE)
	       {
		 R1_real = R1_real +   (Signal.real_part[time+tau-mu]
				      * Signal.real_part[time-tau-mu]
				    +   Signal.imag_part[time+tau-mu]
				      * Signal.imag_part[time-tau-mu])
		                   *    WindowT[half_WindowT_Length+mu];

		 R1_imag = R1_imag +    (Signal.imag_part[time+tau-mu]
				       * Signal.real_part[time-tau-mu]
				      -  Signal.real_part[time+tau-mu]
				       * Signal.imag_part[time-tau-mu])
                                    *    WindowT[half_WindowT_Length+mu];

		 R2_real = R2_real +    (Signal.real_part[time-tau-mu]
				       * Signal.real_part[time+tau-mu]
				      +  Signal.imag_part[time-tau-mu]
				       * Signal.imag_part[time+tau-mu])
		                    *    WindowT[half_WindowT_Length+mu];

		 R2_imag = R2_imag +   (Signal.imag_part[time-tau-mu]
				      * Signal.real_part[time+tau-mu]
                                     -  Signal.real_part[time-tau-mu]
				      * Signal.imag_part[time+tau-mu])
                                    *   WindowT[half_WindowT_Length+mu];
	       }
	     else
	       {
		     R1_real = R1_real +   (Signal.real_part[time+tau-mu]
				      * Signal.real_part[time-tau-mu])
                                     *  WindowT[half_WindowT_Length+mu];

		    R1_imag = 0.0;

		    R2_real = R2_real +   (Signal.real_part[time-tau-mu]
				      * Signal.real_part[time+tau-mu])
                                     *  WindowT[half_WindowT_Length+mu];

		    R2_imag = 0.0;
	       }
        }

         lacf_real[tau] = 0.5*(R1_real * WindowF[half_WindowF_Length+tau]
                              + R2_real * WindowF[half_WindowF_Length-tau]);
         lacf_imag[tau] = 0.5*(R1_imag * WindowF[half_WindowF_Length+tau]
                             +R2_imag  * WindowF[half_WindowF_Length-tau]);
         
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
