/* EXISTS AN INTERFACE PROGRAM TO MATLAB : CTFRPWV.C                          *
 *============================================================================*
 * Name of the function : pwv.c (void)                                        *
 * Authors              : Emmanuel Roy - Manuel DAVY                          *
 * Date of creation     : 10 - 02 - 1999                                      *
 *----------------------------------------------------------------------------*
 * THE ALGORITHM                                                              *
 *                                                                            *
 * Given a signal to analyze in time and frequency, computes the Pseudo-Wigner*
 * Ville Distribution (PWV) :                                                 *
 *                                                                            *
 *                 /                                                          *
 *                 |                                 -j2pi f tau              *
 *     PWV(t,f) =  | h*(tau) x(t+tau/2) x*(t-tau/2) e            dtau         *
 *                 |                                                          *
 *                /                                                           *
 *                                                                            *
 *                                                                            *
 * This function is real valued. Its computation requires a real or complex   *
 * signal, a vector containing time instants, the number of frequency bins    *
 * and a frequency smoothing window.                                          *
 *                                                                            *
 *============================================================================*
 * INPUT VARIABLES                                                            *
 * Name                |              role                                    *
 * Signal              | The signal to analyze. No field modified             *
 *                     |                                                      *
 * Window              | Vector containing the frequency smoothing window     *
 * Window_Length       | Number of points of this window (ODD number !)       *
 *                     |                                                      *
 * tfr                 | Matrix containing the resulting TFR (real)           *
 * tfr.time_instants   | positions of the smoothing window                    *
 * tfr.N_time          | length of '.time_instants' = number of cols.         *
 *                     | in the tfr matrix                                    *
 * tfr.N_freq          | number of frequency bins = number of rows in the tfr *
 *                     | matrix                                               *
 * tfr.is_complex      | must be set to FALSE (a PWV tfr is real-valued)      *
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
 * half_Window_Length  | half-length of the frequency smoothing window        *
 * norm                | normalization factor for the frequency window        *
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
pwv (type_signal Signal,
    double *Window, int Window_Length,
    type_TFR tfr)

{
  int            Nfft, column, row, time;
  int            half_Window_Length;
  int            taumax, tau;
  double        *lacf_real, *lacf_imag;/* local autocorrelation function */
  double         norm;

 /*--------------------------------------------------------------------*/
 /*                      Test the input variables                      */
 /*--------------------------------------------------------------------*/


   if (tfr.is_complex == TRUE)
    {
      printf ("pwv.c : The tfr matrix must be real valued\n");
      exit(0);
    }

  if (tfr.N_freq <= 0)
    {
      printf ("pwv.c : The field tfr.N_freq is not correctly set\n");
      exit(0);
    }

  if (tfr.N_time <= 0)
    {
      printf ("pwv.c : The field tfr.N_time is not correctly set\n");
      exit(0);
    }

  if (ISODD(Window_Length) == 0)
    {
      printf ("pwv.c : The window Length must be an ODD number\n");
      exit(0);
    }

  /*--------------------------------------------------------------------*/
  /*                    Determines some internal constants              */
  /*--------------------------------------------------------------------*/
  half_Window_Length = (Window_Length - 1) / 2;
  norm=Window[half_Window_Length];

  /* normalization of the window */
  for(row=0;row<Window_Length;row++)
    {
      Window[row]=Window[row]/norm;
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
  /* memory allocation and init. of the local autocorrelation fuction   */
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

      /* time instants of interest to compute the tfr  */
      time = ((int) tfr.time_instants[column]) - 1;

      /* taumax limits the range of tau near the edges */
      taumax = MIN (time, (Signal.length - time - 1));
      taumax = MIN (taumax, (tfr.N_freq / 2 - 1));
      taumax = MIN (taumax, half_Window_Length);

      /* The signal is windowed around the current time */
      for (tau = -taumax; tau <= taumax; tau++)
	{
	  row = irem( (tfr.N_freq+tau), tfr.N_freq ) ;
	  if (Signal.is_complex == TRUE)
	    {
	      /* case of complex-valued signal */
	      lacf_real[row] =   (Signal.real_part[time + tau]
                                * Signal.real_part[time - tau]
                               +  Signal.imag_part[time + tau]
                                * Signal.imag_part[time - tau])
                              *   Window[tau+half_Window_Length];

	      lacf_imag[row] =   (Signal.imag_part[time + tau]
                                * Signal.real_part[time - tau]
                               -  Signal.real_part[time + tau]
                                * Signal.imag_part[time - tau])
                              *   Window[tau+half_Window_Length];
	    }
	  else
	    {
	      /* case of real_valued signal */
	      lacf_real[row] =   (Signal.real_part[time + tau]
                                * Signal.real_part[time - tau])
		              *   Window[tau+half_Window_Length];

	      lacf_imag[row] = 0.0;
          }
        }


     tau=floor(tfr.N_freq/2);
       if ((time<=Signal.length-tau-1)&(time>=tau)&(tau<=half_Window_Length))
       {
        if (Signal.is_complex == TRUE)
        {
         lacf_real[tau] =  0.5*( Window[half_Window_Length+tau]
                                *( Signal.real_part[time+tau]
                                * Signal.real_part[time-tau]
                                + Signal.imag_part[time+tau]
                                * Signal.imag_part[time-tau])
                               + Window[half_Window_Length-tau]*(
                               Signal.real_part[time-tau]
                               * Signal.real_part[time+tau]
                               + Signal.imag_part[time-tau]
                               * Signal.imag_part[time+tau])) ;

         lacf_imag[tau] = 0.5*( Window[half_Window_Length+tau]
                                *(Signal.real_part[time-tau]
                                * Signal.imag_part[time+tau]
                                - Signal.real_part[time+tau]
                                * Signal.imag_part[time-tau])
                               + Window[half_Window_Length-tau]*(
                                  Signal.real_part[time+tau]
                                  * Signal.imag_part[time-tau]
                                  + Signal.real_part[time-tau]
                                  * Signal.imag_part[time+tau]));
        }
        else
        {
         lacf_real[tau] =  0.5*(Window[half_Window_Length+tau]
                                 *(Signal.real_part[time+tau]
                                 * Signal.real_part[time-tau])
                               + Window[half_Window_Length-tau]
                               *( Signal.real_part[time-tau]
                               * Signal.real_part[time+tau])) ;
         lacf_imag[tau] = 0;
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
