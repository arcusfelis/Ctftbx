#define SIGNAL           prhs[0]
#define TIME_INST        prhs[1]
#define NFFT             prhs[2]
#define WINDOWG          prhs[3]
#define WINDOWH          prhs[4]

#define MHS_OUT          plhs[0]
#define TIME_OUT         plhs[1]
#define FREQ_OUT         plhs[2]

#include "tftb.h"
#include "divers.c"
#include "create_window.c"
#include "mhs.c"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const
	     mxArray * prhs[])

{
  type_signal    Signal;
  type_TFR       tfr;
  double        *WindowG, *WindowH;
  double        *ptr_time_inst, *ptr_freq_vect, *time_out;
  int            WindowG_Length, WindowH_Length, i;


  /* tests the number of inputs and outputs */
  if ((nrhs < 1) || (nrhs > 5) || (nlhs < 1) || (nlhs >3))
    mexErrMsgTxt
      ("[TFR,t,f]=CTFRMHS(Signal,Time_instants,Number of frequencies,WindowG,WindowH)");



  /* recovery of the signal */
  Signal.length = (int) MAX (mxGetM (SIGNAL), mxGetN (SIGNAL));

  if (Signal.length == 0)
    mexErrMsgTxt ("The signal must not be an empty matrix");

  Signal.real_part = mxGetPr (SIGNAL);
  if (mxIsComplex (SIGNAL))
    {
      Signal.imag_part = mxGetPi (SIGNAL);
      Signal.is_complex = TRUE;
    }
  else
    {
      Signal.is_complex = FALSE;
    }

  /* if there are not the five inputs ... */

  if (nrhs >= 2)		/* the time instants are given */
    {
      tfr.N_time = (int) MAX (mxGetM (TIME_INST), mxGetN (TIME_INST));
      ptr_time_inst = mxGetPr (TIME_INST);
    }
  else
    /* the time instants are no given */
    {
      /* default : Time_instants = 1:Signal_Length */
      tfr.N_time = Signal.length;
      ptr_time_inst = NULL;
    }

  if (nrhs >= 3)		/* the number of frequencies is given */
    {
      tfr.N_freq = (int) mxGetScalar (NFFT);
    }
  else
    /* the number of frequencies is not given */
    {
      /* default : Nfft = length of the signal */
      tfr.N_freq = Signal.length;
    }



  if (nrhs >= 4)		/* the windowG is given */
    {
      WindowG_Length = (int) MAX (mxGetM (WINDOWG), mxGetN (WINDOWG));
      WindowG = mxGetPr (WINDOWG);
    }
  else
    /* the window is not given */
    {
      /* default : gaussian window of length : the next odd
         number after Signal_Length/10 */
      WindowG_Length = (int) (tfr.N_freq / 10.0);
      if (ISODD(WindowG_Length) == 0)
	{
	  WindowG_Length = WindowG_Length + 1;
	}
      WindowG = (double *) ALLOC (WindowG_Length, sizeof (double));

      create_window(HAMMING,WindowG_Length,NULL,0,WindowG);
    }


   if (nrhs >= 5)		/* the windowH is given */
    {
      WindowH_Length = (int) MAX (mxGetM (WINDOWH), mxGetN (WINDOWH));
      WindowH = mxGetPr (WINDOWH);
    }
  else
    /* the window is not given */
    {
      /* default : gaussian window of length : the next odd
         number after Signal_Length/4 */
      WindowH_Length = (int) (tfr.N_freq / 4.0);
      if (ISODD(WindowH_Length) == 0)
	{
	  WindowH_Length = WindowH_Length + 1;
	}
      WindowH = (double *) ALLOC (WindowH_Length, sizeof (double));

      create_window(HAMMING,WindowH_Length,NULL,0,WindowH);
    }

  /* case of the third output : FREQ_OUT vector of frequency bins */
  if (nlhs >= 3)
    {
      FREQ_OUT = mxCreateDoubleMatrix (1, tfr.N_freq, mxREAL);
      ptr_freq_vect = mxGetPr (FREQ_OUT);
    }
  else
    {
      ptr_freq_vect = NULL;
    }


  /* creation of the output real matrix */
  tfr.is_complex = FALSE;
  MHS_OUT = mxCreateDoubleMatrix (tfr.N_freq, tfr.N_time, mxREAL);
  mem_alloc_TFR (&tfr, ptr_freq_vect, ptr_time_inst, mxGetPr (MHS_OUT),
                 mxGetPi (MHS_OUT));

  /* creation of the time_instants vector if not given */
  if (ptr_time_inst == NULL)
    {
      for (i = 0; i < tfr.N_time; i++)
	{
	  tfr.time_instants[i] = i + 1.0;
	}
    }



  /* case of the second output */
  if (nlhs >= 2)
    {
      TIME_OUT = mxCreateDoubleMatrix (1, tfr.N_time, mxREAL);
      time_out = mxGetPr (TIME_OUT);
      for (i = 0; i < tfr.N_time; i++)
	{
	  time_out[i] = tfr.time_instants[i];
	}
    }

  /* ---------------------------------------------------- */
  /*             calls the computation function           */
  /* ---------------------------------------------------- */

  mhs(Signal, WindowG, WindowG_Length,  WindowH, WindowH_Length, tfr);


  /* ---------------------------------------------------- */
  /*                     Free the memory                  */
  /* ---------------------------------------------------- */

  if (nrhs < 2)
    {
      FREE (tfr.time_instants);
    }

  
  if (nlhs < 3)
    {
      FREE (tfr.freq_bins);
    }

}
