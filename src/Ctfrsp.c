#define SIGNAL           prhs[0]
#define TIME_INST        prhs[1]
#define NFFT             prhs[2]
#define WINDOW           prhs[3]

#define SP_OUT           plhs[0]
#define TIME_OUT         plhs[1]
#define FREQ_OUT         plhs[2]
#define NORM_OUT         plhs[3]

#include "tftb.h"
#include "divers.c"
#include "create_window.c"
#include "sp.c"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const
	     mxArray * prhs[])

{
  type_signal    Signal;
  type_TFR       tfr;
  double        *Window, *norm_vector;
  double        *ptr_time_inst, *ptr_freq_vect, *time_out;
  int            Window_Length, i;


  /* tests the number of inputs and outputs */
  if ((nrhs < 1) || (nrhs > 4) || (nlhs < 1) || (nlhs >4))
    mexErrMsgTxt
      ("[TFR,t,f,norm]=CTFRSP(Signal,Time_instants,Number of frequencies,Window)");



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

  /* if there are not the four inputs ... */

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



  if (nrhs >= 4)		/* the window is given */
    {
      Window_Length = (int) MAX (mxGetM (WINDOW), mxGetN (WINDOW));
      Window = mxGetPr (WINDOW);
    }
  else
    /* the window is not given */
    {
      /* default : gaussian window of length : the next odd
         number after Signal_Length/4 */
      Window_Length = (int) (tfr.N_freq / 4.0);
      if (ISODD(Window_Length) == 0)
	{
	  Window_Length = Window_Length + 1;
	}
      Window = (double *) ALLOC (Window_Length, sizeof (double));

      create_window(HAMMING,Window_Length,NULL,0,Window);
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


  if (nlhs >= 4)
    {
      NORM_OUT = mxCreateDoubleMatrix (1, tfr.N_time, mxREAL);
      norm_vector = mxGetPr (NORM_OUT);
    }
  else
    {
      norm_vector = (double *) ALLOC (tfr.N_time, sizeof (double));
    }

  /* creation of the output real matrix */
  tfr.is_complex = FALSE;
  SP_OUT = mxCreateDoubleMatrix (tfr.N_freq, tfr.N_time, mxREAL);
  mem_alloc_TFR (&tfr, ptr_freq_vect, ptr_time_inst, mxGetPr (SP_OUT),
                 mxGetPi (SP_OUT));

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

  sp(Signal, Window, Window_Length, tfr, norm_vector);


  /* ---------------------------------------------------- */
  /*                     Free the memory                  */
  /* ---------------------------------------------------- */

  if (nrhs < 2)
    {
      FREE (tfr.time_instants);
    }

  if (nrhs < 4)
    {
      FREE (Window);
    }

  if (nlhs < 3)
    {
      FREE (tfr.freq_bins);
    }

  if (nlhs < 4)
    {
      FREE (norm_vector);
    }
}
