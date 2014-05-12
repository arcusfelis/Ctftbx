#define SIGNAL           prhs[0]
#define WINDOW           prhs[1] /* WARNING : WINDOW is a MATRIX here */
#define TIME_INST        prhs[2]
#define NFFT             prhs[3]


#define MMCE_OUT         plhs[0]
#define TIME_OUT         plhs[1]
#define FREQ_OUT         plhs[2]

#include "tftb.h"
#include "divers.c"
#include "create_window.c"
#include "mmce.c"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const
	     mxArray * prhs[])

{
  type_signal    Signal;
  type_TFR       tfr;
  double        *Window;
  double        *ptr_time_inst, *ptr_freq_vect, *time_out;
  int            Window_Length, Window_col, i;


  /* tests the number of inputs and outputs */
  if ((nrhs < 1) || (nrhs > 4) || (nlhs < 1) || (nlhs >4))
    mexErrMsgTxt
      ("[TFR,t,f]=CTFRMMCE(Signal,Window[,Time_instants,Number of frequencies])");



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


    if (nrhs < 2)
    mexErrMsgTxt
      ("MMCEMEX(Signal,Window) : At least 2 parameters are required");


    Window_Length = (int) mxGetM (WINDOW);
    Window_col    = (int) mxGetN (WINDOW);
    if (Window_col <= 1)
         mexErrMsgTxt ("WINDOW_H must have at least two columns");
    Window = mxGetPr (WINDOW);


  if (nrhs >= 3)		/* the time instants are given */
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

  if (nrhs >= 4)		/* the number of frequencies is given */
    {
      tfr.N_freq = (int) mxGetScalar (NFFT);
    }
  else
    /* the number of frequencies is not given */
    {
      /* default : Nfft = length of the signal */
      tfr.N_freq = Signal.length;
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


  /* creation of the output complex matrix */
  tfr.is_complex = FALSE;
  MMCE_OUT = mxCreateDoubleMatrix (tfr.N_freq, tfr.N_time, mxREAL);
  mem_alloc_TFR (&tfr, ptr_freq_vect, ptr_time_inst, mxGetPr (MMCE_OUT),
                 mxGetPi (MMCE_OUT));

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

  mmce(Signal, Window, Window_Length, Window_col, tfr);


  /* ---------------------------------------------------- */
  /*                     Free the memory                  */
  /* ---------------------------------------------------- */

  if (nrhs < 3)
    {
      FREE (tfr.time_instants);
    }

  
  if (nlhs < 3)
    {
      FREE (tfr.freq_bins);
    }

  
}
