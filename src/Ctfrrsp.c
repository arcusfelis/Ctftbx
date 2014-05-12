/* interface program between MATLAB and language C
  for Reassign.C */

#define SIGNAL            prhs[0]
#define TIME_INST         prhs[1]
#define NFFT              prhs[2]
#define WINDOW            prhs[3]

#define TFR_REAS_OUT      plhs[0]
#define TFR_NOT_REAS_OUT  plhs[1]
#define FIELD_CPX_OUT     plhs[2]

#include "tftb.h"
#include "divers.c"
#include "create_window.c"
#include "reassign.c"
#include "stft.c"
#include "gradient.c"
#include "reas_spectro.c"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const
	     mxArray * prhs[])

{
  /*-----------------------------------------------------------------*/
  /*                           VARIABLES                             */
  /*-----------------------------------------------------------------*/
  /* inputs */
  type_signal    Signal;
  int            Window_Length;
  double        *Window;

  /* outputs */
  type_TFR       TFR_reassigned, TFR_not_reassigned;
  double        *field_time;
  double        *field_freq;

  /* internal */
  int            time;
  double        *ptr_time_inst, *ptr_tfr_not_reas;

  /* some error cases .... */
 
   /* tests the number of inputs and outputs */
  if ((nrhs < 1) || (nrhs > 4) || (nlhs < 1) || (nlhs > 3))
    mexErrMsgTxt
      ("[TFR_REAS,TFR,FIELD]=Ctfrrsp(signal,[time_instants,Nfft,window])");

  /*-----------------------------------------------------------------*/
  /*                      Recovery of the inputs                     */
  /*-----------------------------------------------------------------*/
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


  /* case of the structures TFR_reassigned and TFR_not_reassigned */


  /* if there are not the four inputs ... */
  if (nrhs >= 2)			/* the time instants are is given */
    {
      /* TFR dimensions recovery */
      TFR_reassigned.N_time = (int) MAX (mxGetM (TIME_INST), mxGetN (TIME_INST));
      /* allocation of memory for the 'time_instants' vectors */
      ptr_time_inst = mxGetPr (TIME_INST);
    }
  else
    /* the time instants are no given */
    {
      /* default : Time_instants = 1:Signal_Length */
      TFR_reassigned.N_time = Signal.length;

      /* allocation of memory for the 'time_instants' vectors */
      ptr_time_inst = NULL;
    }


  if (nrhs >= 3)  /* the Nfft is given */
    {
      TFR_reassigned.N_freq = (int) mxGetScalar (NFFT);
    }
  else
    /* the number of frequencies is not given */
    {
      /* default : Nfft = length of the signal */
      TFR_reassigned.N_freq = Signal.length;
    }


  TFR_reassigned.freq_bins = (double *) ALLOC (TFR_reassigned.N_freq, sizeof(double));
  TFR_reassigned.is_complex = FALSE;



  if (nrhs >= 4)     /* the window is given */
    {
      /* recovery of the window and its size */
      Window = mxGetPr (WINDOW);
      Window_Length = (int) MAX (mxGetM (WINDOW), mxGetN (WINDOW));
    }
  else
    {
      /* default : Hamming window of length : the next odd
         number after Signal_Length/4 */
      Window_Length = (int) (TFR_reassigned.N_freq / 4.0);
      if (ISODD(Window_Length) == 0) /* not odd */
	{
	  Window_Length = Window_Length + 1;
	}
      Window = (double *) ALLOC (Window_Length, sizeof (double));
      create_window(HAMMING,Window_Length,NULL,0,Window);
    }



  /*-----------------------------------------------------------------*/
  /*                  CHECKS FOR THE OUTPUT VARIABLES                */
  /*-----------------------------------------------------------------*/

  /* Creation of the output variable */
  TFR_REAS_OUT = mxCreateDoubleMatrix (TFR_reassigned.N_freq,
				       TFR_reassigned.N_time, mxREAL);

  /* allocation of memory for TFR_reassigned : the frequency bins */
  /* have to be created - the time instants vector depends on wether */
  /* it is given by the user or not */
  mem_alloc_TFR(&TFR_reassigned,
		NULL, ptr_time_inst, mxGetPr(TFR_REAS_OUT), NULL);

  /* creation of the time_instants vector if not given */
  if (ptr_time_inst == NULL)
    {
      for (time = 0; time < TFR_reassigned.N_time; time++)
	{
	  TFR_reassigned.time_instants[time] = time + 1.0;
	}
    }


  /* case of the non-reassigned tfr */
  TFR_not_reassigned.N_freq = TFR_reassigned.N_freq;
  TFR_not_reassigned.N_time = TFR_reassigned.N_time;
  TFR_not_reassigned.is_complex = FALSE;

 
  /* if more than one output is required */
  if (nlhs >= 2) /* the unreassigned spectro is passed in the matlab environment */
    {
      /* creation of the matlab variable */
      TFR_NOT_REAS_OUT = mxCreateDoubleMatrix (TFR_not_reassigned.N_freq,
					       TFR_not_reassigned.N_time, mxREAL);

      /* pointer on the output matrix */
      ptr_tfr_not_reas = mxGetPr (TFR_NOT_REAS_OUT);
    }
  else /* the unreassigned spectro is not required */
    {
      ptr_tfr_not_reas = NULL;
    }

  /* allocation of memory for TFR_not_reassigned : the frequency bins */
  /* are allready created - the time instants vector too */
  /* it is given by the user or not */
  mem_alloc_TFR(&TFR_not_reassigned,
 		TFR_reassigned.freq_bins, TFR_reassigned.time_instants,
		ptr_tfr_not_reas, NULL);

 
  if (nlhs >= 3) /* the field of reassignement is passed in the matlab environment */
    {
      /* creation of the matlab variable */
      FIELD_CPX_OUT = mxCreateDoubleMatrix (TFR_reassigned.N_freq,
					    TFR_reassigned.N_time,
					    mxCOMPLEX);

      /* pointer on the output matrices  */
      field_time = mxGetPr (FIELD_CPX_OUT);
      field_freq = mxGetPi (FIELD_CPX_OUT);
    }
  else
    /* need to create the variables field_time and field_freq */
    {
      field_time = 
	(double *) ALLOC (TFR_reassigned.N_time * TFR_reassigned.N_freq,
			  sizeof (double));
      field_freq = 
	(double *) ALLOC (TFR_reassigned.N_time * TFR_reassigned.N_freq,
			  sizeof (double));

    }
 


  /*-----------------------------------------------------------------*/
  /*                       COMPUTATION                               */
  /*-----------------------------------------------------------------*/
   reas_spectro (Signal,		/*inputs */ 
		Window, Window_Length,  
		TFR_reassigned, TFR_not_reassigned,	/*outputs */ 
		field_time, field_freq); 

  /*-----------------------------------------------------------------*/
  /*                         FREE MEMORY                             */
  /*-----------------------------------------------------------------*/

  if (nrhs < 4)			/* when no window is given */
    {
      FREE (Window);
    }

  if (nrhs < 2) /* when no time instant vector is given */
    {
      FREE (TFR_reassigned.time_instants);
    }

  if (nlhs < 3) /* the reas. field is not required */
    {
      FREE (field_time); 
      FREE (field_freq);
    }

  if (nlhs < 2) /* the unreassigned spectro is not required */
    {
      FREE (TFR_not_reassigned.real_part);
    }

   FREE(TFR_reassigned.freq_bins);
 
}
