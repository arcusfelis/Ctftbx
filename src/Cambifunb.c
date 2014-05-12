/* interface program between Matlab and C language
for AF.C */



#include "tftb.h"

/* inputs */
#define SIGNAL     prhs[0]
#define TAU        prhs[1]
#define N_DOPPLER  prhs[2]

/* outputs */
#define AF_OUT     plhs[0]
#define TAU_OUT    plhs[1]
#define DOPPL_OUT  plhs[2]

#include "divers.c"
#include "af.c"

void
mexFunction (int nlhs, mxArray * plhs[],
	     int nrhs, const mxArray * prhs[])

{
  type_signal    Signal;
  type_AF        AF;
  double        *ptr_delay_vect, *ptr_doppler_vect;
  int            i, rem;
  double        *tau_out;

  if ((nrhs < 1) || (nrhs > 3) || (nlhs < 1) || (nlhs > 3))
    mexErrMsgTxt ("[AF,tau,Xi] = ambimex(X,TAU,N)");

  /* recovery of the signal */
  Signal.length = (int) MAX (mxGetM (SIGNAL), mxGetN (SIGNAL));

  if (Signal.length == 0)
    mexErrMsgTxt ("The signal must not be an empty matrix");

  /* tests whether Signal.length is even or odd */
  rem = ISODD(Signal.length);

  /* recovers data of the signal */
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



  /* recovery of the vector of delay values */
  if (nrhs >= 2)		/* the vector of delays is given */
    {
      AF.N_delay = (int) MAX (mxGetM (TAU), mxGetN (TAU));
      ptr_delay_vect = mxGetPr (TAU);
    }
  else
    /* the vector of delays is not given */
    {
      AF.N_delay = Signal.length - 1 + rem;
      ptr_delay_vect = NULL;
    }


  /* recovery of the number of doppler bins */
  if (nrhs >= 3)		/* the number of doppler bins is given */
    {
      AF.N_doppler = mxGetScalar (N_DOPPLER);
    }
  else
    /* the number of doppler bins is not given */
    {
      AF.N_doppler = Signal.length;
    }

  /* case of the third output : XI -> vector of doppler values */
  if (nlhs >= 3) /* the vector XI is given */ 
    {
      DOPPL_OUT = mxCreateDoubleMatrix (1, AF.N_doppler, mxREAL);
      ptr_doppler_vect = mxGetPr (DOPPL_OUT);
    }
  else /* the vector XI is not given */
    {
      ptr_doppler_vect = NULL;
    }

  /* the case of the compulsary output : the AF matrix */

  AF.is_complex = TRUE;
  AF_OUT = mxCreateDoubleMatrix (AF.N_doppler, AF.N_delay, mxCOMPLEX);


  /* allocation of memory for the AF matrix according to the pointers required */
  mem_alloc_AF (&AF, ptr_doppler_vect, ptr_delay_vect, mxGetPr (AF_OUT),
		mxGetPi (AF_OUT));




  /* creation of the vector of delay values in case of no specification */
  if (nrhs < 2)			/* the vector of delays is not given */
    {
      for (i = 0; i < AF.N_delay; i++)
	{
	  AF.delay_bins[i] = -(AF.N_delay - 1.0) / 2.0 + i;
	}
    }

  /* case of the second output */
  if (nlhs >= 2) /* the vector of delay bins is a required output */
    {
      TAU_OUT = mxCreateDoubleMatrix (1, AF.N_delay, mxREAL);
      tau_out = mxGetPr (TAU_OUT);
      for (i = 0; i < AF.N_delay; i++)
	{
	  tau_out[i] = AF.delay_bins[i];
	}
    }

  /* ---------------------------------------------------- */
  /*             calls the computation function           */
  /* ---------------------------------------------------- */
  af (Signal, AF);

  /* ---------------------------------------------------- */
  /*                     Free the memory                  */
  /* ---------------------------------------------------- */

  if (nrhs < 2)			/* the vector of delays is locally used */
    {
      FREE (AF.delay_bins);
    }

  if (nlhs < 3)			/* the vector of doppler is locally used */
    {
      FREE (AF.doppler_bins);
    }
}
