/* interface program between MATLAB and language C
for the program AF2TFR.C */


#include "tftb.h"

#define AF               prhs[0]
#define KERNEL           prhs[1]

#define TFR_OUT          plhs[0]

#include "divers.c"
#include "af2tfr.c"


void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const
	     mxArray * prhs[])
{
  type_AF        ambif, kernel;
  type_TFR       tfr;

  if (!(nrhs == 2))
    mexErrMsgTxt ("tfr=Caf2tfr(AF,kernel)");

  if (!(mxIsComplex (AF)))
    mexErrMsgTxt ("AF must be complex valued");
  /* recovery of the AF */
  ambif.real_part = mxGetPr (AF);
  ambif.imag_part = mxGetPi (AF);
  ambif.N_doppler = mxGetN (AF);
  ambif.N_delay = mxGetM (AF);
  ambif.is_complex = TRUE;

  kernel.real_part = mxGetPr (KERNEL);
  kernel.N_doppler = mxGetN (KERNEL);
  kernel.N_delay = mxGetM (KERNEL);
  kernel.is_complex = FALSE;

  if ((kernel.N_doppler != ambif.N_doppler) || (kernel.N_delay != ambif.N_delay))
    mexErrMsgTxt
      ("The Ambiguity function and the kernel must be the same size");


  tfr.is_complex = FALSE;
  tfr.N_time = ambif.N_doppler;
  tfr.N_freq = ambif.N_delay;
  TFR_OUT = mxCreateDoubleMatrix (tfr.N_freq, tfr.N_time, mxREAL);

  /* pointer on the results dist */
  tfr.real_part = mxGetPr (TFR_OUT);
  /* computation of the distance */
  af2tfr (ambif, kernel, tfr);
}
