/* interface program between MATLAB and language C
  for Reassign.C */
#include "tftb.h"

#define TFR             prhs[0]
#define FIELD_X         prhs[1]
#define FIELD_Y         prhs[2]

#define TFR_REAS_OUT    plhs[0]

#include "divers.c"
#include "reassign.c"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const
	     mxArray * prhs[])

{
  double        *field_x, *field_y;
  type_TFR       TFR_to_reassign, TFR_reassigned;

  if (!(nrhs == 3))
    mexErrMsgTxt ("TFR_REAS=Ctfrreas(TFR,field_x,field_y)");

  if ((mxGetM (TFR) != mxGetM (FIELD_X)) || (mxGetM (TFR) != mxGetM
					     (FIELD_Y)) || (mxGetM (FIELD_X)
							    != mxGetM (FIELD_Y)))
    mexErrMsgTxt ("Time frequency matrix and field  must be the same size");

  if ((mxGetN (TFR) != mxGetN (FIELD_X)) || (mxGetN (TFR) != mxGetN
					     (FIELD_Y)) || (mxGetN (FIELD_X)
							    != mxGetN (FIELD_Y)))
    mexErrMsgTxt ("Time frequency matrix and field  must be the same size");


  /* Recovery of TFR to reassign */
  TFR_to_reassign.real_part = mxGetPr (TFR);
  TFR_to_reassign.is_complex = FALSE;

  /* recovery of the field */
  field_x = mxGetPr (FIELD_X);
  field_y = mxGetPr (FIELD_Y);

  /* TFR dimensions recovery */
  TFR_to_reassign.N_time = mxGetN (TFR);
  TFR_to_reassign.N_freq = mxGetM (TFR);




  /* Creation of the output variable */
  TFR_reassigned.N_time = TFR_to_reassign.N_time;
  TFR_reassigned.N_freq = TFR_to_reassign.N_freq;
  TFR_REAS_OUT = mxCreateDoubleMatrix (TFR_reassigned.N_freq,
				       TFR_reassigned.N_time, mxREAL);


  TFR_reassigned.is_complex = FALSE;
  /* pointer on the results dist */
  TFR_reassigned.real_part = mxGetPr (TFR_REAS_OUT);


  /* Reassignement of the TFR matrix according to the field of vectors */
  reassign (TFR_to_reassign, field_x, field_y, TFR_reassigned);

}
