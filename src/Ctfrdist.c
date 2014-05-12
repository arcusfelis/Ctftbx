/* Interface program between MATLAB and language C
for the progrma DISTANCE.C */


#include "tftb.h"

#define TFR1             prhs[0]
#define TFR2             prhs[1]
#define DIST_NAME        prhs[2]
#define DIST_COEF        prhs[3]

#define DISTANCE_OUT     plhs[0]

#include "divers.c"
#include "distance.c"

/* at the moment , complex TFRs are not accepted */

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const
	     mxArray * prhs[])

{
  type_TFR       first_TFR, second_TFR;
  double         coef, *dist;
  int            dist_name_length;
  int            N_time, N_freq, dist_name;
  char          *dist_name_string;

  if ((nrhs < 3) || (nrhs > 4))
    mexErrMsgTxt ("Dist=Ctfrdist(TFR1,TFR2,distance_name,distance_coef)");


   
  /* Recovery of TFRs */
  first_TFR.real_part = mxGetPr (TFR1);
  second_TFR.real_part = mxGetPr (TFR2);
  first_TFR.is_complex = FALSE;
  second_TFR.is_complex = FALSE;


  /* TFRs dimensions recovery */
  first_TFR.N_time = mxGetM (TFR1);
  first_TFR.N_freq = mxGetN (TFR1);
  second_TFR.N_time = mxGetM (TFR2);
  second_TFR.N_freq = mxGetN (TFR2);

  if ((first_TFR.N_time != second_TFR.N_time) ||
      (first_TFR.N_freq != second_TFR.N_freq))
    mexErrMsgTxt ("Time frequency matrices must be the same size");


  /* recovery of the ditance name */
  
  if(!mxIsChar(DIST_NAME))
    mexErrMsgTxt("Variable distance_name must contain a string.\n");
  
  dist_name_length = mxGetN(DIST_NAME)+1;
  dist_name_string = (char *) ALLOC (dist_name_length,sizeof(char));
  mxGetString(DIST_NAME,dist_name_string,dist_name_length);


  /* distance name recovery */
  dist_name = 0;
 

  if((!(strcmp(dist_name_string,"Lq"))) ||
     (!(strcmp(dist_name_string,"LQ"))) ||
     (!(strcmp(dist_name_string,"lq"))))
    {
      dist_name = LQ;
    }

  if((!(strcmp(dist_name_string,"Quadratic"))) ||
     (!(strcmp(dist_name_string,"QUADRATIC"))) ||
     (!(strcmp(dist_name_string,"quadratic"))))
     {
       dist_name = QUADRATIC;
     }

  if((!(strcmp(dist_name_string,"Correlation"))) ||
     (!(strcmp(dist_name_string,"CORRELATION"))) ||
     (!(strcmp(dist_name_string,"correlation"))))
     {
       dist_name = CORRELATION;
     }

  if((!(strcmp(dist_name_string,"Kolmogorov"))) ||
     (!(strcmp(dist_name_string,"KOLMOGOROV"))) ||
     (!(strcmp(dist_name_string,"kolmogorov"))))
     {
       dist_name = KOLMOGOROV;
     }
    
  if((!(strcmp(dist_name_string,"Kullback"))) ||
     (!(strcmp(dist_name_string,"KULLBACK"))) ||
     (!(strcmp(dist_name_string,"kullback"))))
     {
       dist_name = KULLBACK;
     }
  if((!(strcmp(dist_name_string,"Chernoff"))) ||
      (!(strcmp(dist_name_string,"CHERNOFF"))) ||
      (!(strcmp(dist_name_string,"chernoff"))))
     {
       dist_name = CHERNOFF;
     }
     
   if((!(strcmp(dist_name_string,"Matusita"))) ||
      (!(strcmp(dist_name_string,"MATUSITA"))) ||
      (!(strcmp(dist_name_string,"matusita"))))
     {
       dist_name = MATUSITA;
     }
 
 
   if((!(strcmp(dist_name_string,"NLq"))) ||
      (!(strcmp(dist_name_string,"NLQ"))) ||
      (!(strcmp(dist_name_string,"nlq"))))
     {
       dist_name = NLQ;
     }


   if((!(strcmp(dist_name_string,"Lsd"))) ||
      (!(strcmp(dist_name_string,"LSD"))) ||
      (!(strcmp(dist_name_string,"lsd"))))
     {
       dist_name = LSD;
     }
     
  if((!(strcmp(dist_name_string,"Jensen"))) ||
     (!(strcmp(dist_name_string,"JENSEN"))) ||
     (!(strcmp(dist_name_string,"jensen"))))
     {
       dist_name = JENSEN;
     }
 
  if (dist_name == 0)
    {
      mexErrMsgTxt ("Unknown distance name");
    }


  /* some tests about the distence coef */
  /* to have or not to have a coefficient ? */

  if ((dist_name == LQ)       || (dist_name == CHERNOFF) ||
      (dist_name == MATUSITA) ||(dist_name == NLQ) || 
      (dist_name == LSD)      || (dist_name == JENSEN))
    {
      if (nrhs == 3) /* cases where there should be a coef given */
	{
	  mexErrMsgTxt ("A coefficient is required for this distance");
	}
    }
 


  /* Distance coefficient recovery */
  if (nrhs == 4)
    {
      coef = mxGetScalar (DIST_COEF);
    }
  else
    {
      coef = 1.0;
    }


  /* Error cases */
  if ((dist_name == LQ) && (coef <= 0))
    mexErrMsgTxt ("distance_coef > 0 for Lq distance");

  if ((dist_name == CHERNOFF) && ((coef > 1) || (coef < 0)))
    mexErrMsgTxt ("0 <= distance_coef <= 1 for Chernoff distance");

  if ((dist_name == MATUSITA) && (coef < 1))
    mexErrMsgTxt ("distance_coef >= 1 for Generalized Matusita distance");

  if ((dist_name == NLQ) && (coef <= 0))
    mexErrMsgTxt ("distance_coef > 0 for Normalized Lq distance");

  if ((dist_name == LSD) && (coef <= 0))
    mexErrMsgTxt ("distance_coef > 0 for q-Log Spectral Deviation");

  if ((dist_name == JENSEN) && (coef <= 0))
    mexErrMsgTxt ("distance_coefq > 0 for q-Jensen divergence");

  if ((dist_name != LQ) && (dist_name != QUADRATIC) &&
      (dist_name != CORRELATION) && (dist_name != KOLMOGOROV) &&
      (dist_name != KULLBACK) && (dist_name != CHERNOFF) &&
      (dist_name != MATUSITA) && (dist_name != NLQ) &&
      (dist_name != LSD) && (dist_name != JENSEN))
    mexErrMsgTxt ("Bad distance identification number");


  /* Creation of the output variable */
  DISTANCE_OUT = mxCreateDoubleMatrix (1, 1, mxREAL);

  /* pointer on the results dist */
  dist = mxGetPr (DISTANCE_OUT);


  /* computation of the distance */
  distance (first_TFR, second_TFR, dist_name, coef, dist);

  FREE (dist_name_string);

}

