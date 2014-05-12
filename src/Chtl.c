/* interface program between MATLAB and kernel.c
SEE kernel.c for comments */

#include "tftb.h"


/* input parameters */
#define TFR              prhs[0]
#define N_THETA          prhs[1]
#define N_RHO            prhs[2]
/* output - result */
#define HOUGH_OUT        plhs[0]
#define RHO_OUT          plhs[1]
#define THETA_OUT        plhs[2]





#include "divers.c"
#include "hough.c"



void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const
	     mxArray * prhs[])

{
  int            nb_theta, nb_rho;
  type_TFR       tfr;
  double        *transfo_hough, *rho_vect, *theta_vect;


  /* checks the number of inputs */
  if ((nrhs > 3)||(nrhs<1)||
      (nlhs > 3)||(nlhs<1))
    mexErrMsgTxt ("[HT,RHO,THETA]=Chtl(tfr,nb_theta,nb_rho)");


  /* Recovery of the tfr parameter */
  tfr.N_freq = mxGetN(TFR);
  tfr.N_time = mxGetM(TFR);
  tfr.real_part = mxGetPr(TFR);
  tfr.is_complex = FALSE;

  /* recovery of the other inputs */
  if(nrhs>1)
    {
      nb_theta = mxGetScalar(N_THETA);
    }
  else
    {
      nb_theta = tfr.N_time;
    }
 if(nrhs>2)
    {
      nb_rho = mxGetScalar(N_RHO);
    }
  else
    {
      nb_rho = tfr.N_freq;
    }

 /* recovery of the outputs */
 HOUGH_OUT = mxCreateDoubleMatrix(nb_rho,nb_theta,mxREAL);
 transfo_hough = mxGetPr(HOUGH_OUT);
 
 /* recovery or creation of the vectors of indices of the matrix */
 if(nlhs>1) /* rho_vect is required as an output */
   {
     RHO_OUT = mxCreateDoubleMatrix(1,nb_rho,mxREAL);
     rho_vect = mxGetPr(RHO_OUT);
   }
 else
   {
     rho_vect = ALLOC(nb_rho,sizeof(double));
   }

 if(nlhs>2) /* rho_vect is required as an output */
   {
     THETA_OUT = mxCreateDoubleMatrix(1,nb_theta,mxREAL);
     theta_vect = mxGetPr(THETA_OUT);
   }
 else
   {
     theta_vect = ALLOC(nb_theta,sizeof(double));
   }

 hough (tfr,nb_theta,nb_rho,transfo_hough,rho_vect,theta_vect);

 /* free the memory */
 if(nlhs<3) 
   {
     FREE(theta_vect);
   }
 if(nlhs<2) 
   {
     FREE(rho_vect);
   }

}
