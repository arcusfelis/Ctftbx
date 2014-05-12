/* interface program between MATLAB and kernel.c
SEE kernel.c for comments */

#include "tftb.h"


/* input parameters */
#define N_DOPPLER        prhs[0]
#define N_DELAY          prhs[1]
#define KERNEL_NAME      prhs[2]
#define PARAMETERS       prhs[3]
/* output - result */
#define KERNEL_OUT       plhs[0]




#include "divers.c"
#include "af.c"
#include "kernel.c"



void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const
	     mxArray * prhs[])

{
  double        *parameters;
  type_AF        ker;
  int            kernel_name, nb_param, kernel_name_length;
  char          *kernel_name_string;

  /* checks the number of inputs */
  if (!(nrhs == 4))
    mexErrMsgTxt ("ker=Ctfrker(NDoppler,NDelay,kernel_name,parameters)");


  /* Recovery of the kernel dimensions */
  ker.N_delay = mxGetScalar (N_DELAY);
  ker.N_doppler = mxGetScalar (N_DOPPLER);

  if ((ker.N_delay < 1) || ( ker.N_doppler < 1))
    {
      mexErrMsgTxt("Invalid number of rows/columns.");
    }
  ker.is_complex = FALSE;

  /* Which kernel shape ? */
  if(!mxIsChar(KERNEL_NAME))
    {
      mexErrMsgTxt("Variable kernel_name must contain a string.\n");
    }

  kernel_name_length = mxGetN(KERNEL_NAME)+1;
  kernel_name_string = (char *) ALLOC (kernel_name_length,sizeof(char));
  mxGetString(KERNEL_NAME,kernel_name_string,kernel_name_length);

  /* kernel name recovery */
  kernel_name = 0;
 
  if((!(strcmp(kernel_name_string,"MTEK"))) ||
     (!(strcmp(kernel_name_string,"mtek"))) ||
     (!(strcmp(kernel_name_string,"Mtek"))))
    {
      kernel_name = MTEK;
    }
  if((!(strcmp(kernel_name_string,"GMCWK"))) ||
     (!(strcmp(kernel_name_string,"gmcwk"))) ||
     (!(strcmp(kernel_name_string,"Gmcwk"))))
    {
      kernel_name = GMCWK;
    }
  if((!(strcmp(kernel_name_string,"RGK"))) ||
     (!(strcmp(kernel_name_string,"rgk"))) ||
     (!(strcmp(kernel_name_string,"Rgk"))))
    {
      kernel_name = RGK;
    }
  if((!(strcmp(kernel_name_string,"WV"))) ||
     (!(strcmp(kernel_name_string,"Wv"))) ||
     (!(strcmp(kernel_name_string,"wv"))))
    {
      kernel_name = WIGNER;
    }
  if((!(strcmp(kernel_name_string,"SPECTRO"))) ||
     (!(strcmp(kernel_name_string,"spectro"))) ||
     (!(strcmp(kernel_name_string,"Spectro"))))
    {
      kernel_name = SPECTRO;
    }

  if (kernel_name == 0)
    {
      mexErrMsgTxt ("Unknown kernel type");
    }



  /* recovers the parameters from inputs */
  parameters = mxGetPr (PARAMETERS);
  nb_param = (int) MAX (mxGetM (PARAMETERS), mxGetN (PARAMETERS));

  /* error cases */
  if (kernel_name == MTEK)
    {
      if (nb_param != NB_PARAM_MTEK)
	mexErrMsgTxt ("MTEK: 7 parameters are required");
      if (!(((BETA == 1) && (GAMMA == 1)) || ((BETA == 2) && (GAMMA == 0.5))))
	mexErrMsgTxt ("MTEK: beta=gamma=1 or beta=2 and gamma=0.5");
    }
  if (kernel_name == GMCWK)
    {
      if (nb_param < 2)
	mexErrMsgTxt ("GMCWK: at least one branch is required");
    }
  if (kernel_name == RGK)
    {
      if (nb_param < 3)
	mexErrMsgTxt
	  ("RGK: at least one pair of Fourier desciptors required");
    }
  if (kernel_name == SPECTRO)
    {
      if ((int) (nb_param / 2.0) != (nb_param / 2.0))
	mexErrMsgTxt ("Spectrogramme : the window length must be even");
      if (nb_param > ker.N_delay)
	mexErrMsgTxt
	  ("Spectrogramme : the window length must smaller than N_Delay ");


    }

  /* Creation of the output variable (matrix) */
  KERNEL_OUT = mxCreateDoubleMatrix (ker.N_doppler, ker.N_delay, mxREAL);

  /* pointer on the results matrix */
  ker.real_part = mxGetPr (KERNEL_OUT);


  /* computation of the kernel - call to kernel.c */
  kernel (kernel_name, parameters, nb_param, ker);


  FREE ( kernel_name_string);
}
