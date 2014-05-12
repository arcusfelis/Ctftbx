/* Interface program between MATLAB and language C
for the program CREATE_WINDOW.C */


#include "tftb.h"

#define LENGTH           prhs[0]
#define WINDOW_NAME      prhs[1]
#define PARAM1           prhs[2]
#define PARAM2           prhs[3]

#define WINDOW_OUT       plhs[0]

#include "divers.c"
#include "create_window.c"

/* at the moment , complex TFRs are not accepted */

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const
	     mxArray * prhs[])

{
 int       Window_length, Window_name, nb_param;
 int       Window_name_length;
 double   *param, *Window;
 char     *Window_name_string;


  /* checks the number of inputs */
  if ((nrhs < 2) || (nrhs > 4))
    mexErrMsgTxt ("Window=Ctfrwindow(length,'name'[,param1, param2])");


   
  /* recovery of the length */
  Window_length = mxGetScalar(LENGTH);


  Window_name_length = mxGetN(WINDOW_NAME)+1;
  Window_name_string = (char *) ALLOC (Window_name_length,sizeof(char));
  mxGetString(WINDOW_NAME,Window_name_string,Window_name_length);

  /* window name recovery */
  Window_name = 0;
 
  if((!(strcmp(Window_name_string,"HAMMING"))) ||
     (!(strcmp(Window_name_string,"hamming"))) ||
     (!(strcmp(Window_name_string,"Hamming"))))
    {
      Window_name = HAMMING;
    }
  if((!(strcmp(Window_name_string,"HANNING"))) ||
     (!(strcmp(Window_name_string,"hanning"))) ||
     (!(strcmp(Window_name_string,"Hanning"))))
    {
      Window_name = HANNING;
    }
  if((!(strcmp(Window_name_string,"NUTTALL"))) ||
     (!(strcmp(Window_name_string,"nuttall"))) ||
     (!(strcmp(Window_name_string,"Nuttall"))))
    {
      Window_name = NUTTALL;
    }
  if((!(strcmp(Window_name_string,"PAPOULIS"))) ||
     (!(strcmp(Window_name_string,"papoulis"))) ||
     (!(strcmp(Window_name_string,"Papoulis"))))
    {
      Window_name = PAPOULIS;
    }
  if((!(strcmp(Window_name_string,"HARRIS"))) ||
     (!(strcmp(Window_name_string,"harris"))) ||
     (!(strcmp(Window_name_string,"Harris"))))
    {
      Window_name = HARRIS;
    }
  if((!(strcmp(Window_name_string,"RECT"))) ||
     (!(strcmp(Window_name_string,"rect"))) ||
     (!(strcmp(Window_name_string,"Rect"))))
    {
      Window_name = RECTANG;
    }
  if((!(strcmp(Window_name_string,"TRIANG"))) ||
     (!(strcmp(Window_name_string,"triang"))) ||
     (!(strcmp(Window_name_string,"Triang"))))
    {
      Window_name = TRIANG;
    }
  if((!(strcmp(Window_name_string,"BARTLETT"))) ||
     (!(strcmp(Window_name_string,"bartlett"))) ||
     (!(strcmp(Window_name_string,"Bartlett"))))
    {
      Window_name = BARTLETT;
    }
  if((!(strcmp(Window_name_string,"BARTHANN"))) ||
     (!(strcmp(Window_name_string,"barthann"))) ||
     (!(strcmp(Window_name_string,"BartHann"))))
    {
      Window_name = BARTHANN ;
    }
  if((!(strcmp(Window_name_string,"BLACKMAN"))) ||
     (!(strcmp(Window_name_string,"blackman"))) ||
     (!(strcmp(Window_name_string,"Blackman"))))
    {
      Window_name = BLACKMAN ;
    }
  if((!(strcmp(Window_name_string,"GAUSS"))) ||
     (!(strcmp(Window_name_string,"gauss"))) ||
     (!(strcmp(Window_name_string,"Gauss"))))
    {
      Window_name = GAUSS;
    }
  if((!(strcmp(Window_name_string,"PARZEN"))) ||
     (!(strcmp(Window_name_string,"parzen"))) ||
     (!(strcmp(Window_name_string,"Parzen"))))
    {
      Window_name = PARZEN;
    }
  if((!(strcmp(Window_name_string,"DOLPH"))) ||
     (!(strcmp(Window_name_string,"dolph"))) ||
     (!(strcmp(Window_name_string,"Dolph"))))
    {
      Window_name = DOLPH;
    }
  if((!(strcmp(Window_name_string,"DOLF"))) ||
     (!(strcmp(Window_name_string,"dolf"))) ||
     (!(strcmp(Window_name_string,"Dolf"))))
    {
      Window_name = DOLF;
    }
  if((!(strcmp(Window_name_string,"HANNA"))) ||
     (!(strcmp(Window_name_string,"hanna"))) ||
     (!(strcmp(Window_name_string,"Hanna"))))
    {
      Window_name = HANNA;
    }
   if((!(strcmp(Window_name_string,"NUTBESS"))) ||
     (!(strcmp(Window_name_string,"nutbess"))) ||
     (!(strcmp(Window_name_string,"Nutbess"))))
    {
      Window_name = NUTBESS;
    }
  if((!(strcmp(Window_name_string,"SPLINE"))) ||
     (!(strcmp(Window_name_string,"spline"))) ||
     (!(strcmp(Window_name_string,"Spline"))))
    {
      Window_name = SPLINE;
    }
 if (Window_name == 0)
    {
      mexErrMsgTxt ("Unknown window type");
    }



  nb_param = 0;
  param = NULL;

  /* recovery of the parameters */ 
  if (nrhs == 4)
    {
      nb_param = 2;
      param = (double *) ALLOC (2, sizeof(double));
      param[0] = mxGetScalar(PARAM1);
      param[1] = mxGetScalar(PARAM2);

    }
 if (nrhs == 3)
    {
      nb_param = 1;
      param = (double *) ALLOC (1, sizeof(double));
      param[0] = mxGetScalar(PARAM1);
    }

 

  /* Creation of the output variable */
  WINDOW_OUT = mxCreateDoubleMatrix (Window_length, 1, mxREAL);

  /* pointer on the results dist */
  Window = mxGetPr (WINDOW_OUT);


  /* computation of the distance */
  create_window(Window_name, Window_length, param, nb_param, Window);

    
 if (nrhs > 2)
    {
      FREE (param);
    }

}

