/*====================================================================*
 * this file contains several programs used in a signal processing    *
 * context. The programs are :                                        *
 * - fft       : Computes the fast Fourier transform of a complex     *
 *               signal                                               *
 * - po2       : Computes the next higher power of two                *
 * - idx       : In a matrix stored like a vector, computes the       *
 *               vector index corresponding to given line and column  *
 *               numbers in the matrix                                *
 * - ifft      : Computes the inverse fft                             *
 * - sqr       : square of a number                                   *
 * - powof     : Computes x to the power of alpha                     *
 * - Recover_Signal : In programs used inside the Matlab              *
 *                environement, recovers a signal from a matlab       *
 *                function                                            *
 * - gradient  : Computes the bidimensional gradient of a surface     *
 *               (called 'field of potential')                        *
 *====================================================================*/


/*====================================================================*
 * Name of the function : powof (double)                               *
 * Author               : Manuel DAVY - IRCYN                         *
 * Date of creation     : 02 - 02 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 *                                                                    *
 *  Computes x^alpha and takes the case x=0 into account              *
 *  !!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!              *
 *  x must be non-negative when alpha is not an interger              *
 *  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              *
 *====================================================================*/

double
powof (double x, double alpha)
{
  double         resu;



  if (x >= 0) /* in this case, no problem */
    {
      if (x == 0)
	resu = 0.0;
      else
	resu = exp (alpha * log (x));
    }
  else /* there may be problems */
    {
      if (alpha == (int)(alpha)) /* if alpha is an integer */
	{
	  /* then x^alpha is real-valued */ 
	  resu = powof ( -x, alpha);
	  /* and the sign of the results depends on
	     wether alpha is ODD or EVEN */
	  if (ISODD(alpha) == 1)
	    {
	      /* alpha is ODD */
	      resu = -resu;
	    }
	}
      else
	{
	  printf ("Attempt to compute x^alpha with x<0 : complex valued result\n");
	  exit(0);
	}
    }
  return resu;
}


/*                                                                    *
 *====================================================================*
 * Name of the function :                			      *
 * Author               :                                             *
 * Date of creation     :                                             *
 *--------------------------------------------------------------------*
 * THE ALGORITHM             				              *
								      *
								      *
 *====================================================================*
 * INPUT VARIABLES   					              *
 * Name        | type  |              role               	      *
 *             |       | kernel 				      *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES    						      *
 * Name        | type  |              role                	      *
 *								      *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES 						      *
 * Name        | type  |              role                 	      *
 *								      *
 *====================================================================*
 * SUBROUTINES USED HERE				      	      *
 *--------------------------------------------------------------------*
 * Name   |                                          	     	      *
 * Action |							      *
 * Place  | 							      *
 *====================================================================*/

void
fft (int Signal_Length, int Nfft, double *sig_real, double *sig_imag)
{
  /*------------------------------------------------------------------*/
  /*          when the signal length is a power of two                */
  /*------------------------------------------------------------------*/
  if (Signal_Length == (int) powof (2.0, Nfft) + 1)
    {
      /* local variables */
      int            i, j, k, n, n2;
      double         c, s, e, a, t1, t2;

      j = 0;			/* bit-reverse  */
      n2 = Signal_Length / 2;
      for (i = 1; i < Signal_Length - 1; i++)
	{
	  n = n2;
	  while (j >= n)
	    {
	      j = j - n;
	      n = n / 2;
	    }
	  j = j + n;

	  if (i < j)
	    {
	      t1 = sig_real[i];
	      sig_real[i] = sig_real[j];
	      sig_real[j] = t1;
	      t1 = sig_imag[i];
	      sig_imag[i] = sig_imag[j];
	      sig_imag[j] = t1;
	    }
	}


      n = 0;			/* FFT  */
      n2 = 1;

      for (i = 0; i < Nfft; i++)
	{
	  n = n2;
	  n2 = n2 + n2;
	  e = -6.283185307179586 / n2;
	  a = 0.0;

	  for (j = 0; j < n; j++)
	    {
	      c = cos (a);
	      s = sin (a);
	      a = a + e;

	      for (k = j; k < Signal_Length; k = k + n2)
		{
		  t1 = c * sig_real[k + n] - s * sig_imag[k + n];
		  t2 = s * sig_real[k + n] + c * sig_imag[k + n];
		  sig_real[k + n] = sig_real[k] - t1;
		  sig_imag[k + n] = sig_imag[k] - t2;
		  sig_real[k] = sig_real[k] + t1;
		  sig_imag[k] = sig_imag[k] + t2;
		}
	    }
	}
    }
  /*------------------------------------------------------------------*/
  /*        when the signal length is NOT a power of two              */
  /*            Calls the matlab subroutine fft                       */
  /*------------------------------------------------------------------*/
  else
    {
      int            num_out, num_in, i;
      mxArray       *outputArray[1];
      mxArray       *inputArray[1];
      mxArray       *array_ptr;


      num_out = 1;
      num_in = 1;

      /* recopy the real and imag parts of the signal in matrices */
      array_ptr = mxCreateDoubleMatrix (1, Signal_Length, mxCOMPLEX);
      memcpy (mxGetPr (array_ptr), sig_real, Signal_Length * sizeof (double));
      memcpy (mxGetPi (array_ptr), sig_imag, Signal_Length * sizeof (double));
      inputArray[0] = array_ptr;

      /* calls the MATLAB function */
      mexCallMATLAB (num_out, outputArray, num_in, inputArray, "fft");
      memcpy (sig_real, mxGetPr (outputArray[0]), Signal_Length * sizeof (double));

      if (mxIsComplex (outputArray[0]))
	{
	  memcpy (sig_imag, mxGetPi (outputArray[0]), Signal_Length * sizeof (double));
	}
      else
	{
	  for (i = 0; i < Signal_Length; i++)
	    sig_imag[i] = 0;
	}
      /* free memory */
      mxDestroyArray (outputArray[0]);
      mxDestroyArray (inputArray[0]);
    }
  return;
}

/*====================================================================*
 * Name of the function : po2 (int)                                   *
 * Author               :                                             *
 * Date of creation     :                                             *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * Computes the next power of to of an integer n, i.e.                *
 * the smallest k such that m=2^k>=n                                  *
 *====================================================================*/

int
po2 (int n)
{
  int            nextpo2, m;

  m = 1;
  nextpo2 = 0;
  while (m < n)
    {
      ++nextpo2;
      m = m + m;
    }

  return (nextpo2);
}

/*====================================================================*
 * Name of the function : sinc (double)                               *
 * Author               : Emmanuel ROY                                *
 * Date of creation     : 22 - 06 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * Computes the sinc function of a number, defined by :               *
 *                                                                    *
 *   for any x in R,                                                  *
 *                                                                    *
 *             /                                                      *
 *             | 1 if x = 0                                           *
 *             |                                                      *
 *   sinc(x) = | sin (pi * x)                                         *
 *             | ------------ if x != 0                               *
 *             |   pi * x                                             *
 *             \                                                      *
 *                                                                    *
 *====================================================================*/




double sinc (double x)
{
  double         resu;

  if (x == 0)
   	resu = 1.0;
  else
	resu = sin(pi*x)/(pi*x);

  return resu;
}


/*====================================================================*
 * Name of the function : irem (double x, double y)                   *
 * Author               : Emmanuel ROY                                *
 * Date of creation     : 22 - 06 -1999                               *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 *                                                                    *
 * Computes the remainder after Euclidean division of two doubles     *
 *                                                                    *
 *====================================================================*/

int irem( double x, double y)
{
 int result;

 if (y != 0)
   {
     result =  x-y* (int)(x/y);
   }
 else
   {
     result = 0;
     printf("irem.c : attempt to divide by 0\n");
   }

 return result;
}

/*====================================================================*
 * Name of the function : idx (int i_row, int j_col, int nb_row)      *
 * Author               : Emmanuel ROY                                *
 * Date of creation     :                                             *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 *                                                                    *
 * The matrices are stored as vectors, column by column               *
 * This program computes the vector index corresponding to the        *
 * specified line number (line), the column number (col) and the      *
 * number of lines   (nb_line)                                        *
 *====================================================================*/

int
idx (int i_row, int j_col, int nb_row)
{
  if (i_row >= nb_row)
    {
      printf("idx : incorrect row number\n");
      return 0;
      exit(0);
    }
  else
    {
      return (i_row + (j_col * nb_row));
    }
}

/*                                                		      *
 *====================================================================*
 * Name of the function :                			      *
 * Author               :                                             *
 * Date of creation     :                                             *
 *--------------------------------------------------------------------*
 * THE ALGORITHM             				              *
								      *
								      *
 *====================================================================*
 * INPUT VARIABLES   					              *
 * Name        | type  |              role               	      *
 *             |       | kernel 				      *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES    						      *
 * Name        | type  |              role                	      *
 *								      *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES 						      *
 * Name        | type  |              role                 	      *
 *								      *
 *====================================================================*
 * SUBROUTINES USED HERE				      	      *
 *--------------------------------------------------------------------*
 * Name   |                                          	     	      *
 * Action |							      *
 * Place  | 							      *
 *====================================================================*/

void
ifft (int Signal_Length, int Nfft, double *sig_real, double *sig_imag)
{
  /*------------------------------------------------------------------*/
  /*          when the signal length is a power of two                */
  /*------------------------------------------------------------------*/
  if (Signal_Length == (int) powof (2, Nfft) + 1)
    {
      int            i, j, k, n, n2;
      double         c, s, e, a, t1, t2;

      j = 0;			/* bit-reverse  */
      n2 = Signal_Length / 2;
      for (i = 1; i < Signal_Length - 1; i++)
	{
	  n = n2;
	  while (j >= n)
	    {
	      j = j - n;
	      n = n / 2;
	    }
	  j = j + n;

	  if (i < j)
	    {
	      t1 = sig_real[i];
	      sig_real[i] = sig_real[j];
	      sig_real[j] = t1;
	      t1 = sig_imag[i];
	      sig_imag[i] = sig_imag[j];
	      sig_imag[j] = t1;
	    }
	}


      n = 0;			/*IFFT  */
      n2 = 1;

      for (i = 0; i < Nfft; i++)
	{
	  n = n2;
	  n2 = n2 + n2;
	  e = 6.283185307179586 / n2;
	  a = 0.0;

	  for (j = 0; j < n; j++)
	    {
	      c = cos (a);
	      s = sin (a);
	      a = a + e;

	      for (k = j; k < Signal_Length; k = k + n2)
		{
		  t1 = c * sig_real[k + n] - s * sig_imag[k + n];
		  t2 = s * sig_real[k + n] + c * sig_imag[k + n];
		  sig_real[k + n] = sig_real[k] - t1;
		  sig_imag[k + n] = sig_imag[k] - t2;
		  sig_real[k] = sig_real[k] + t1;
		  sig_imag[k] = sig_imag[k] + t2;
		}
	    }
	}
      /* divide by Signal_Length */
      for (k = 0; k < Signal_Length; k++)
	{
	  sig_real[k] = sig_real[k] / Signal_Length;
	  sig_imag[k] = sig_imag[k] / Signal_Length;
	}
    }
  /*------------------------------------------------------------------*/
  /*        when the signal length is NOT a power of two              */
  /*            Calls the matlab subroutine ifft                      */
  /*------------------------------------------------------------------*/
  else
    {
      int            num_out, num_in, i;
      mxArray       *outputArray[1];
      mxArray       *inputArray[1];
      mxArray       *array_ptr;


      num_out = 1;
      num_in = 1;

      /* recopy the real and imag parts of the signal in matrices */
      array_ptr = mxCreateDoubleMatrix (1, Signal_Length, mxCOMPLEX);
      memcpy (mxGetPr (array_ptr), sig_real, Signal_Length * sizeof (double));
      memcpy (mxGetPi (array_ptr), sig_imag, Signal_Length * sizeof (double));
      inputArray[0] = array_ptr;

      /* calls the MATLAB function */
      mexCallMATLAB (num_out, outputArray, num_in, inputArray, "ifft");

      /* recovers the output */
      memcpy (sig_real, mxGetPr (outputArray[0]), Signal_Length * sizeof (double));
      if (mxIsComplex (outputArray[0]))
	{
	  memcpy (sig_imag, mxGetPi (outputArray[0]), Signal_Length * sizeof (double));
	}
      else
	{
	  for (i = 0; i < Signal_Length; i++)
	    sig_imag[i] = 0;
	}

      /* free memory */
      mxDestroyArray (outputArray[0]);
      mxDestroyArray (inputArray[0]);
    }
  return;
}



/*====================================================================*
 * Name of the function : sqr (double)                                *
 * Author               : Manuel DAVY - IRCYN                         *
 * Date of creation     : 02 - 02 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 *								      *
 * Computes the square value of x                                     *
 *====================================================================*/

double
sqr (double x)
{
  return (x * x);
}



/*====================================================================*
 * Name of the function : Recover_Signal (void)                       *
 * Author               : Manuel DAVY - IRCYN                         *
 * Date of creation     : 02 - 02 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * In Matlab environment, recovers a signal given by a signo.m program*
 * the function signo.m must be :                                     *
 * signal=signo(N,class,number);                                      *
 * N     = length of the signal                                       * 
 * class = class of the signal                                        *
 * number = number of the signal in its class                         *
 *====================================================================*/

void
Recover_Signal (int class, int signal_number,
		type_signal Signal)
{
  int            num_out, num_in;
  mxArray       *outputArray[1];
  mxArray       *inputArray[3];
  mxArray       *array_ptr;
  double        *inter_ptr;
  double         inter;


  /* arguments of the matlab function */
  num_out = 1;
  num_in = 3;

  inter = (double) Signal.length;
  array_ptr = mxCreateDoubleMatrix (1, 1, mxREAL);
  memcpy (mxGetPr (array_ptr), &inter, sizeof (double));
  inputArray[0] = array_ptr;


  inter = (double) class;
  array_ptr = mxCreateDoubleMatrix (1, 1, mxREAL);
  memcpy (mxGetPr (array_ptr), &inter, sizeof (double));
  inputArray[1] = array_ptr;

  inter = (double) signal_number;
  array_ptr = mxCreateDoubleMatrix (1, 1, mxREAL);
  memcpy (mxGetPr (array_ptr), &inter, sizeof (double));
  inputArray[2] = array_ptr;

  /* calls the MATLAB function */
  mexCallMATLAB (num_out, outputArray, num_in, inputArray, "signo3");

  /* recovers the output */
  inter_ptr = mxGetPr (outputArray[0]);
  memcpy (Signal.real_part, inter_ptr, Signal.length * sizeof (double));
  inter_ptr = mxGetPi (outputArray[0]);
  memcpy (Signal.imag_part, inter_ptr, Signal.length * sizeof (double));

  /* free memory */
  mxDestroyArray (outputArray[0]);
  mxDestroyArray (inputArray[0]);
  mxDestroyArray (inputArray[1]);
  mxDestroyArray (inputArray[2]);
}




 /*===================================================================*
 * Name of the function : transpose      			      *
 * Author               : Manuel Davy                                 *
 * Date of creation     : 10 - 02 - 1999                              *
 *--------------------------------------------------------------------*
 * THE ALGORITHM             				              *
 *								      *
 *								      *
 *====================================================================*
 * INPUT VARIABLES   					              *
 * Name        | type  |              role               	      *
 *             |       | kernel 				      *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES    						      *
 * Name        | type  |              role                	      *
 *								      *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES 						      *
 * Name        | type  |              role                 	      *
 *								      *
 *====================================================================*
 * SUBROUTINES USED HERE				      	      *
 *--------------------------------------------------------------------*
 * Name   |                                          	     	      *
 * Action |							      *
 * Place  | 							      *
 *====================================================================*/

void
transpose (int N_line, int N_col, double *matrix)
{
  int            line, col, index;
/* checks if the the matrix is not reduced to a single element */
  {
    if ((N_line > 1) && (N_col > 1))

      if (N_line == N_col)	/* the matrix is square */
	{
	  /* requires an intermediary element */
	  double         inter;
	  int            index_1, index_2;

	  for (line = 1; line < N_line; line++)
	    {
	      for (col = 0; col < line; col++)
		{
		  index_1 = idx (line, col, N_line);	/* in the under triangle */
		  index_2 = idx (col, line, N_col);	/* in the upper triangle */

		  inter = matrix[index_1];
		  matrix[index_1] = matrix[index_2];
		  matrix[index_2] = inter;
		}
	    }
	}
      else
	/* the matrix is not square */
	{
	  /* requires an intermediary matrix */
	  double        *inter;

	  inter = (double *) ALLOC (N_line * N_col, sizeof (double));


	  /* recopy in a transpose matrix */
	  for (line = 0; line < N_line; line++)
	    {
	      for (col = 0; col < N_col; col++)
		{
		  inter[idx (col, line, N_col)] = matrix[idx (line, col, N_line)];
		}
	    }
	  /* recopy in the original matrix */
	  for (index = 0; index < (N_line * N_col); index++)
	    {
	      matrix[index] = inter[index];

	    }

	  FREE (inter);
	}
  }
}

/*====================================================================*
 * Name of the function : fftshift                                    *
 * Date of creation     : 02 - 06 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * swaps the first and second halves of a vector. Example             *
 * [1 2 3 4 5 6 7 8 9 10 ] becomes [6 7 8 9 10 1 2 3 4 5]             *
 * The parameters to pass are :                          	      *
 *   - the input vector		                            	      *
 *   - the output vector	                            	      *
 *   - its length                                                     *
 * if the length is odd, example [1 2 3 4 5] becomes [4 5 1 2 3]      *
 *====================================================================*/
void
fftshift (double *vector_in, double *vector_out, int vector_length)

{
  double inter1, inter2;
  int i, half_length;


  /* computation of the half length in case of odd or even length */
  half_length = (int) (vector_length/2.0);


  /* case where the length is odd */
  if (ISODD(vector_length)==1)
    {
      inter2=vector_in[half_length];
      for (i=0; i<half_length; i++)
	{
	  inter1 = vector_in[i];
	  vector_out[i] = vector_in[half_length+i+1];
	  vector_out[half_length + i ] = inter1;
	}      
      vector_out[vector_length-1]=inter2;
    }
  /* case where the length is even */
  else
    {
      for (i=0; i<half_length; i++)
	{
	  inter1 = vector_in[half_length + i ];
	  vector_out[half_length + i] = vector_in[i];
	  vector_out[i] = inter1;
	}
    }
  /* fftshifting of the vector */
 


}

/*====================================================================*
 * Name of the function : mem_alloc_signal                            *
 * Date of creation     : 06 - 04 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * memory allocation for the type 'type_signal'                       *
 * the fields :                                                       *
 *   - length		                                 	      *
 *   - is_complex		                            	      *
 * must be previously initialized                                     *
 *====================================================================*/

void
mem_alloc_signal (type_signal *Signal, double *ptr_time_instant,
		  double *ptr_real_part, double *ptr_imag_part)
{
  /* some tests to make sure that all the fields are initialized */
  if ((*Signal).length <= 0)
    {
      printf ("mem_alloc_signal : Signal.length incorrect\n");
      exit(0);
    }


  if (((*Signal).is_complex != TRUE) && ((*Signal).is_complex != FALSE))
    {
      printf ("mem_alloc_signal : Signal.is_complex incorrect\n");
      exit(0);
    }


  /* memory allocation for the field 'time_instants' */
  if (ptr_time_instant == NULL)
    {
      (*Signal).time_instants = (double *) ALLOC ((*Signal).length, sizeof (double));
    }
  else
    {
      (*Signal).time_instants = ptr_time_instant;
    }

  if ((*Signal).time_instants == NULL)
    {
      printf ("mem_alloc_signal : memory allocation error\n");
      exit(0);
    }

  /* memory allocation for the field 'real_part' */

  if (ptr_real_part == NULL)
    {
      (*Signal).real_part = (double *) ALLOC ((*Signal).length, sizeof (double));
    }
  else
    {
      (*Signal).real_part = ptr_real_part;
    }

  if ((*Signal).real_part == NULL)
    {
      printf ("mem_alloc_signal : memory allocation error\n");
      exit(0);
    }

  /* memory allocation for the field 'imag_part' */

  if ((*Signal).is_complex == TRUE)
    {
      if (ptr_imag_part == NULL)
	{

	  (*Signal).imag_part = (double *) ALLOC ((*Signal).length, sizeof (double));
	}
      else
	{
	  (*Signal).imag_part = ptr_imag_part;
	}

      if ((*Signal).imag_part == NULL)
	{
	  printf ("mem_alloc_signal : memory allocation error\n");
	  exit(0);
	}
    }

}

/*====================================================================*
 * Name of the function : mem_free_signal                             *
 * Date of creation     : 06 - 04 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * free memory for variables of the type 'type_signal'                *
 *====================================================================*/

void
mem_free_signal (type_signal *Signal)
{

  FREE ((*Signal).time_instants);
  FREE ((*Signal).real_part);


  if ((*Signal).is_complex == TRUE)
    {
      FREE ((*Signal).imag_part);
    }
}


/*====================================================================*
 * Name of the function : mem_alloc_AF                                *
 * Date of creation     : 06 - 04 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * memory allocation for the type 'type_AF'                           *
 * the fields :                                                       *
 *   - N_doppler	                                 	      *
 *   - N_delay  	                                 	      *
 *   - is_complex		                            	      *
 * must be previously initialized                                     *
 *====================================================================*/

void
mem_alloc_AF (type_AF *ambig_func,
	      double *ptr_doppler_bins, double *ptr_delay_bins,
	      double *ptr_real_part, double *ptr_imag_part)
{

  /* some tests to make sure that all the fields are initialized */
  if ((*ambig_func).N_doppler <= 0)
    {
      printf ("mem_alloc_AF : AF.N_doppler incorrect\n");
      exit(0);
    }

  if ((*ambig_func).N_delay <= 0)
    {
      printf ("mem_alloc_AF : AF.N_delay incorrect\n");
      exit(0);
    }

  if (((*ambig_func).is_complex != TRUE) && ((*ambig_func).is_complex != FALSE))
    {
      printf ("mem_alloc_AF : AF.is_complex incorrect\n");
      exit(0);
    }


  /* memory allocation for the field 'doppler_bins' */
  if (ptr_doppler_bins == NULL)
    {
      (*ambig_func).doppler_bins = (double *) ALLOC
            ((*ambig_func).N_doppler * (*ambig_func).N_delay, sizeof (double));

    }
  else
    {
      (*ambig_func).doppler_bins = ptr_doppler_bins;
    }

  if ((*ambig_func).doppler_bins == NULL)
    {
      printf ("mem_alloc_AF : memory allocation error\n");
      exit(0);
    }

  /* memory allocation for the field 'delay_bins' */
  if (ptr_delay_bins == NULL)
    {
      (*ambig_func).delay_bins = (double *) ALLOC ((*ambig_func).N_doppler
						    * (*ambig_func).N_delay,
						    sizeof (double));

    }
  else
    {
      (*ambig_func).delay_bins = ptr_delay_bins;
    }

  if ((*ambig_func).delay_bins == NULL)
    {
      printf ("mem_alloc_AF : memory allocation error\n");
      exit(0);
    }


  /* memory allocation for the field 'real_part' */

  if (ptr_real_part == NULL)
    {
      (*ambig_func).real_part = (double *) ALLOC ((*ambig_func).N_doppler *
						   (*ambig_func).N_delay,
						   sizeof (double));

    }
  else
    {
      (*ambig_func).real_part = ptr_real_part;
    }

  if ((*ambig_func).real_part == NULL)
    {
      printf ("mem_alloc_AF : memory allocation error\n");
      exit(0);
    }

  /* memory allocation for the field 'imag_part' */
  if ((*ambig_func).is_complex == TRUE)
    {
      if (ptr_imag_part == NULL)
	{
	  (*ambig_func).imag_part = (double *) ALLOC((*ambig_func).N_doppler
						       *
						       (*ambig_func).N_delay,
						       sizeof (double));

	}
      else
	{
	  (*ambig_func).imag_part = ptr_imag_part;
	}

      if ((*ambig_func).imag_part == NULL)
	{
	  printf ("mem_alloc_AF : memory allocation error\n");
	  exit(0);
	}
    }
}

/*====================================================================*
 * Name of the function : mem_free_AF                                 *
 * Date of creation     : 06 - 04 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * free memory for variables of the type 'type_AF'                    *
 *====================================================================*/

void
mem_free_AF (type_AF *ambig_func)
{
  FREE ((*ambig_func).doppler_bins);
  FREE ((*ambig_func).delay_bins);
  FREE ((*ambig_func).real_part);


  if ((*ambig_func).is_complex == TRUE)
    {
      FREE ((*ambig_func).imag_part);
    }
}


/*====================================================================*
 * Name of the function : mem_alloc_TFR                               *
 * Date of creation     : 06 - 04 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * memory allocation for the type 'type_TFR'                          *
 * the fields :                                                       *
 *   - N_freq   	                                 	      *
 *   - N_time   	                                 	      *
 *   - is_complex		                            	      *
 * must be previously initialized                                     *
 *====================================================================*/

void
mem_alloc_TFR (type_TFR *tfr,
	       double *ptr_freq_bins, double *ptr_time_instants,
	       double *ptr_real_part, double *ptr_imag_part)
{

  /* some tests to make sure that all the fields are initialized */
  if ((*tfr).N_freq <= 0)
    {
      printf ("mem_alloc_TFR : TFR.N_freq incorrect\n");
      exit(0);
    }

  if ((*tfr).N_time <= 0)
    {
      printf ("mem_alloc_TFR : TFR.N_time incorrect\n");
      exit(0);
    }

  if (((*tfr).is_complex != TRUE) && ((*tfr).is_complex != FALSE))
    {
      printf ("mem_alloc_TFR : TFR.is_complex incorrect\n");
      exit(0);
    }


  /* memory allocation for the field 'freq_bins' */
  if (ptr_freq_bins == NULL)
    {
      (*tfr).freq_bins = (double *) ALLOC ((*tfr).N_freq * (*tfr).N_time,
					    sizeof (double));
    }
  else
    {
      (*tfr).freq_bins = ptr_freq_bins;
    }

  if ((*tfr).freq_bins == NULL)
    {
      printf ("mem_alloc_TFR : memory allocation error\n");
      exit(0);
    }

  /* memory allocation for the field 'time_instants' */
  if (ptr_time_instants == NULL)
    {
      (*tfr).time_instants = (double *) ALLOC ((*tfr).N_freq *
						(*tfr).N_time, sizeof (double));
    }
  else
    {
      (*tfr).time_instants = ptr_time_instants;
    }

  if ((*tfr).time_instants == NULL)
    {
      printf ("mem_alloc_TFR : memory allocation error\n");
      exit(0);
    }


  /* memory allocation for the field 'real_part' */

  if (ptr_real_part == NULL)
    {
      (*tfr).real_part = (double *) ALLOC ((*tfr).N_freq * (*tfr).N_time,
					    sizeof (double));
    }
  else
    {
      (*tfr).real_part = ptr_real_part;
    }

  if ((*tfr).real_part == NULL)
    {
      printf ("mem_alloc_TFR : memory allocation error\n");
      exit(0);
    }

  /* memory allocation for the field 'imag_part' */
  if ((*tfr).is_complex == TRUE)
    {
      if (ptr_imag_part == NULL)
	{
	  (*tfr).imag_part = (double *) ALLOC ((*tfr).N_freq *
						(*tfr).N_time, sizeof (double));
	}
      else
	{
	  (*tfr).imag_part = ptr_imag_part;
	}

      if ((*tfr).imag_part == NULL)
	{
	  printf ("mem_alloc_TFR : memory allocation error\n");
	  exit(0);
	}
    }
}

/*====================================================================*
 * Name of the function : mem_free_TFR                                *
 * Date of creation     : 06 - 04 - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 * free memory for variables of the type 'type_TFR'                   *
 *====================================================================*/

void
mem_free_TFR (type_TFR *tfr)
{
  FREE ((*tfr).freq_bins);
  FREE ((*tfr).time_instants);
  FREE ((*tfr).real_part);


  if ((*tfr).is_complex == TRUE)
    {
      FREE ((*tfr).imag_part);
    }
}
