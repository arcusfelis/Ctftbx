/* EXISTS AN INTERFACE PROGRAM TO MATLAB : KERMEX.C		      *
 *====================================================================*
 * Name of the function : kernel (void)				      *
 * Author               : Manuel DAVY - IRCYN			      *
 * Date of creation     : 20 - 01 - 1999			      *
 *--------------------------------------------------------------------*
 * THE ALGORITHM          					      *
 * This function generates various kernel given a shape and the params*
 * The kernel is computed in the AMBIGUITY PLANE (doppler-delay). The *
 * kernel is a matrix of given size. The kernels are :		      *
 * - MTEK  : multiform tiltable exponential kernel 		      *
 * - RGK   : Radially gaussian kernel				      *
 * - GMCWK : Generalized marginals Choi-Williams kernel		      *
 * - WV    : Wigner-Ville kernel				      *
 * See the references given in the comments below.		      *
 *								      *
 * Dans le cadre de l'analyse temps-frequence, cette fonction genere  *
 * des noyaux de Representations temps-frequence (RTF) dans le plan   *
 * des ambiguites (retard - doppler). Le noyau est mis sous la forme  *
 * d'une matrice de dimensions passees a la fonction. Chaque noyau est*
 * defini par un type et des parametres.                              *
 *								      *
 *====================================================================*
 * INPUT VARIABLES     						      *
 * Name           |                   role                  	      *
 * N_Doppler      | Number of doppler bins i.e. number of rows in the *
 *                | kernel matrix                            	      *
 * N_Delay        | Number of delays bins i.e. number of columns in   *
 *                | the kernel matrix                        	      *
 * type           | Kernel type (see the definitions in 'tftb.h')     *
 * parameters     | Parameter vector for the kernel      )    	      *
 * nb_param       | Number of parameters                     	      *
 * ker.N_doppler  | Number of lines in the kernel matrix              *
 * ker.N_delay    | Number of columns in the kernel matrix            *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES                                                   *
 * Name           |                   role                  	      *
 * ker            | Matrix containing the computed kernel    	      *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES						      *
 * Name           |                   role                  	      *
 * line,col       | Line and columns index in the matrix              * 
 * doppler        | Value of the doppler (resp. delay) in the current *
 * delay          |  line (resp. column)  			      *
 * inter          | Intermediary value for a computation	      *
 *====================================================================*/



void
kernel (int type, double *parameters, int nb_param, type_AF ker)
{
  int            line, col;
  double         doppler, delay;
  double         inter;


  /* some tests */

  if ((ker.N_doppler <1) || (ker.N_delay < 1))
    {
      printf("kernel.c : invalid number of lines / columns in the kernel matrix \n");
      exit(0);
    }

  switch (type)
    {
      /***************************************************************
       *              Multiform Tiltable Exponentiel Kernel          *
       ***************************************************************
       * The parameters vector is made of                            *
       * alpha, beta, gamma, r, tau0, nu0, lambda]                   *
       *   see the reference :                                       * 
       *	H. Costa and G.F. Boudreaux-Bartels,                 *
       *	Design of Time-Frequency Representations Using a     *
       *		  Multiform, Tiltable Exponential Kernel     *
       *	IEEE Trans. on Signal Processing                     *
       *	October 1995, vol. 43, no 10, pp 2283-2301           *
       ***************************************************************/
      /* LOCAL VARIABLES                                             *
       * Name           |                   role                     *
       * A              | intermediary in the computation of the     *
       *                | MTEK kernel                                *
       ***************************************************************/
   case MTEK:
      {
	/* local variables for the MTEK */

	double         A;

	/* Kernel Construction */
	for (line = 0; line < ker.N_doppler; line++)
	  {
	    doppler = (line - ker.N_doppler / 2.0 + 1.0) / ker.N_doppler;
	    for (col = 0; col < ker.N_delay; col++)
	      {
		delay = col - ker.N_delay / 2.0 + 1.0;
		A = (doppler * delay) / (TAU_0 * NU_0);

		/* case of the symmetrical kernel */
		if ((BETA == 2) && (GAMMA == 0.5))
		  {
		    A = ABS (A);
		  }

		/* case where the marginals do not have to be verified */
		if (ALPHA == 0)
		  {
		    inter = sqr (delay / TAU_0) + sqr (doppler / NU_0)
		      + 2.0
		      * R * A;
		  }
		/* case where the marginals have to be verified */
		else
		  {
		    inter = sqr (delay / TAU_0) * powof (sqr (doppler /
							      NU_0), ALPHA)
		      + sqr (doppler / NU_0) * powof (sqr (delay /
							   TAU_0), ALPHA)
		      + 2.0 * R * A;
		  }
		/* test to avoid the computation of log(0) */
		if (inter == 0)
		  {
		    ker.real_part[idx (line, col, ker.N_doppler)] = 1.0;
		  }
		else
		  {
		    ker.real_part[idx (line, col, ker.N_doppler)] =
		      exp (-pi * powof (sqr (inter), LAMBDA));
		  }
		inter = 0;
		A = 0;
	      }
	  }
      }
      break;
      /***************************************************************
       *                 Radially Gaussian kernel                    *
       ***************************************************************
       * The parameters vector is made according to the rule :       *
       * if the order of the kernel is p, the vector is              *
       * [ c ,a1 , ... , ap, b1,... ,bp] where c is the constant     *
       * ai are the cosine coefficients and bi the sine              *
       * coefficients in the Fourier series decomposition of the     *
       * contour. See the reference :                                *
       *     M. Davy and C. Doncarli,                                *
       *    Optimal Kernels of Time-Frequency Representations for    *
       *    Signal Classification,                                   *
       *	  TFTS 1998, pp 581-584.                             *
       ***************************************************************/
      /* LOCAL VARIABLES                                             *
       * Name           |                   role                     *
       * order          | Fixes the maximum p in the vector of params*
       * p              | Current parameter p in the vector of params*
       * a,b            | Vectors containing the parameters ap, bp   *
       *                | from the vector of params                  *
       * c              | The minimum value of the contour function  *
       * inter, mini    | Intermediary values in the computations    *
       * phi            | angle parameter in the polar coordinates   *
       * rho2           | square radius parameter in polar coord.    *         
       ***************************************************************/


    case RGK:
      {
	/* local variables for the RGK */
	int            order, p;
	double         *a, *b;
	double         c, inter, mini;
	double         phi, rho2;


	/* some error cases to avoid ... */
	if (ISODD(nb_param) == 0)
	  {
	    printf("kernel.c : the number of RGK parameters must be ODD\n");
	    exit(0);
	  } 


	order = (nb_param - 1) / 2;


	/* memory allocation for a and b */
	a = (double *) ALLOC ( order , sizeof(double) );
	b = (double *) ALLOC ( order , sizeof(double) );

	/* variables recovery */
	c = parameters[0];
	for (p = 0; p < order; p++)
	  {
	    a[p] = parameters[p + 1];
	    b[p] = parameters[order + p + 1];
	  }
	/*-----------------------------------------------*/
	/*             Kernel Construction               */
	/*-----------------------------------------------*/

	/* minimum value of the contour function */
	mini = 0;	


	/* construction of the matrix of the contour function */
	for (line = 0; line < ker.N_doppler; line++)
	  {
	    /* current doppler value */
	    doppler = (line - ker.N_doppler / 2.0 + 1.0) / ker.N_doppler;
	    /* normalization of the doppler to have angles in radians */
	    doppler = doppler * sqrt (ker.N_delay);

	    for (col = 0; col < ker.N_delay; col++)
	      {
		/* currrent delay value */
		delay = col - ker.N_delay / 2.0 + 1.0;
		/* normalization of the delay to have angles in radians */
		delay = delay / sqrt (ker.N_delay);

		/* computation of the angles in the ambiguity plane */
		if (((delay > 0) && (doppler > 0))
		    || ((delay < 0) && (doppler < 0)))
		  {
		    phi = atan (doppler / delay);
		  }
		else
		  {
		    phi = atan (doppler / delay) + pi;
		  }

		inter = 0;
		for (p = 0; p < order; p++)
		  {
		    inter = inter + a[p] * cos (2.0 * (p + 1) * phi) +
		      b[p]
		      * sin (2.0 * (p + 1) * phi);
		  }
		/* look for the minimum */
		if (inter < mini)
		  {
		    mini = inter;
		  }
		/* matrix of the contour function : each element in this */
		/* matrix contains the value of the contour function for */
		/* the corresponding delay and doppler values */
		ker.real_part[idx (line, col, ker.N_doppler)] = inter;

	      }
	  }

	/* construction of the RGK matrix */
	for (line = 0; line < ker.N_doppler; line++)
	  {
	    /* current normalized doppler */
	    doppler = (line - ker.N_doppler / 2.0 + 1.0) / ker.N_doppler;
	    doppler = doppler * sqrt (ker.N_delay);

	    for (col = 0; col < ker.N_delay; col++)
	      {
		/* current normalized delay */
		delay = col - ker.N_delay / 2.0 + 1.0;
		delay = delay / sqrt (ker.N_delay);
		/* Square Polar radius */
		rho2 = sqr (doppler) + sqr (delay);
		ker.real_part[idx (line, col, ker.N_doppler)] =
		  exp (-(rho2 / 2.0) / 
		       sqr (ker.real_part[idx (line, col, ker.N_doppler)]
			    - mini + c));

		/* case of the center of the ambiguity plane */
		if ((delay == 0) && (doppler == 0))
		  ker.real_part[idx (line, col, ker.N_doppler)] = 1.0;
	      }
	  }
	/* free the memory used here */
	FREE (a);
	FREE (b);
      }
      break;

      /***************************************************************
       *          Generalized Marginals Choi-Williams                *
       ***************************************************************
       * The parameters vector is made of the width of the           *
       * branches "sigma" and the angle of each branch theta_i       *
       * (0 to Pi). The number of branches is given by the number    *
       * of coefficients                                             *
       * [ sigma theta_1 theta_2  ... theta_p]                       *
       *                                                             *
       * see the reference :                                         *
       * X.-G. Xia and Y. Owechko and B. H. Soffer and R. M. Matic,  *
       * Generalized-Marginal Time-Frequency Distributions,          *
       * TFTS 1996, pp. 509-51                                       *
       ***************************************************************/
      /* LOCAL VARIABLES                                             *
       * Name           |                   role                     *
       * N_Branch       | Number of kernel branches                  *
       * branch         | Current branch number                      *
       * angle          | vector of branches angles                  *
       * sigma          | value of the parameter sigma (kernel width)*
       * inter          | intermediary variable                      *
       ***************************************************************/


    case GMCWK:
      {
	/* local variables for the RGK */
	int            N_branch, branch;
	double         *angle, sigma, inter;



	/* some error cases to avoid ... */
	if (nb_param < 2)
	  {
	    printf("kernel.c : at least 2 GMCWK parameters required\n");
	    exit(0);
	  } 


	/* variables recovery */
	sigma = parameters[0];
	N_branch = nb_param - 1;



	/* recovery of the angles */
	angle = &(parameters[1]);

	/* Kernel Construction  */
	for (line = 0; line < ker.N_doppler; line++)
	  {
	    doppler = (line - ker.N_doppler / 2.0 + 1.0) 
	      / ker.N_doppler;
	    doppler = doppler * sqrt (ker.N_delay);

	    for (col = 0; col < ker.N_delay; col++)
	      {
		delay = col - ker.N_delay / 2.0 + 1.0;
		delay = delay / sqrt (ker.N_delay);
		inter = 1;
		for (branch = 0; branch < N_branch; branch++)
		  {
		    inter = inter 
		      * sqr (doppler * cos (angle[branch])
			     + delay * sin (angle[branch]));
		  }
		ker.real_part[idx (line, col, ker.N_doppler)] = 
		  exp (-inter / sigma);
	      }
	  }
      }
      break;
      /***************************************************************
       *                Wigner-Ville kernel                          *
       ***************************************************************
       * No parameter is required. The matrix is equal to            *
       * one in each point                                           *
       ***************************************************************/


    case WIGNER:
      {



	/* Kernel Construction  */
	for (line = 0; line < ker.N_doppler; line++)
	  {
	    for (col = 0; col < ker.N_delay; col++)
	      {
		ker.real_part[idx (line, col, ker.N_doppler)] = 1.0;
	      }
	  }
      }
      break;
      /***************************************************************
       *                   Spectrogram  kernel                       *
       ***************************************************************
       *  The paremters are the window points in the time            *
       * domain. The kernel is the ambiguity function of the window  *
       ***************************************************************/
      /* LOCAL VARIABLES                                             *
       * Name           |                   role                     *
       * window         | structure containing the window of the     *
       *                | spectrogram                                *
       * AF_h           | Ambiguity function of th window            *
       * min_delay |    | in order to pad with zeros, defines the    *
       * max_delay |    | indices where the zero padding begins/ends *
       * index          | indice to access the elements in matrices  *
       *                | (stored as vectors)                        *
       ***************************************************************/



    case SPECTRO:
      {
	/* local variable */
	type_signal    window;
	type_AF        AF_h;
	int            min_delay, max_delay, index;

	/* some error cases to avoid ... */
	if (ISODD(nb_param) == 1)
	  {
	    printf("kernel.c : the window length must be EVEN for SPECTRO\n");
	    exit(0);
	  } 



	/* initialization of the local variables : window */
	window.length = ker.N_delay;
	window.is_complex = TRUE;
	/* creation of memory */
	window.real_part = (double *) ALLOC (window.length, sizeof (double));
	window.imag_part = (double *) ALLOC (window.length, sizeof (double));


	min_delay = (int) (window.length / 2 - nb_param / 2);
	max_delay = (int) (window.length / 2 + nb_param / 2 - 1);

	/* initialization of the imaginary part for the window 
	   and zero padding out of [min_delay;max_delay] */
	for (col = 0; col < window.length; col++)
	  {
	    if ((col >= min_delay) && (col <= max_delay))
	      {
		window.real_part[col] = parameters[col - min_delay];
	      }
	    else
	      {
		window.real_part[col] = 0;
	      }
	    window.imag_part[col] = 0;
	  }


	/* initialization of the local variables : AF_h */
	AF_h.N_doppler = ker.N_doppler;
	AF_h.N_delay = ker.N_delay;
	AF_h.is_complex = TRUE;

	mem_alloc_AF (&AF_h, NULL, NULL, NULL, NULL);
 	for (index = 0; index < AF_h.N_delay; index++)
	  {
	    AF_h.delay_bins[index] = -AF_h.N_delay / 2.0 + 1.0 + index;
	  }




	/* computation of the AF of the window */
	af (window, AF_h);

	/* Warning : if ker.N_delay != ker.N_doppler,
	   one must add or suppress lines in this AF */
	for (line = 0; line < ker.N_doppler; line++)
	  {
	    for (col = 0; col < ker.N_delay; col++)
	      {
		index = idx (line, col, ker.N_delay);
		ker.real_part[idx (line, col, ker.N_doppler)] =
		  sqrt (sqr (AF_h.real_part[index])
			+ sqr (AF_h.imag_part[index]));
	      }
	  }
	/* free memory !! */

	mem_free_AF (&AF_h);
	FREE (window.real_part);
	FREE (window.imag_part);

      }
      break;
    }
}
