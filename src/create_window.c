/*====================================================================*
 * Name of the function : window.c                                    *
 * Author               : Manuel DAVY                                 *
 * Date of creation     : 03 - 06 - 1999                              *
 *--------------------------------------------------------------------*
 * THE ALGORITHM                                                      *
 *                                                                    *
 * creates a window of given shapes and length.                       *
 * possible shapes and name :                                         *
 *     Rectangular -> RECTANG                                         *
 *     Hamming     -> HAMMING                                         *
 *     Hanning     -> HANNING                                         *
 *     Kaiser      -> KAISER            (1 optional parameter)        *
 *     Nuttal      -> NUTTALL                                          *
 *     Blackman    -> BLACKMAN                                       *
 *     Harris      -> HARRIS                                          *
 *     Triangular  -> BARTLETT, TRIANG                                *
 *     Barthann    -> BARTHANN                                        *
 *     Papoulis    -> PAPOULIS                                        *
 *     Gauss       -> GAUSS             (1 optional parameter)        *
 *     Parzen      -> PARZEN                                          *
 *     Hanna       -> HANNA             (1 optional parameter)        *
 *     Dolph (Dolf)-> DOLPH, DOLF                                     *
 *     Nutbess     -> NUTBESS           (2 optional parameters)       *
 *     Spline      -> SPLINE            (1 compulsary and 1 optional  *
 *                                       parameter)                   *
 *====================================================================*
 * INPUT VARIABLES                                                    *
 * Name           |                   role                            *
 * Window_type    | Type of window to compute (one of the choices     *
 *                | above)                                            *
 * Window_length  | Length of the window                              *
 * nb_param       | Number of parameters passed to compute the window *
 *                | used only for certain distances  (0<=nb_param<=2) *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES                                                   *
 * Name           |                   role                            *
 * i              | index in the window                               *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES                                                 *
 * Name           |                   role                            *
 * window         | vector of length 'Window_length' containing the   *
 *                | computed window                                   *
 *====================================================================*
 * SUBROUTINES USED HERE                                              *
 *--------------------------------------------------------------------*
 * Name   |                                                           *
 * Action |                                                           *
 * Place  |                                                           *
 *====================================================================*/

void 
create_window (int Window_type, int Window_length, double* param,
	int nb_param, double *Window)
{
  /* variables */
  int i;
  

  /* tests the input variable */
  if (Window_length <= 0)
    {
      printf("create_window.c : Bad window length\n");
      exit(0);
    }

  if ((nb_param != 0) && (nb_param != 1) && (nb_param != 2))
     {
      printf("create_window.c : Bad number of parameters\n");
      exit(0);
    }



  /* computation according to the window type */
  switch (Window_type)
    {
      /* -----------------------------------------------------*/
    case RECTANG:
      /* -----------------------------------------------------*/
      if (nb_param >0)
	{
	  printf("create_window.c : no param. required for RECTANG window\n");
	}
      for (i = 0; i < Window_length; i++)
	{
	  Window[i] = 1;
	}
      break;

      /* -----------------------------------------------------*/
    case HAMMING:
      /* -----------------------------------------------------*/
      if (nb_param >0)
	{
	  printf("create_window.c : no param. required for HAMMING window\n");
	}

   for (i = 0; i < Window_length; i++)
	{
	  Window[i] = 0.54 - 0.46 * cos ((2.0 * pi * (i + 1.0)) /
				     (Window_length + 1.0));
	}
      break;
      /* -----------------------------------------------------*/
    case HANNING:
      /* -----------------------------------------------------*/
      if (nb_param >0)
	{
	  printf("create_window.c : no param. required for HANNING window\n");
	}

      for (i = 0; i < Window_length; i++)
	{
	  Window[i] = 0.50 - 0.50 * cos ((2.0 * pi * (i + 1.0)) /
				     (Window_length + 1.0));
	}
      break;

      /* -----------------------------------------------------*/
    case KAISER:
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double beta;
      

	if (nb_param >1)
	  {
	    printf("create_window.c : maximum one param. required for KAISER window\n");
	  }
	
	if (nb_param == 1)
	  {
	    beta =  param[0];
	  }
	else
	  {
	    beta= 3.0 * pi ;
	  }
	for (i = 0; i < Window_length; i++)
	  {
	    Window[i] = 0;
	  }
	printf("create_window : Window KAISER not implemented yet\n");
      }
      break;

      /* -----------------------------------------------------*/
    case NUTTALL:
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double   indice;

	if (nb_param >0)
	  {
	    printf("create_window.c : no param. required for NUTTAL window\n");
	  }
	
	for (i = 0; i < Window_length; i++)
	  {
	    indice = ((-(Window_length - 1.0)/2.0 + i) * 2.0 * pi) / Window_length;
	    Window[i] = 0.3635819 + 
	      0.4891775*cos(indice) + 
	      0.1363995*cos(2.0*indice) +
	      0.0106411*cos(3.0*indice) ;
	  }
      }
      break;

      /* -----------------------------------------------------*/
    case BLACKMAN:
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double   indice;

	if (nb_param >0)
	  {
	    printf("create_window.c : no param. required for BLACKMANN window\n");
	  }
	
	for (i = 0; i < Window_length; i++)
	  {
	    indice = ((-(Window_length - 1.0)/2.0 + i) * 2.0 * pi) / Window_length;
	    Window[i] = 0.42 + 0.5 * cos(indice) +
	      0.08 * cos (2.0 * indice);
	  }
      }
    break;

      /* -----------------------------------------------------*/
    case HARRIS:
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double   indice;

	if (nb_param >0)
	  {
	    printf("create_window.c : no param. required for HARRIS window\n");
	  }

	for (i = 0; i < Window_length; i++)
	  {
	    indice = (2.0 * pi * (i + 1.0)) /(Window_length + 1.0);
	    Window[i] = 0.35875   - 
	      0.48829 * cos(indice) +
	      0.14128 * cos(2.0*indice) -
	      0.01168 *cos(3.0*indice);
	  }
      }
      break;

      /* -----------------------------------------------------*/
    case BARTLETT: /* or case TRIANG */
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double   indice;

	if (nb_param >0)
	  {
	    printf("create_window.c : no param. required for BARTLETT window\n");
	  }
	
	for (i = 0; i < Window_length; i++)
	  {
	    indice = MIN(i + 1.0 , Window_length - i);
	    Window[i] = (2.0 * indice) / (Window_length + 1.0);
	  }
      }
      break;


      /* -----------------------------------------------------*/
    case BARTHANN:
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double   indice;
      
	if (nb_param >0)
	  {
	    printf("create_window.c : no param. required for BARTHANN window\n");
	  }

	for (i = 0; i < Window_length; i++)
	  {
	    indice = MIN(i + 1.0 , Window_length - i);
	    Window[i] = 0.38 * (1.0 - cos ((2.0 * pi * (i + 1.0)) /
					   (Window_length + 1.0)))
	      +
	      0.48 * indice / (Window_length + 1.0);
	  }
      }
      break;

      /* -----------------------------------------------------*/
    case PAPOULIS:
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double   indice;
      
	if (nb_param >0)
	  {
	    printf("create_window.c : no param. required for PAPOULIS window\n");
	  }

	for (i = 0; i < Window_length; i++)
	  {
	    indice = ((i + 1.0) * pi) / (Window_length + 1.0);
	    Window[i] = sin (indice);
	  }
      }
    break;

      /* -----------------------------------------------------*/
    case GAUSS:
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double   indice, K;

	if (nb_param >1)
	  {
	    printf("create_window.c : maximum one param. required for GAUSS window\n");
	  }
     
	if (nb_param == 1)
	  {
	    K = param[0];
	  }
	else
	  {
	    K = 0.005;
	  }
	
	for (i = 0; i < Window_length; i++)
	  {
	    indice = -1.0 + (2.0*i) / (Window_length - 1.0);
	    Window[i] = exp((indice * indice) * log(K));
	  }
      }
    break;

      /* -----------------------------------------------------*/
    case PARZEN:
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double   indice , temp;

	if (nb_param >0)
	  {
	    printf("create_window.c : no param. required for PARZEN window\n");
	  }
	
	for (i = 0; i < Window_length; i++)
	  {
	    indice = ABS(((-(Window_length - 1.0)/2.0 + i)) * 2.0 / Window_length);
	    temp = 2.0 * powof(1.0 - indice,3);
	    Window[i] = MIN(temp - (1.0 - 2.0*indice) * (1.0 - 2.0*indice) * (1.0 - 2.0*indice) , temp );
	  }
      }
    break;

      /* -----------------------------------------------------*/
    case HANNA:
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double   L;
      
	if (nb_param >1)
	  {
	    printf("create_window.c : maximum one param. required for HANNA window\n");
	  }
	
	if (nb_param == 1)
	  {
	    L = param[0];
	  }
	else
	  {
	    L = 1;
	  }
	
	for (i = 0; i < Window_length; i++)
	  {
	    Window[i] = powof(
			      sin(((2.0 * i + 1.0) * pi) / (2.0 * Window_length)),
			      2.0 * L);
	  }
      }
    break;

      /* -----------------------------------------------------*/
    case DOLPH: /* or case DOLF */
      /* -----------------------------------------------------*/
      /* local variable */
      {     
	printf("create_window : window DOLPH not implemented yet\n");
      }
	break;


      /* -----------------------------------------------------*/
    case NUTBESS:
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double   indice, beta, nu;

	beta= 3.0 * pi;
	nu = 0.5 ;

	if (nb_param >2)
	  {
	    printf("create_window.c : maximum two param. required for NUTBESS window\n");
	  }
     
	if (nb_param == 1)
	  {
	    beta = param[0];
	  }
	if (nb_param == 2)
	  {
	    nu = param[1];
	  }
	
	for (i = 0; i < Window_length; i++)
	  {
	    indice = ((-(Window_length - 1.0) + i) * pi) / Window_length;
	    Window[i] = 0;
	  }
	printf("create_window : window NUTBESS not implemented yet\n");
      }
    break;


    /* -----------------------------------------------------*/
    case SPLINE:
      /* -----------------------------------------------------*/
      /* local variable */
      {
	double     nfreq,p;
	double indice, inter;
	
	if ((nb_param != 1) && (nb_param != 2))
	  {
	  printf("create_window.c : One/two parameter required for SPLINE window\n");
	  exit(0);
	  }
	else
	  {
	    nfreq = param[0];
	  }
	if (nb_param == 2)
	  {
	    p = param[1];
	  }
	else
	  {
	    p = pi * Window_length * nfreq / 10.0;
	  }
	
	
	for (i = 0; i < Window_length; i++)
	  {
	    indice = -(Window_length - 1.0)/2.0 + i;
	    inter = (0.5 * nfreq / p) * indice ;

	    if (inter != 0.0)
	      {
		inter = sin(inter * pi)/(inter * pi);
		Window[i] = powof(inter,p);
	      }
	    else
	      {
		Window[i] = 1.0;
	      }
	  }
      }
    break;


    default :
      {
      printf("create_window.c : Unknowm window type\n");
      exit(0);
      }
    break;

    }


}
