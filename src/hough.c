/* EXISTS AN INTERFACE PROGRAM TO MATLAB : ---MEX.c                   *
 *====================================================================*
 * Name of the function :                                             *
 * Author               :                                             *
 * Date of creation     :                                             *
 *--------------------------------------------------------------------*
 * THE ALGORITHM                                                      *
 *                                                                    *
 *====================================================================*
 * INPUT VARIABLES                                                    *
 * Name           |                   role                            *
 *                |                                                   *
 *                |                                                   *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES                                                   *
 * Name           |                   role                            *
 *                |                                                   *
 *                |                                                   *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES                                                 *
 * Name           |                   role                            *
 *                |                                                   *
 *                |                                                   *
 *====================================================================*
 * SUBROUTINES USED HERE                                              *
 *--------------------------------------------------------------------*
 * Name   |                                                           *
 * Action |                                                           *
 * Place  |                                                           *
 *====================================================================*/

/* entete reduit pour les sous programmes courts */
/*====================================================================*
 * Name of the function :                                             *
 * Author               : Manuel DAVY - IRCYN                         *
 * Date of creation     :    -    - 1999                              *
 *--------------------------------------------------------------------*
 * Action of the function                                             *
 *                                                                    *
 *								      *
 *====================================================================*/
void
hough (type_TFR tfr, double nb_theta,double  nb_rho,
       double* transfo_hough, double* rho_vect, double* theta_vect)
{
  int      time,freq, index, index2;
  int      row_min, row_mid, row_max, row;
  int      column_min, column_mid, column_max, column;

  double   rho_max, step_rho, step_theta;
  double   inter, max_tfr;
  double   rho,theta;


  /* the tfr matrix must be real valued */
  if(tfr.is_complex == TRUE)
    {
      printf ("hough.c : the input tfr must be real-valued \n");
      exit(0);
    }

  rho_max = sqrt(sqr(tfr.N_freq)+sqr(tfr.N_time))/2.0;
  step_rho = rho_max / (nb_rho-1.0);
  step_theta = 2.0*pi / nb_theta;

  /* construction of the vectors rho_vect */
  for(index=0;index<nb_rho;index++)
    {
      rho_vect[index] = step_rho * index;
    }
 /* construction of the vectors theta_vect */
  for(index=0;index<nb_theta;index++)
    {
      theta_vect[index] = step_theta * index;
    }



  /* research of the maximum value of the tfr */
  inter = -1.0* 1e100;
  for(time=0;time<tfr.N_time;time++)
    {
      for(freq=0;freq<tfr.N_freq;freq++)
	{
	  inter = MAX(tfr.real_part[idx(freq,time,tfr.N_freq)],inter);
	}
    }
  max_tfr = inter;

  /* initialization of the output matrix */
  for(index=0;index<nb_theta*nb_rho;index++)
    {
      transfo_hough[index] = 0.0;
    }

  /* determines the min and max bounds for rho et theta */
  /* rows (rho) */
  if((tfr.N_freq)%2 == 1) /* the number of rows is odd */
    {
      row_mid = (tfr.N_freq+1)/2;
      row_min = 1-row_mid;
      row_max = row_mid-1;
    }
  else
    {
     row_mid = tfr.N_freq/2;
     row_min = 1-row_mid;
     row_max = row_mid;
    }
  /* columns (theta)*/
 if((tfr.N_time)%2 == 1) /* the number of columns is odd */
    {
      column_mid = (tfr.N_time+1)/2;
      column_min = 1-column_mid;
      column_max = column_mid-1;
    }
  else
    {
     column_mid = tfr.N_time/2;
     column_min = 1-column_mid;
     column_max = column_mid;
    }
 /* computation of the hough transform */
 for(row=row_min;row<=row_max;row++)
   for(column=column_min;column<=column_max;column++)
     for(index=0; index<nb_theta; index++)
       {
	 theta=theta_vect[index];
	 rho = row*cos(theta)-column*sin(theta);
	 if((rho>0)&&(rho<=rho_max))
	   {
	     index2=idx(ROUND(rho/step_rho),index,nb_rho);
	     transfo_hough[index2]=transfo_hough[index2]+
	       tfr.real_part[idx(row+row_mid-1,column+column_mid-1,tfr.N_freq)];
	   }
       }
}
