/* EXISTS AN INTERFACE PROGRAM TO MATLAB : reasSP.C		      *
 *====================================================================*
 * Name of the function : reas_spectro   			      *
 * Author               : Manuel DAVY                                 *
 * Date of creation     : 15 - 03 - 1999                              *
 *--------------------------------------------------------------------*
 * THE ALGORITHM             				              *
 *								      *
 * Given a signal, a smoothing window and the time instants where to  *
 * compute, the reassigned and the classical spectrogram are computed.*
 * Furthermore, the reassignment field is given as an output. It is   *
 * possible, for certain applications, to weigh the field of          *
 * reassignment with a weighting_field.				      *
 *								      *
 * The metod consists in computing the 3 spectrograms of the signal   *
 * with the window, its derivative and the time times the window.     *
 *  STFT_H, STFT_DH, STFT_TH					      *
 *								      *
 *								      *
 * Then, the reassignment field is computed according to the classical*
 * formulas (in each time-frequency point) : 			      *
 *								      *
 *								      *
 *	       -              -					      *
 *	       |     STFT_TH  | 				      *
 *	       |  Re -------  | 				      *
 *	       |     STFT_H   |					      *
 *    Field =  | 	      |					      *
 *	       |     STFT_DH  |					      *
 *	       |  Im -------  |					      *
 *	       |     STFT_H   |					      *
 *	       -	      -					      *
 *								      *
 * Rmq : The subroutine 'stft' computes the Short time Fourier Transf.*
 *	 With a ponderation that depends on wether the window covers  *
 *	 only points of the signal or zeros padding points. It is     *
 *	 necessary to suppress this normalization here, because it    *
 *	 depends on the window (and then is different for e.g. the    *
 *	 derivative of the window)				      *
 *								      *
 *								      *
 * The additional weighting field is optional.			      *
 *====================================================================*
 * INPUT VARIABLES   					              *
 * Name              |                role               	      *
 * Signal            | Signal for which the reassigned spectro is     *
 *                   | computed                                       *
 * Signal.length     | Length of the input signal		      *
 * Signal.is_complex | Indicates wether the signal is complex valued  *
 *                   | or real_valued                                 *
 * Signal.real_part  | Real and part of the signal	              *
 * Signal.imag_part  | Imaginary part of the signal (only if          *
 *                   | 'Signal.is_complex == TRUE)                    *
 *                   |                                                *
 * Window            | The window 'h' in the STFT, with ODD length    *
 * Window_Length     | Length of the window			      *
 *                   |                                                *
 * weighting_field   | Optional weigthing field (matrix size equals   *
 *                   | the size of the reassigned spectrogram         *
 *use_weigthing_field| TRUE if the preceding field is employed        *
 * TFR_reassigned :  | tfr that is reassigned after the computation   *
 *   .N_time         | Number of positions for the smoothing window   *
 *                   | = Number of columns in the final reas. spectro *
 *   .N_freq         | Number of frequency bins in the reas. spectro  *
 *                   | = Number of rows in the final reas. spectro    *
 *   .freq_bins      | MUST BE ALLOCATED                              * 
 *   .time_instants  | time instants where the smoothing window is    *
 *                   | positioned (length = TFR_reassigned.N_time)    *
 *                   | its points must be REGULARLY spaced            *
 *   .is_complex     | must be set to FALSE here (the spectro is real *
 *                   | valued)                                        *
 * TFR_not_reassigned| final spectrogram (not reassigned)             *
 *   .N_time         | Number of positions for the smoothing window   *
 *                   | = Number of columns in the final spectrogram   *
 *   .N_freq         | Number of frequency bins in the spectrogram    *
 *                   | = Number of rows in the final spectrogram      * 
 *   .freq_bins      | MUST BE ALLOCATED                              * 
 *   .time_instants  | time instants where the smoothing window is    *
 *                   | positioned (length = TFR_not_reassigned.N_time)*
 *                   | its points must be REGULARLY spaced            *
 *   .is_complex     | must be set to FALSE here (the spectro is real *
 *                   | valued)                                        *
 *--------------------------------------------------------------------*
 * OUTPUT VARIABLES    						      *
 * Name              |                  role                	      *
 * TFR_reassigned :  | tfr that is reassigned after the computation   *
 *  .freq_bins       | vector of frequency bins                       *
 *  .real_part       | matrix of the reassigned spectrogram           *
 *                   |                                                *
 * TFR_not_reassigned| final spectrogram (not reassigned)             *
 *  .freq_bins       | vector of frequency bins                       *
 *  .real_part       | matrix of the reassigned spectrogram           *
 *                   |                                                *
 * field_time        |   Field of reassignement employed, decomposed  *
 * field_freq        |   on the basis vectors parallel to the time    *
 *                   |   axis and the freq. axis		      *
 *--------------------------------------------------------------------*
 * INTERNAL VARIABLES 						      *
 * Name              |                role                 	      *
 * normh             |   Vectors of the normalization used at each    *
 * normth            |   time instant in the computation of a STFT    *
 * normdt            |   size : TFR_reassigned.N_time		      *
 *                   |             				      *
 * twindow           |   The window times a vector of time	      *
 * dwindow           |   The derivative of the window		      *
 *                   |             				      *
 * stft_h            | the STFT of the signal computed with the       *
 * and its fields    | window 'window'			              * 
 *                   |             				      *
 * stft_th           | the STFT of the signal computed with the       *
 * and its fields    | window 'twindow'			              * 
 *                   |             				      *
 * stft_dh           | the STFT of the signal computed with the       *
 * and its fields    | window 'dwindow'			              * 
 *                   |             				      *
 * time, freq        |   current location in the TFR matrices         *
 * index             |   index of the preceding location in the vector*
 *                   |   containing the matrix.			      *
 *                   |             				      *
 * module            |   value of the spectrogram at the position     *
 *                   |   (time,freq)				      *
 *                   |             				      *
 * factor            | intermediary variable corresponding to a factor*
 *                   | of normalization of the reassignment field     *
 * step              | used to correct the derivative of the window at*
 *                   | its edges                                      *
 * step_time         | displacement step of the window                *
 * flag              | boolean variable used when a test concerning   *
 *                   | the input variables fails.
 *====================================================================*
 * SUBROUTINES USED HERE				      	      *
 *--------------------------------------------------------------------*
 * Name   | int idx(int line, int row, int nb_row)                    *
 * Action | computes the vector index for an element in a matrix given* 
 *        | the line and column indices and the number of columns     *
 * Place  | divers.c                                                  *
 *--------------------------------------------------------------------*
 * Name   | void stft (double *Sig_real, double *Sig_imag,            *
 *        |            int Signal_Length, double *Time_inst,          *
 *        |            int Time_inst_Length, int Nfft, double *Window,*
 *        |            int Window_Length, double *stft_real,          *
 *        |            double *stft_imag, double *norm_vector)        *
 * Action | Computes the Short time Fourier transform of a signal     *
 * Place  | stft.c						      *
 *--------------------------------------------------------------------*
 * Name   | void gradient (double *matrix, int size_x, int size_y,    *
 *        |             double step_x, double step_y, double *grad_x, * 
 *        |             double *grad_y)                               *
 * Action | Computes the gradient of a matrix			      *
 * Place  | gradient.c						      *
 *--------------------------------------------------------------------*
 * Nom    | double sqr(x)   				              *
 * Action | computes the square of x				      *
 * Lieu   | divers.c						      *
 *--------------------------------------------------------------------*
 * Name   | void reassign (double *TFR_to_reassign, double *field_x,  *
 *        |                double *field_y, int N_time, int N_freq,   *
 *        |	           double *TFR_reassigned)                    *
 * Action | moves the pixels in a TFR according to a field of reassig.*
 * Place  | reassign.c						      *
 *====================================================================*/

void
reas_spectro (type_signal Signal,
	      double *Window, int Window_Length,
	      type_TFR TFR_reassigned, type_TFR TFR_not_reassigned,
	      double *field_time, double *field_freq)

{
  /*-----------------------------------------------------------------*/
  /*                           VARIABLES                             */
  /*-----------------------------------------------------------------*/
  double        *normh, *normth, *normdh;
  double        *twindow, *dwindow;
  type_TFR       stft_h, stft_th, stft_dh;
  int            time, freq, index;
  double         module, step;
  double         step_time;
  unsigned char  flag;
  /*-----------------------------------------------------------------*/
  /*                       SOME TESTS ...                            */
  /*-----------------------------------------------------------------*/

 /*                         TFR_reassigned                           */

   if (TFR_reassigned.is_complex == TRUE)
    {
      printf ("reas_spectro.c : The reassigned tfr matrix must be real valued\n");
      exit(0);
    }
   if (TFR_reassigned.N_time <0)
    {
      printf ("reas_spectro.c : The field TFR_reassigned.N_time is not correctly set\n");
      exit(0);
    }
   if (TFR_reassigned.N_freq <0)
    {
      printf ("reas_spectro.c : The field TFR_reassigned.N_freq is not correctly set\n");
      exit(0);
    }
  if (TFR_reassigned.time_instants[0] <0)
    {
      printf ("reas_spectro.c : The field TFR_reassigned.time_instants is not correctly set\n");
      exit(0);
    }
 /*                       TFR_not_reassigned                         */

   if (TFR_not_reassigned.is_complex == TRUE)
    {
      printf ("reas_spectro.c : The tfr matrix must be real valued\n");
      exit(0);
    }
   if (TFR_not_reassigned.N_time <0)
    {
      printf ("reas_spectro.c : The field TFR_not_reassigned.N_time is not correctly set\n");
      exit(0);
    }
   if (TFR_not_reassigned.N_freq <0)
    {
      printf ("reas_spectro.c : The field TFR_not_reassigned.N_freq is not correctly set\n");
      exit(0);
    }

  if (TFR_not_reassigned.time_instants[0] <0)
    {
      printf ("reas_spectro.c : The field TFR_not_reassigned.time_instants is not correctly set\n");
      exit(0);
    }
 /*checks the compatibility between TRF_reassigned ans TFR_not_reassigned*/

   if ((TFR_not_reassigned.N_time != TFR_reassigned.N_time ) ||
       (TFR_not_reassigned.N_freq != TFR_reassigned.N_freq ) )
    {
      printf ("reas_spectro.c : incompatible fields in TFR_reassigned and TFR_not_reassigned\n");
    }

   flag=FALSE;
   for (time=0 ; time <TFR_reassigned.N_time ; time++)
     {
       if (TFR_not_reassigned.time_instants[time]
                 != 
            TFR_reassigned.time_instants[time])
	 {
	   flag=TRUE;
	 }
     }

   if (flag == TRUE)
     {
       printf ("reas_spectro.c : incompatible fields in TFR_reassigned and TFR_not_reassigned\n");
       exit(0);
     }

  step_time = TFR_reassigned.time_instants[1] - TFR_reassigned.time_instants[0];

   flag=FALSE;
   for (time=0 ; time <(TFR_reassigned.N_time-1) ; time++)
     {
       if ((TFR_reassigned.time_instants[time+1]-
            TFR_reassigned.time_instants[time])
           != step_time)
	 {
	   flag=TRUE;
	 }
     }

  if (flag == TRUE)
     {
       printf ("reas_spectro.c : time instants not regularly spaced\n");
       exit(0);
     }

 /*                   checks that the window length is odd           */

  if (ISODD(Window_Length) == 0)
    {
      printf ("stft.c : The window Length must be an ODD number\n");
      exit(0);
    }

 /*                     tests the weighing fiels                     */



 /*--------------------------------------------------------------------*/
 /*           creation of the vector of frequency bins  (output)       */
 /*--------------------------------------------------------------------*/
  for (time = 0; time < TFR_reassigned.N_freq; time++)
    {
      TFR_reassigned.freq_bins[time] = (double) time / TFR_reassigned.N_freq;
    }
  for (time = 0; time < TFR_not_reassigned.N_freq; time++)
    {
      TFR_not_reassigned.freq_bins[time] = (double) time / TFR_not_reassigned.N_freq;
    }

  /*-----------------------------------------------------------------*/
  /*                        MEMORY ALLOCATION                        */
  /*-----------------------------------------------------------------*/
  /* allocation of memory for the norm vectors */
  normh = (double *) ALLOC (TFR_reassigned.N_time, sizeof (double));
  normth = (double *) ALLOC (TFR_reassigned.N_time, sizeof (double));
  normdh = (double *) ALLOC (TFR_reassigned.N_time, sizeof (double));

  /* allocation of memory for the windows */
  twindow = (double *) ALLOC (Window_Length, sizeof (double));
  dwindow = (double *) ALLOC (Window_Length, sizeof (double));



  /*-----------------------------------------------------------------*/
  /*                   INITIALIZATION of the OTHER FIELDS            */
  /*-----------------------------------------------------------------*/
  stft_h.is_complex = TRUE;
  stft_th.is_complex = TRUE;
  stft_dh.is_complex = TRUE;
  stft_h.N_time = TFR_reassigned.N_time;
  stft_h.N_freq = TFR_reassigned.N_freq;
  stft_th.N_time = TFR_reassigned.N_time;
  stft_th.N_freq = TFR_reassigned.N_freq;
  stft_dh.N_time = TFR_reassigned.N_time;
  stft_dh.N_freq = TFR_reassigned.N_freq;

  /* allocation of memory for the STFTs matrices */

  /* in order to save memory, 'stft_h.real_part' is stored in the same table 
     as 'TFR_not_reassigned.real_part' */

  mem_alloc_TFR (&stft_h, NULL, TFR_reassigned.time_instants,
		 TFR_not_reassigned.real_part, NULL);
  /* in order to save memory, 'stft_th.real_part is stored in the same table 
     as 'field_time' */
  mem_alloc_TFR (&stft_th, NULL, TFR_reassigned.time_instants,
		 field_time, NULL);
  /* in order to save memory, 'stft_dh.real_part'is stored in the same table 
     as 'field_freq' */
  mem_alloc_TFR (&stft_dh, NULL, TFR_reassigned.time_instants,
		 field_freq, NULL);

  /*-----------------------------------------------------------------*/
  /*                       COMPUTATION                               */
  /*-----------------------------------------------------------------*/
  step_time = TFR_reassigned.time_instants[1] - TFR_reassigned.time_instants[0];

  /* computation of the STFT of the signal with 'window' */

  stft (Signal, Window, Window_Length, stft_h, normh);

  /* dwindow is the derivative of window */
  /* (here twindow is used as an compulsary but unusefull output) */
  gradient (Window, Window_Length, 1, 1, 1, dwindow, twindow);
  step = (Window[0] + Window[Window_Length - 1]) / 2.0;
  dwindow[0] = dwindow[0] + step;
  dwindow[Window_Length - 1] = dwindow[Window_Length - 1] - step;

  /* twindow is the window multiplied by a vector of time */
  /* -(Window_Length/2-1) .... 0 ... Window_Length/2) */
  for (time = 0; time < Window_Length; time++)
    {
      twindow[time] = Window[time] * (time - (Window_Length - 1.0) / 2.0);
    }


  /* computation of the STFT of the signal with 'twindow' */
  stft (Signal, twindow, Window_Length, stft_th, normth);

  /* computation of the STFT of the signal with 'dwindow' */
  stft (Signal, dwindow, Window_Length, stft_dh, normdh);

  /* computation of the reassignement fields */
  /* in order to save memory, the first field is stored in stft_th_real
     and the other in stft_dh_real */

  for (time = 0; time < TFR_reassigned.N_time; time++)
    {
      for (freq = 0; freq < TFR_reassigned.N_freq; freq++)
	{
	  index = idx (freq, time, TFR_reassigned.N_freq);

	  /* module = spectrogram (not reassigned) at the place 'index' */
	  module = sqr (stft_h.real_part[index]) + sqr (stft_h.imag_part[index]);
	  if  (module > EPS)
	    {
	      /* the first field equals REAL(STFT_TH / STFT_H) */
	      field_time[index] = (stft_th.real_part[index] *
				   stft_h.real_part[index]
				   + stft_th.imag_part[index] *
				   stft_h.imag_part[index])
		/ module;

	      /* the second field equals - IMAG(STFT_DH / STFT_H) */
	      field_freq[index] = -(stft_dh.imag_part[index] *
				    stft_h.real_part[index]
				    - stft_dh.real_part[index] *
				    stft_h.imag_part[index])
		/ module;


	      /* normalization of the fields of reassignement */

	      field_time[index] = field_time[index]
		* ( normh[time]/ normth[time]) / step_time;

	      field_freq[index] = field_freq[index]
		* (normh[time] / normdh[time])
		* (0.5*TFR_reassigned.N_freq / pi);

	    }
	  else
	    /* when the spectrogram is nearly null : no reassignement */
	    {
	      field_time[index] = 0.0;
	      field_freq[index] = 0.0;
	    }
	
	  /* the original (not reassigned) spectrogram is stored in stft_th_imag */
	  TFR_not_reassigned.real_part[index] = module;
 
	}
    }
  
  /* Reassignement of the spectrogram according to the field of vectors */
   reassign (TFR_not_reassigned, field_time, field_freq, TFR_reassigned); 
  

  /*-----------------------------------------------------------------*/
  /*                         FREE MEMORY                             */
  /*-----------------------------------------------------------------*/
 
   FREE (dwindow); 
   FREE (twindow); 
   FREE (normdh); 
   FREE (normth);
   FREE (normh);
   FREE (stft_dh.imag_part);
   FREE (stft_th.imag_part);
   FREE (stft_h.imag_part);
   FREE (stft_dh.freq_bins);
   FREE (stft_th.freq_bins);
   FREE (stft_h.freq_bins);
 
}
