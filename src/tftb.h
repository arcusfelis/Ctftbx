#include <math.h>
#include <stdio.h>
#include "mex.h"

#define UNIX               1
#define WINDOWS            2


#include "system.h"




/* WARNING WARNING WARNING WARNING WARNING WARNING */
/* WARNING WARNING WARNING WARNING WARNING WARNING */
/* WARNING WARNING WARNING WARNING WARNING WARNING */
/* WARNING WARNING WARNING WARNING WARNING WARNING */
/* WARNING WARNING WARNING WARNING WARNING WARNING */

/*------------------------------------------------*/
/*          YOU SHOULD NOT EDIT THIS FILE         */
/*------------------------------------------------*/

/* WARNING WARNING WARNING WARNING WARNING WARNING */
/* WARNING WARNING WARNING WARNING WARNING WARNING */
/* WARNING WARNING WARNING WARNING WARNING WARNING */
/* WARNING WARNING WARNING WARNING WARNING WARNING */
/* WARNING WARNING WARNING WARNING WARNING WARNING */


#if SYSTEME == UNIX
   #define ALLOC(A,B)     calloc((A),(B))
   #define FREE(A)        cfree((A)) 
#else
  #define ALLOC(A,B)      malloc((A)*(B))
  #define FREE(A)         free((A)) 
#endif

/*   local functions   */
#define MAX(A, B)       ((A) > (B) ? (A) : (B))
#define MIN(A, B)       ((A) < (B) ? (A) : (B))
#define sgn(A)          ((A) > 0.0 ? 1.0 : -1.0)
#define ABS(a)          (((a) >= (0)) ? (a) : (-a))
#define SWAP(a,b)       {temp = (a); (a)=(b); (b)=temp;}
#define ROUND1(x)       (((((x)-(int)(x))>=0) &&(((x)-(int)(x))<0.5)) ? ((int)(x)) : ((int)(x+1)))
#define ROUND(x)        ((int)(sgn((x))*ROUND1((x)*sgn((x)))))
#define ISODD(x)        ((x/2.0)== ((int)(x/2)) ? 0 : 1)


/* local constants */
#define pi               3.141592653589793
#define EPS              0.0000000001

#define TRUE             1
#define FALSE            0

/* definition of the window identifiers */
#define RECTANG          1
#define HAMMING          2
#define HANNING          3
#define KAISER           4
#define NUTTALL          5
#define BLACKMAN         6
#define HARRIS           7
#define BARTLETT         8
#define TRIANG           8
#define BARTHANN         9
#define PAPOULIS         10
#define GAUSS            11
#define PARZEN           12
#define HANNA            13
#define DOLPH            14
#define DOLF             14
#define NUTBESS          15
#define SPLINE           16


/* definition of the distance identifiers */
#define LQ               1
#define QUADRATIC        2
#define CORRELATION      3
#define KOLMOGOROV       4
#define KULLBACK         5
#define CHERNOFF         6
#define MATUSITA         7
#define NLQ              8
#define LSD              9
#define JENSEN           10


/* definition of the kernel shapes */
#define MTEK             1
#define RGK              2
#define GMCWK            3
#define WIGNER           4
#define SPECTRO          5


/* parametres for the MTEK */
#define NB_PARAM_MTEK    7
#define ALPHA            parameters[0]
#define BETA             parameters[1]
#define GAMMA            parameters[2]
#define R                parameters[3]
#define TAU_0            parameters[4]
#define NU_0             parameters[5]
#define LAMBDA           parameters[6]

/*-------------------------------------------*/
/*  definition of the structures and types   */
/*-------------------------------------------*/

/* Signal structure */
typedef struct SIG
  {
    int            length;	/* Length of the signal in points */
    double         sample_freq;	/* Sample frequency of the signal */
    double        *time_instants; /* instants of sample for the signal */
    unsigned char  is_complex;	/* TRUE if there exists an imag part */
    double        *real_part;	/* real part of the signal */
    double        *imag_part;	/* imaginary part of the signal */
  }
type_signal;

typedef struct Time_freq_rep
  {
    int            N_freq;	/* number of freq bins in the TFR matrix */
    int            N_time;	/* number of time_bins in the TFR matrix */
    double        *freq_bins;	/* fqs for each line of the matrix */
    double        *time_instants; /* instant for each column of the TFR */
    unsigned char  is_complex;	/* TRUE if there exists an imag part */
    double        *real_part;	/* real part of the TFR */
    double        *imag_part;	/* imaginary part of the TFR */
  }
type_TFR;

typedef struct Ambi_func
  {
    int            N_doppler;	/* number of doppler bins in the AF */
    int            N_delay;	/* number of delay bins in the AF matrix */
    double        *doppler_bins; /* doppler bin for each line of the AF */
    double        *delay_bins;	/* delay bin for each column of the AF */
    unsigned char  is_complex;	/* TRUE if there exists an imag part */
    double        *real_part;	/* real part of the AF */
    double        *imag_part;	/* imaginary part of the AF */
  }
type_AF;
