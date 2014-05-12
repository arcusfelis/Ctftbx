#%CTFRKER Time frequency representation kernel 
% kernel = Ctfrker(N_Doppler,N_Delay,kernel_type,optional parameters)
%
% Computes the Cohen's group Time Frequency Representation's (TFR)
% kernel, in the ambiguity plane.
%
% N_doppler   = Number of doppler bins 
%             = number of rows in the kernel matrix
% N_delay     = Number of delay bins
%             = number of columns in the kernel matrix
% kernel      = Computed kernel
%
% kernel_type = one of the following types:
%         
%  MTEK    :  Multiform tiltable exponential kernel
%      7 parameters required :
%      [alpha beta gamma r tau_0 nu_0 lambda]
%  RGK     : Radially gaussian kernel :
%      an odd number of parameters required, the
%      Fourier descriptors of the contour function
%      [c a1 ... ap b1 ... bp]
%  GMCWK   : Generalized marginals Choi-Williams kernel
%      at least two parameters required : the kernel width
%      and the angles of the branches in [0 pi].
%      [sigma angle_1 .. angle_2]
%  WV      : Wigner-Ville kernel : no parameters required
%  SPECTRO : Spectrogram kernel
%      The parameters are the window with EVEN length
%
% note:  more kernels will be implemented in future versions
%
% Example :
%   AF = Cambifunb(hilbert(sin(2*pi*0.25*(1:128))),-63:64,128);
%   ker = Ctfrker(128,128,'rgk',[1 1 1]);
%   TFR = Caf2tfr(AF,ker);
%   imagesc(1:128,linspace(0,0.5,128),TFR); axis xy;
%
% SEE ALSO : Cambifunb, Caf2tfr