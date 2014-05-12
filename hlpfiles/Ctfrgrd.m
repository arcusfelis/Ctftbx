%CTFRGRD Generalized rectangular time-frequency distribution.
% [TFR,T,F]=TFRGRD(X,T,N,G,H,RS,MOVERN)
%
% computes the Generalized Rectangular distribution
% of a signal X. 
%
% X      = Analyzedd signal
% T      = time instant(s)         (default : 1:length(X)).
% N      = number of frequency bins (default : length(X))
% G      = time smoothing window    (default : Hamming(N/10)). 
% H      = frequency smoothing window (default : Hamming(N/4)). 
% RS     = kernel width            (default : 1).
% MOVERN = dissymmetry ratio       (default : 1).
%
% TFR    = time-frequency representation
% F      = vector of normalized frequencies.
% 
% Example :
%     x = hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128));  
% 	 g=Cwindow(9,'hamming'); h=window(27,'hamming'); 
% 	 t=1:128; [tfr,T,F]=Ctfrgrd(x,t,128,g,h,36,1/5);
%        imagesc(T,F,tfr); axis xy; xlabel('time'); ylabel('frequency')
