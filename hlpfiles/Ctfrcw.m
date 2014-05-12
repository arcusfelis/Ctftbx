%CTFRCW   Choi-Williams time-frequency distribution.
% [TFR,T,F]=TFRCW(X,T,N,G,H,SIGMA)
%
% computes the Choi-Williams distribution of a signal X. 
%  
% X     = Analyzed signal.
% T     = time instant(s)          (default : 1:length(X)).
% N     = number of frequency bins (default : length(X)).
% G     = time smoothing window,   (default : Hamming(N/10)). 
% H     = frequency smoothing window, in the time-domain,
%                                  (default : Hamming(N/4)). 
% SIGMA = kernel width
%
% TFR   = time-frequency representation.
% F     = vector of normalized frequencies.
%
% Example :
%   x = hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128)); 
%   g=Cwindow(9,'Hamming'); h=Cwindow(27,'Hamming'); 
%   t=1:128; [tfr,T,F]=Ctfrcw(x,t,128,g,h,1);
%   imagesc(T,F,tfr); axis xy; xlabel('time'); ylabel('frequency')
%
% SEE ALSO Ctfrwv, Ctfrsp

