%CTFRRIDH  Reduced Interference Distribution with Hanning kernel.
% [TFR,T,F]=TFRRIDH(X,T,N,G,H)
%
% Computes the Reduced Interference Distribution with a kernel
%  based on the Hanning window, for a signal X. 
%  
% X     = Analyzed signal.
% T     = time instant(s)          (default : 1:length(X)).
% N     = number of frequency bins (default : length(X)).
% G     = time smoothing window,   (default : Hamming(N/10)). 
% H     = frequency smoothing window (default : Hamming(N/4)). 
%
% TFR   = time-frequency representation.
% F     = vector of normalized frequencies.
%
% Example :
%   x = hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128));  
%   g = Cwindow(31,'rect'); h=Cwindow(63,'rect'); 
%   t=1:128; [tfr,T,F]=Ctfrridh(x,t,128,g,h);
%   imagesc(T,F,tfr); axis xy; xlabel('time'); ylabel('frequency')
%
% SEE ALSO : CTRRIDB, CTFRRIDT, CTFRRIDBN





