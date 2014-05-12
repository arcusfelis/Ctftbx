%CTFRRIDB  Reduced Interference Distribution with Bessel kernel.
% [TFR,T,F]=TFRRIDB(X,T,N,G,H)
%
% Computes the Reduced Interference Distribution with a kernel
%  based on the Bessel function of the first kind, for a signal X. 
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
%   t=1:128; [tfr,T,F]=Ctfrridb(x,t,128,g,h);
%   imagesc(T,F,tfr); axis xy; xlabel('time'); ylabel('frequency')
%
% SEE ALSO : CTRRIDH, CTFRRIDT, CTFRRIDBN





