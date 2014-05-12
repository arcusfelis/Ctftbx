%CTFRPWV	Pseudo Wigner-Ville time-frequency distribution.
% [TFR,T,F]=TFRPWV(X,T,N,H)
%
% computes the Pseudo Wigner-Ville distribution of a  signal X. 
%  
% X     = Analyzed signal.
% T     = time instant(s)          (default : 1:length(X)).
% N     = number of frequency bins (default : length(X)).
% H     = frequency smoothing window, in the time-domain,
%                                  (default : Hamming(N/4)). 
%
% TFR   = time-frequency representation.
% F     = vector of normalized frequencies.
%
% Example :
%   x = hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128));  
%   t=1:128; [tfr,T,F]=Ctfrpwv(x,t,128);
%   imagesc(T,F,tfr); axis xy; xlabel('time'); ylabel('frequency')
%
% SEE ALSO : Ctfrwv, Ctfrspwv
