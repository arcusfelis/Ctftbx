%CTFRPMH Pseudo Margenau-Hill time-frequency distribution.
% [TFR,T,F]=TFRPMH(X,T,N,H)
%
% computes the Pseudo Margenau-Hill distribution of a  signal X. 
%  
% X     = Analyzed signal.
% T     = time instant(s)          (default : 1:length(X)).
% N     = number of frequency bins (default : length(X)).
% H     = frequency smoothing window (default : Hamming(N/4)). 
%
% TFR   = time-frequency representation.
% F     = vector of normalized frequencies.
%
% Example :
%   x = hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128));  
%   t=1:128; [tfr,T,F]=Ctfrpmh(x,t,128);
%   imagesc(T,F,tfr); axis xy; xlabel('time'); ylabel('frequency')
%
% SEE ALSO : Ctfrmh, Ctfrpwv
