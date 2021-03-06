%CTFRMH	 Margenau-Hill time-frequency distribution.
% [TFR,T,F]=TFRMH(X,T,N)
%
% computes the Margenau-Hill distribution of a  signal X. 
%  
% X     = Analyzed signal.
% T     = time instant(s)          (default : 1:length(X)).
% N     = number of frequency bins (default : length(X)).
%
% TFR   = time-frequency representation.
% F     = vector of normalized frequencies.
%
% Example :
%   x = hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128));  
%   t=1:128; [tfr,T,F]=Ctfrmh(x,t,128);
%   imagesc(T,F,tfr); axis xy; xlabel('time'); ylabel('frequency')
%
% SEE ALSO : Ctfrpmh, Ctfrwv, Ctfrri
