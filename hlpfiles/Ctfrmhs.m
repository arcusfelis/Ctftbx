%CTFRMHS Margenau-Hill-Spectrogram time-frequency distribution.
% [TFR,T,F]=TFRMHS(X,T,N,G,H)
%
% computes the Margenau-Hill-Spectrogram distribution of a  signal X. 
%  
% X     = Analyzed signal.
% T     = time instant(s)          (default : 1:length(X)).
% N     = number of frequency bins (default : length(X)).
% G     = time smoothing window,   (default : Hamming(N/10)). 
% H     = frequency smoothing window, in the time-domain,
%                                  (default : Hamming(N/4)). 
%
% TFR   = time-frequency representation.
% F     = vector of normalized frequencies.
%
% Example :
%   x = hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128));
%   g=Cwindow(21,'Hamming');	 h=Cwindow(63,'Hamming'); 
%   t=1:128; [tfr,T,F]=Ctfrmhs(x,t,128,g,h);
%   imagesc(T,F,tfr); axis xy; xlabel('time'); ylabel('frequency')
%
% SEE ALSO : Ctfrmh, Ctfrpmh, Ctfrsp