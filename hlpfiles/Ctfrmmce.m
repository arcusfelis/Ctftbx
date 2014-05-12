%CTFRMMCE Minimum mean cross-entropy combination of spectrograms.
% [TFR,T,F]=TFRMMCE(X,H,T,N)
%
% computes  the minimum mean cross-entropy combination of
% spectrograms using aswindows the columns of the matrix H
%  
% X     = Analyzed signal.
% H     = frequency smoothing windows, stored in a matrix.
%         Each row is a window
% T     = time instant(s)          (default : 1:length(X)).
% N     = number of frequency bins (default : length(X))
%
% TFR   = time-frequency representation.
% F     = vector of normalized frequencies.
%
% Example :
%   x = hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128));
%   h(10+(-5:5),1)=Cwindow(11,'gauss'); h(10+(-7:7),2)=Cwindow(15,'gauss');  
%   h(10+(-9:9),3)=Cwindow(19,'gauss'); [tfr,T,F]=Ctfrmmce(x,h);
%   imagesc(T,F,tfr); axis xy; xlabel('time'); ylabel('frequency')
%
% SEE ALSO : Ctfrsp