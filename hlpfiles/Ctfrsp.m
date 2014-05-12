%CTFRSP Spectrogram Time Frequency distribution.
% [SP,T,F,norm] = Ctfrsp(X,T,N,H)
%
% Computes the  spectrogram of a signal X
%
% X        = Analyzed signal
% T        = the time instant(s)   (default : 1:length(X)).
% N        = number of frequency bins (default : length(X)).
% H        = frequency smoothing window, (default : Hamming(N/4)).
%
% SP       = Spectrogram
% F         = Vector of frequency bins
% NORM      = Vector of normalization applied at each time instant
%
% Example :
%    x = hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128));
%    [sp,T,F] = Ctfrsp(x,1:128,128);
%    imagesc(T,F,sp); axis xy
%
% SEE ALSO : Ctfrrsp, Ctfrstft, Ctfrwv, Cwindow
