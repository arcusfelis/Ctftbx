%CTFRSTFT Short time Fourier transform
% [STFT,T,F,NORM] = Ctfrstft(X,T,N,H);
%
% Computes the Short Time Fourier Transform (STFT) of Signal
%
% X         = signal for which the STFT is computed
% T         = time instant(s)          (default : 1:length(X)).
% N         = number of frequency bins (default : length(X)).
% H         = frequency smoothing window (default : Hamming(N/4)). 
%
% STFT      = Computed Short time Fourier Transform
% F         = Vector of frequency bins
% NORM      = Vector of normalization applied at each time instant
%
% Example : Spectrogram of a tone embedded in noise
%    x = hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128));
%    [stft,T,F] = Ctfrstft(x,1:128,128);
%    imagesc(T,F,abs(stft).^2); axis xy
%
% SEE ALSO : Ctfrsp, Ctfrrsp, Cwindow