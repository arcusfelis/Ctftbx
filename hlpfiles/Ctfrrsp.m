%CTFRRSP Reassigned Spectrogram
% [SP_reas,SP,FIELD] = Ctfrrsp(X,T,N,H)
%
% Computes the reassigned spectrogram, the spectrogram and the
% field of reassignment vectors
%
% X        = Analyzed signal
% T        = the time instant(s)   (default : 1:length(X)).
% N        = number of frequency bins (default : length(X)).
% H        = frequency smoothing window, (default : Hamming(N/4)).
%
% SP_reas  = Reassigned spectrogram 
% SP       = Spectrogram (not reassigned)
% FIELD    = Field of reassignment vectors
%
% Example: SP_reas=Ctfrrsp(hilbert(sin(2*pi*0.25*(1:128)))+0.5*randn(1,128));
%          imagesc(SP_reas); axis xy;
%
% SEE ALSO : Ctfrreas, Cwindow
