%CAMBIFUNB Ambiguity function
% [NAF,TAU,XI] = Cambifunb(X,TAU,N);
%
% Computes the narrow-band ambiguity function
%
% X           = Signal for which the AF should be computed
% TAU         = Vector of lag values     (default : -Nx/2:Nx/2).
% N           = Number of frequency bins (default : length(X)).
% NAF         = doppler-lag representation, with the doppler bins stored 
% 	        in the rows and the time-lags stored in the columns.
% XI          = vector of doppler values.
%
% Note: The cross ambiguity function is not implemented yet
%
% Example :
%
%  [AF,TAU,XI]=Cambifunb(hilbert(sin(2*pi*0.25*(1:128)))+0.5*randn(1,128));
%  imagesc(XI,TAU,abs(AF)); axis xy; 
%  xlabel('Time lag'); ylabel('Frequency lag')