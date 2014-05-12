%CAF2TFR From ambiguity plane to time frequency plane
% TFR = Caf2tfr(AF,kernel)
%
% Computes the Cohen's group Time Frequency Representation (TFR)
% of a signal given its ambiguity function and the kernel
% of the TFR.
%
% AF      Ambiguity function of the signal
% kernel  kernel of the TFR
%
% Warning : AF and kernel are matrices of the same size
%
% Example : Wigner TFR of a tone embedded in noise
%   AF = Cambifunb(hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128)));
%   ker = Ctfrker(128,128,'Wigner');
%   Wigner_TFR = Caf2tfr(AF,ker);
%   imagesc(1:128,linspace(0,0.5,128),Wigner_TFR); axis xy;
%
% SEE ALSO : Cambifunb, Ctfrker