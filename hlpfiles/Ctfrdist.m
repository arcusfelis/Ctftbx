%CTFRDIST Time frequency distance 
% distance = Ctfrdist(TFR1,TFR2,dist_name,dist_coef);
%
% Computes distances between Time Frequency Representations (TFRs)
%
% TFR1        = first TFR
% TFR2        = second TFR, the same size as previous one
% dist_name   = identifier of the distance to compute
% dist_coef   = optional parameter necessary for some distance
%               measures
%
% dist_name is one of the following :
%
% -------- distances without a normalization of the TFRs ----------
% 'Lq'         : Lq distance, 'dist_coef > 0' is the parameter q 
% 'Quadratic'  : Quadratic distance, no 'dist_coef' required
% 'Correlation': Corelation distance,  no 'dist_coef' required
% --------- distances with a normalization of the TFRs ------------
% 'Kolmogorov' : Kolmogorov distance, no 'dist_coef' required
% 'Kullback'   : Kullback distance, no 'dist_coef' required
% 'Chernoff'   : Chernoff distance, one 'dist_coef in [0;1]' required
% 'Matusita'   : Matusita distance, one 'dist_coef >= 1' required
% 'NLq'        : Normalized Lq distance,'dist_coef > 0' is q
% 'LSD'        : Log Spectral deviation, one 'dist_coef > 0' required
% 'Jensen'     : Jensen distance, one 'dist_coef > 0' required
%
%
% Example : AF1 = Cambifunb(hilbert(sin(2*pi*0.25*(1:128)))+0.5*randn(1,128));
%           AF2 = Cambifunb(hilbert(sin(2*pi*0.15*(1:128)))+0.5*randn(1,128));
%           ker = Ctfrker(128,128,'Wigner');
%           Wigner_TFR1 = Caf2tfr(AF1,ker);
%           Wigner_TFR2 = Caf2tfr(AF2,ker);
%           dist_kolm = Ctfrdist(Wigner_TFR1,Wigner_TFR2,'Kolmogorov');
%
% SEE ALSO : Cambifunb, Ctfrker, Caf2tfr