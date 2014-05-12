%CWINDOW Window generation.
% H=WINDOW(N,NAME,PARAM,PARAM2)
%
% Creates a window of length N with a given shape.
%
% N      = length of the window
% NAME   = name of the window shape
% PARAM  = optional parameter
% PARAM2 = second optional parameters
% 
% Possible names are :
% 'Hamming', 'Hanning', 'Nuttall',  'Papoulis', 'Harris',
% 'Rect',    'Triang',  'Bartlett', 'BartHann', 'Blackman'
% 'Gauss',   'Parzen', 'Dolph',    'Hanna', 'Nutbess', 'spline'
% 
% 	For the gaussian window, an optionnal parameter K
% 	sets the value at both extremities. The default value is 0.005
% 
% 	For the Spline windows, h=window(N,'spline',nfreq,p)
% 	yields a spline weighting function of order p and frequency
% 	bandwidth proportional to nfreq.
% 
% 	Example : h=window(256,'Gauss',0.005); plot(h);