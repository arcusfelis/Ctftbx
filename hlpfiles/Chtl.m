%CHTL	Hough transform for detection of lines in images.
% [HT,RHO,THETA]=CHTL(IM,M,N).
%
% From an image IM, computes the integration of the values
% of the image over all the lines. The lines are parametrized 
% using polar coordinates. The origin of the coordinates is fixed
% at the center of the image, and theta is the angle between the
% VERTICAL axis and the perpendicular (to the line) passing through 
% the origin.
%
% IM    = image to be analyzed (size Xmax x Ymax).
% M     = desired number of samples along the radial axis 
%         (default : Xmax).
% N     = desired number of samples along the azimutal (angle) axis
%         (default : Ymax). 
%
% HT    = output matrix (MxN matrix).
% RHO   = sequence of samples along the radial axis.
% THETA = sequence of samples along the azimutal axis.
%
% Example :
%   x = hilbert(sin(2*pi*0.25*(1:128))+0.5*randn(1,128));
%   tfr=Ctfrwv(x); [HT,RHO,THETA]=Chtl(tfr);
%   imagesc(THETA,RHO,HT); axis xy; ylabel('rho'); xlabel('theta')
%
% SEE ALSO : Ctfrsp

