%CTFRREAS Time Frequency Representation reassignment
% TFR_reas = Ctfrreas(TFR,field_time,field_freq)
%
% Reassigns the pixels of a TFR, given a field of reassignment
% vectors.
%
% TFR         = TFR to be reassigned  
% field_time  = Time component of the field of reassignment
%               vectors
% field_freq  = Frequency component of the field of
%               reassignment vectors
% 
% TFR_reas    = Reassigned TFR
%
% Example:
% [SP_reas,SP,field]=Ctfrrsp(hilbert(sin(2*pi*0.25*(1:128)))+0.5*randn(1,128));
% SP_reas2=Ctfrreas(SP,real(field),imag(field));
% imagesc(SP_reas); axis xy;
% figure; imagesc(SP_reas2); axis xy;
%
% SEE ALSO : Ctfrrsp
