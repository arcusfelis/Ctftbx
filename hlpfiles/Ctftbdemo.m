function Ctftbdemo();
%-----------------------------------
% ANSI C Time Frequency Toolbox
%-----------------------------------
%        Matlab demo file
%-----------------------------------
close all
echo on

disp('---------------------------------------------')
disp(' Welcome to the ANSI C Time frequency toolbox')
disp(' demonstration file')
disp('---------------------------------------------')

disp('--------------------------------------------')
disp(' First of all, we create FM a signal')
disp(' with additive noise')
disp('---------------------------------------------')
t=1:128;
sig=sin(2*pi*(0.3+0.2*t+0.001*t.^2))+0.5*randn(1,128);
sig=hilbert(sig);

disp('--------------------------------------------')
disp(' this signal looks like this')
disp('---------------------------------------------')
figure(1)
plot(1:128,real(sig),'b',1:128,imag(sig),'r');
xlabel('Time (points)')
ylabel('Amplitude')
title('Analyzed signal (blue=real part, red=imag part)')
axis([0 128 -2.5 2.5]);

disp(' Press any key to continue')
pause

disp('--------------------------------------------')
disp(' Let''s have a look to a spectrogram')
disp('---------------------------------------------')
[TFR,T,F]=Ctfrsp(sig,1:128,128,Cwindow(31,'hamming'));
figure(1)
imagesc(T,F,TFR);
axis xy
xlabel('Time (points)');
ylabel('Normalized frequency');
title('Spectrogram');

disp(' Press any key to continue')
pause


disp('--------------------------------------------')
disp(' What about a  Reassigned spectrogram ?')
disp('---------------------------------------------')
TFR=Ctfrrsp(sig,1:128,128,Cwindow(31,'hamming'));
figure(1)
imagesc(T,F,TFR);
axis xy
xlabel('Time (points)');
ylabel('Normalized frequency');
title('Reassigned Spectrogram');


disp(' Press any key to continue')
pause

disp('--------------------------------------------')
disp(' Or a Choi-Williams Representation ?')
disp('--------------------------------------------')
[TFR,T,F]=Ctfrcw(sig);
figure(1)
imagesc(T,F,TFR);
axis xy
xlabel('Time (points)');
ylabel('Normalized frequency');
title('Choi-Williams representation');


disp(' Press any key to continue')
pause

disp('--------------------------------------------')
disp(' We could try a lot more TFR kernels, ')
disp(' such as the following')
disp('--------------------------------------------')
echo off
figure(1);
[TFR,T,F]=Ctfrbj(sig);
subplot(2,3,1);
imagesc(T,F,TFR);
axis xy
xlabel('Time (points)');
ylabel('Normalized frequency');
title('Born-Jordan')

[TFR,T,F]=Ctfrbud(sig);
subplot(2,3,2);
imagesc(T,F,TFR);
axis xy
xlabel('Time (points)');
ylabel('Normalized frequency');
title('Butterworth')

[TFR,T,F]=Ctfrgrd(sig);
subplot(2,3,3);
imagesc(T,F,TFR);
axis xy
xlabel('Time (points)');
ylabel('Normalized frequency');
title('Generalized rectangular')

[TFR,T,F]=Ctfrmh(sig);
subplot(2,3,4);
imagesc(T,F,TFR);
axis xy
xlabel('Time (points)');
ylabel('Normalized frequency');
title('Generalized Rectangular')
[TFR,T,F]=Ctfrmh(sig);

[TFR,T,F]=Ctfrwv(sig);
subplot(2,3,5);
imagesc(T,F,TFR);
axis xy
xlabel('Time (points)');
ylabel('Normalized frequency');
title('Wigner-Ville')

[TFR,T,F]=Ctfrpwv(sig);
subplot(2,3,6);
imagesc(T,F,TFR);
axis xy
xlabel('Time (points)');
ylabel('Normalized frequency');
title('Pseudo Wigner-Ville')
echo on
disp(' Press any key to continue')
pause
subplot(1,1,1)
disp(' Now, why not try do adapt the TFR kernel')
disp(' to the signals. This is done in 4 steps')
disp('--------------------------------------------')
disp(' 1 - Compute the ambiguity function of the signal')
disp('--------------------------------------------')
figure(1)
[AF,TAU,XI]=Cambifunb(sig,-63:64,128);
imagesc(TAU,XI,abs(AF));
axis xy
xlabel('Time lag (points)');
ylabel('Normalized frequency lag');
title('Ambiguity function (module)')

disp(' Press any key to continue')
pause

disp('--------------------------------------------')
disp(' 2 - Choose a kernel shape among the possible')
disp(' choices (see Ctfrker). We choose a Generalized')
disp(' Marginals Choi-Williams kernel with one branch')
disp('--------------------------------------------')
kernel=Ctfrker(128,128,'gmcwk',[1 0.2]);
imagesc(TAU,XI,kernel);
axis xy
xlabel('Time lag (points)');
ylabel('Normalized frequency lag');
title('Initial one branch kernel')

disp(' Press any key to continue')
pause

disp('---------------------------------------------')
disp(' 3 - Optimize the parameters of this kernel')
disp(' to fit the analyzed signal')
disp('---------------------------------------------')
d_min=1e100;
sigma=0.1;
theta_min=0;
for theta=linspace(0,2*pi,30);
  kernel=Ctfrker(128,128,'gmcwk',[sigma theta]);
  d=Ctfrdist(abs(AF),kernel,'kullback');
  if (d<d_min),
    d_min=d;
    theta_min=theta;
  end
end

disp('---------------------------------------------')
disp(' 4 - Display the optimal kernel')
disp('---------------------------------------------')
kernel=Ctfrker(128,128,'gmcwk',[sigma theta_min]);
imagesc(TAU,XI,kernel);
axis xy
xlabel('Time lag (points)');
ylabel('Normalized frequency lag');
title('One branch kernel')

disp(' Press any key to continue')
pause
disp('---------------------------------------------')
disp(' 4 - here is the optimal TFR for this signal')
disp('---------------------------------------------')

TFR=Caf2tfr(AF,kernel);
imagesc(1:128,linspace(0,0.5,128),TFR);
axis xy
xlabel('Time(points)');
ylabel('Normalized frequency');
title('Optimal One branch kernel TFR')

disp('---------------------------------------------')
disp('Kernel optimization is simple to implement !!')
disp('---------------------------------------------')
