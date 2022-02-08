function pwrplt(X);
% Convert a Gaussian pulse from the time domain to the frequency domain.
% 
% Define signal parameters and a Gaussian pulse, |X|.
 
% Copyright 2015 The MathWorks, Inc.
 
Fs = 15.5;           % Sampling frequency
t = -1:1/Fs:6;  % Time vector 
L = length(t);      % Signal length
 
% X = traces2(:,3);
%% 
% Plot the pulse in the time domain.
figure('Position',[100 100 900 400])
subplot(1,2,1)
plot(t,X)
title('Signal in Time Domain')
xlabel('Time (t)')
ylabel('X(t)')
%% 
% To use the |fft| function to convert the signal to the frequency domain, 
% first identify a new input length that is the next power of 2 from the original 
% signal length. This will pad the signal |X| with trailing zeros in order to 
% improve the performance of |fft|.
n = 2^nextpow2(L);
%% 
% Convert the Gaussian pulse to the frequency domain.
Y = fft(X,n);
%% 
% Define the frequency domain and plot the unique frequencies.
f = Fs*(0:(n/2))/n;
P = abs(Y/n);

subplot(1,2,2)
plot(f,P(1:n/2+1)) 
[v,p]=findpeaks(P);
title('Gaussian Pulse in Frequency Domain')
xlabel('Frequency (f)')
ylabel('|P(f)|')
