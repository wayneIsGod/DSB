% MATLAB script for Illustrative Problem 3.1.
% Demonstration script for DSB-AM. The message signal is 
% +1 for 0 < t < t0/3, -2 for t0/3 < t < 2t0/3, and zero otherwise.
echo on
clear;clc;close all;
load('�W�д���.mat');
load('�ɰ����.mat');
t0=10;                                 % signal duration
fc=10000;                                 % carrier frequency
ts=1/(4*fc);                               % sampling interval
snr=20;                                 % SNR in dB (logarithmic)
fs=1/ts;                                % sampling frequency
df=0.3;                                 % desired freq. resolution
t=[ts:ts:t0];                            % time vector
snr_lin=10^(snr/10);                    % linear SNR
% message signal
%m=[ones(1,t0/(3*ts)),-2*ones(1,t0/(3*ts)),zeros(1,t0/(3*ts)+1)];
c=cos(2*pi*fc.*t);  % carrier signal
y=y.';
u=y.*c;                                 % modulated signal
[M,y,df1]=fftseq(y,ts,df);              % Fourier transform 
M=M/fs;                                 % scaling                   
[U,u,df1]=fftseq(u,ts,df);              % Fourier transform 
U=U/fs;                                 % scaling
[C,c,df1]=fftseq(c,ts,df);              % Fourier transform
f=[0:df1:df1*(length(y)-1)]-fs/2;       % freq. vector
signal_power=spower(u(1:length(t)));    % power in modulated signal
noise_power=signal_power/snr_lin;       % Compute noise power.
noise_std=sqrt(noise_power);            % Compute noise standard deviation.
noise=noise_std*randn(1,length(u));     % Generate noise.
r=u+noise;                              % Add noise to the modulated signal.
[R,r,df1]=fftseq(r,ts,df);              % spectrum of the signal+noise 
R=R/fs;                                 % scaling
%pause  % Press a key to show the modulated signal power.
%{
signal_power
pause  % Press any key to see a plot of the message.
clf
subplot(2,2,1)
plot(t,y(1:length(t)))
xlabel('Time')
title('The message signal')
pause  % Press any key to see a plot of the carrier.
subplot(2,2,2)
plot(t,c(1:length(t)))
xlabel('Time')
title('The carrier')
pause  % Press any key to see a plot of the modulated signal.
subplot(2,2,3)
plot(t,u(1:length(t)))
xlabel('Time')
title('The modulated signal')
pause   % Press any key to see plots of the magnitude of the message and the
    % modulated signal in the frequency domain.
subplot(2,1,1)
plot(f,abs(fftshift(M)))
xlabel('Frequency')
title('Spectrum of the message signal')
subplot(2,1,2)
plot(f,abs(fftshift(U)))
title('Spectrum of the modulated signal')
xlabel('Frequency')
pause  % Press a key to see a noise sample.
subplot(2,1,1)
plot(t,noise(1:length(t)))
title('Noise sample') 
xlabel('Time')
pause  % Press a key to see the modulated signal and noise.
subplot(2,1,2)
plot(t,r(1:length(t)))
title('Signal and noise')
xlabel('Time')
pause  % Press a key to see the modulated signal and noise in freq. domain.
subplot(2,1,1)
plot(f,abs(fftshift(U)))
title('Signal spectrum')
xlabel('Frequency')
subplot(2,1,2)
plot(f,abs(fftshift(R))) 
title('Signal and noise spectrum')
xlabel('Frequency')
%}

