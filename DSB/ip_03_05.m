% MATLAB script for Illustrative Problem 3.5.
% Demonstration script for DSB-AM demodulation. The message signal
% is +1 for 0 < t < t0/3, -2 for t0/3 < t < 2t0/3, and zero otherwise.
echo on
clear;clc;close all;
load('ÀWÃÐ´ú¸ê.mat');
load('®É°ì´ú¸ê.mat');
fc=10000;                                 % carrier frequency
t0=10;                                 % signal duration
ts=1/(4*fc);                              % sampling interva
fs=1/ts;                                % sampling frequency
t=[ts:ts:t0];                            % time vector
df=0.3;                                 % desired frequency resolution
snr=0;
snr_lin=10^(snr/10); 
% message signal
m=[ones(1,t0/(3*ts)),-2*ones(1,t0/(3*ts)),zeros(1,t0/(3*ts)+1)];
c=cos(2*pi*fc.*t);                      % carrier signal
m=y.';
u=m.*c;                                 % modulated signal                          

signal_power=spower(u(1:length(t)));    % power in modulated signal
noise_power=signal_power/snr_lin;       % Compute noise power.
noise_std=sqrt(noise_power);            % Compute noise standard deviation.
noise=noise_std*randn(1,length(u));     % Generate noise.
r=u+noise; 
y=r.*c;                                 % mixing
[M,m,df1]=fftseq(m,ts,df);              % Fourier transform 
M=M/fs;                                 % scaling
[U,u,df1]=fftseq(u,ts,df);
[R,r,df1]=fftseq(r,ts,df);              % Fourier transform 
R=R/fs;                                 % scaling
[Y,y,df1]=fftseq(y,ts,df);              % Fourier transform
Y=Y/fs;                                 % scaling                              
f_cutoff=10000;                           % cutoff freq. of the filter
n_cutoff=floor(f_cutoff/df1);                % Design the filter.
f=[0:df1:df1*(length(y)-1)]-fs/2;
H=zeros(size(f));                    
H(1:n_cutoff)=2*ones(1,n_cutoff);    
H(length(f)-n_cutoff+1:length(f))=2*ones(1,n_cutoff);
DEM=H.*Y;                   % spectrum of the filter output
dem=real(ifft(DEM))*fs;             % filter output
pause % Press a key to see the effect of mixing.
sound(y,fs);
length(y)
clf
subplot(3,1,1)
plot(f,fftshift(abs(M)))
title('Spectrum of the Message Signal')
xlabel('Frequency')
subplot(3,1,2)
plot(f,fftshift(abs(R)))
title('Spectrum of the Modulated Signal')
xlabel('Frequency')
subplot(3,1,3)
plot(f,fftshift(abs(Y)))
title('Spectrum of the Mixer Output')
xlabel('Frequency')
pause % Press a key to see the effect of filtering on the mixer output.
clf
subplot(3,1,1)
plot(f,fftshift(abs(Y)))
title('Spectrum of the Mixer Output')
xlabel('Frequency')
subplot(3,1,2)
plot(f,fftshift(abs(H)))
title('Lowpass Filter Characteristics')
xlabel('Frequency')
subplot(3,1,3)
plot(f,fftshift(abs(DEM)))
title('Spectrum of the Demodulator output')
xlabel('Frequency')
pause % Press a key to compare the spectra of the message and the received signal.
clf
subplot(2,1,1)
plot(f,fftshift(abs(M)))
title('Spectrum of the Message Signal')
xlabel('Frequency')
subplot(2,1,2)
plot(f,fftshift(abs(DEM)))
title('Spectrum of the Demodulator Output')
xlabel('Frequency')
pause % Press a key to see the message and the demodulator output signals.
subplot(2,1,1)
plot(t,m(1:length(t)))
title('The Message Signal')
xlabel('Time')
subplot(2,1,2)
plot(t,dem(1:length(t)))
title('The Demodulator Output')
xlabel('Time')
