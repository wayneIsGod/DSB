%{
clear;clc;close all;
Fs=40000;
N=Fs*10;
FFTsize=1024;
y=wavrecord(N,Fs);
Y=spectrum(y,FFTsize);
 
Freq=[0:Fs/FFTsize:Fs/2];
Time=[1:N]/Fs;
subplot(2,1,1),plot(Time,y),
ylabel('Amplitude'),xlabel('Time(s)');
subplot(2,1,2),plot(Freq,10*log10(Y/max(Y))),
ylabel('Spectrum(dB)'),xlabel('Frequency(Hz)');
sound(y,Fs)
%}

y=load('办代戈.mat');
audiowrite('办代戈.wav',y,Fs);
audiowrite('繵眯代戈.wav',Y,Fs);
Y=load('繵眯代戈.mat');
sound(y,Fs); 
%{
y=load("办代戈.mat");
Y=load("繵眯代戈.mat");
%}