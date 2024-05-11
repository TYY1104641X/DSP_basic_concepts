%function [outputArg1,outputArg2] = plot_a_sinusoid(inputArg1,inputArg2)
% Function: plot a sin function:
%   x(t)=Acos(2pi f t+alpha)   \to   x(n)=Acos(2pi f nTs+alpha)
%   x(n)=Acos(2pi f/Fs n+alpha)=Acos(w0 n+alpha), where w0 in [0,2pi]
%   Author: Yuanyuan Tang
%   Data: May 9, 2024


clc;
close all;
clear all;

%% Initialize parameters

    n=1:100;
    A=10;
    f=2*10^2;    % Frequency 
    Fs=10^3;     % Sampling rate
    w0=2*pi*f/Fs;
    alpha=0.1;   % Phase offset

%% Sin function

    x=A*cos(w0*n+alpha);

%% Plot figure
    figure(1);
    plot(n/Fs,x);
    xlabel('time [s]'),ylabel('Amplitude');
    title('x(n)=2*10 cos(2 \pi*f/Fs+\alpha)');




%% FFT
    X=fft(x);
    k=1:1:length(x);   % Number of samples
    F=k*Fs/length(x);   % Frequency range
    figure(2);
    
    subplot(3,1,1);
    plot(F,abs(X));
    xlabel('Frequency[Hz]'),ylabel('Amplitude');
    subtitle('fft(x)');
    title('x(n)=2*10 cos(2 \pi*200/1000+\alpha)');
    subplot(3,1,2);
    X1=fftshift(X);
    k=-length(x)/2:1:length(x)/2-1;   % Number of samples
    F1=k*Fs/length(X1);   % Frequency range
    plot(F1,abs(X1));
    xlabel('Frequency[Hz]'),ylabel('Amplitude');
    subtitle('fftshift(x)');
    y=ifft(fftshift(X1));
    subplot(3,1,3);
    plot(x,'r--*');
    hold on;
    plot(y,'b-o');
    hold off;
    legend('x','y=ifft(fft(x))')
    xlabel('time[Hz]'),ylabel('Amplitude');
    subtitle('x(n)=2*10 cos(2 \pi*200/1000+\alpha)');


    %% Convolution channel 
    h=[1,0,0,0.5,0,-0.2,0,0.1]; % impulse response
    %h=[1,0,0,0,0,0.2,0,0]; % impulse response
    %length(h);
    n=0:200; 
    x=2*cos(0.1*pi*n); % input signal
    y=conv(x,h); % output signal
    figure(3);
    subplot(2,1,1);
    plot(x);
    subtitle('x(n)');
    subplot(2,1,2);
    plot(y);  %(length(h)/2-2:1:length(y))
    subtitle('y(n)=x(n)*h(n)');




    %% Autocorrelation
    % data
     figure(4);
    sigx=sqrt(2);
    x=sigx*randn(1,1000);
    subplot(2,1,1);
    plot(x);
    subtitle('x(n)');
    % autocorrelation
    N=length(x);
    rx=xcorr(x)/N;
    max_lag=length(x)-1;
    m=-max_lag:max_lag;
    subplot(2,1,2);
    plot(m,rx);
    subtitle('Autocorrelation x(n)');




    %% Upsampling and downsampling
    n=1:40;
    A=10;
    f=1*10^2;    % Frequency 
    Fs=10^3;     % Sampling rate
    w0=2*pi*f/Fs;
    alpha=0.1;   % Phase offset

%% Sin function

    x=A*cos(w0*n+alpha);
    % Upsampling 
    y = upsample(x,4);
    


    %% FFT 
    Fx=fftshift(fft(x));
    Fy=fftshift(fft(y));
    k=-length(x)/2:1:length(x)/2-1;   % Number of samples
    F1=k*Fs/length(x);   % Frequency range
    figure(6);
    subplot(5,1,1);
    plot(F1,Fx);
    subtitle('Fx');
    subplot(5,1,2);
    k=-length(y)/2:1:length(y)/2-1;   % Number of samples
    F2=k*Fs*4/length(y);   % Frequency range
    plot(F2,Fy);
    subtitle('Fy');
    subplot(5,1,3);
    lpf=rectangularPulse(-200,200, F2);
    plot(F2,lpf);
    subtitle('low-pass filter');
    subplot(5,1,4);
    Fy=Fy.*lpf;
    plot(F2,Fy);
    subtitle('Fy1');
    y1=ifft(fftshift(Fy));
    %subplot(5,1,5);
    
    % Downsampling
    c = downsample(y1,4)*4;
    k=-length(c)/2:1:length(c)/2-1;   % Number of samples
    Fc=k*Fs/length(c);   % Frequency range
    subplot(5,1,5);
    plot(Fc, fftshift(fft(c)));
    subtitle('Fc');

    
    % plots the data sequence
    figure(5);
    subplot(4,1,1)
    stem(x);
    subtitle('x(n)=2*10 cos(2 \pi*100/1000+\alpha)')
    subplot(4,1,2)
    stem(y);
    subtitle('y=upsampling(x,4)')
    subplot(4,1,3)
    stem(y1);
    subtitle('y1=lpf(y) in [-200,200]');

    subplot(4,1,4)
    stem (c,'*');
    hold on;
    stem(x,'o');
    hold off;
    legend('x','c');
    subtitle('c=4*downsampling(y1,4)');


    
    




    

%end