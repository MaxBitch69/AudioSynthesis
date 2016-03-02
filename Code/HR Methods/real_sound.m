%% Script for real audio sounds
clear all;
close all;
clc;

%% Read signal
[x, Fs] = audioread('ClocheA.WAV');
%soundsc(x, Fs);
N = length(x);

% Display original sound
figure();
plot(real(x));
title('Real sound');

%% Periodogram - Question 4.1
% Nfft = 512 
Nfft = 512;
perio = abs(fft(x, Nfft)).^2/N;
f = (0:Nfft-1)/Nfft;

% Display results 
figure();
plot(f, db(perio)/2); %db(x) = 20*log(x) 
title('Nfft = 512');

% NB: really hard to see the right frequencies, that's why we're going to
% use the HR methods

%% HR Method on real sound x - Question 4.2
n = 512;
K = 54;
N = 1535;
l = 2*n;

% Extraction of the signal
x = x(10000:10000+N);

% Method ESPRIT + LeastSquares
[delta_e, f_e] = ESPRIT(x, n, K);
[a_e, phi_e] = LeastSquares(x, delta_e, f_e);

%% Re-synthesis of the signal
length = 15*N; % durée plus longue pour mettre en évidence les résonances
s = Synthesis(length, delta_e, f_e, a_e, phi_e);

% Display results
figure();
plot(real(s));
title('Sound synthesised');

% Listen to result
soundsc(real(s), Fs);


