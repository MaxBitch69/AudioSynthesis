%% Script for synthesised sound
clear all;
close all;
clc;

%% Parameters
N = 63;
f0 = 1/4;
f1 = f0 + 1/N;
a0 = 1;
a1 = 10;
delta0 = 0;
delta1 = -0.05;

%% Synthesis - Question 3
% Vectors
delta_s = [delta0 delta1]';
f_s = [f0 f1]';
a_s = [a0 a1]';
phi_s = 2*pi*rand(2, 1)-pi;

% Synthesis
s = Synthesis(N, delta_s, f_s, a_s, phi_s);

% Display synthesis
figure();
plot(abs(s));

%% Periodogram - Question 3.1
% Nfft = N
Nfft = N;
x = abs(fft(s, Nfft)).^2/N;
figure();
subplot(1, 2, 1);
f = (0:Nfft-1)/Nfft;

% Display results 
plot(f, db(x)/2); % db(x) = 20*log(x)
str = sprintf('Nfft = N = %d', N); % title
title(str);

% Nfft = 1024 > N
Nfft = 1024;
x = abs(fft(s, Nfft)).^2/N;
subplot(1, 2, 2);
f = (0:Nfft-1)/Nfft;

% Display results 
plot(f, db(x)/2); %db(x) = 20*log(x)       
title('Nfft = 1024 > N');

% NB: en augmentant Nfft, on fait du zero-padding et on augmente donc la précision en fréquence
% Ainsi, pour N=Nfft, on ne voit qu'un pic alors qu'en augmentant la précision, on arrive à déceler un second pic
% juste à côté.
% Les 2 pics sont très proches en fréquence. Il n'est pas impossible que les 2 pics soient mal placés en plus. 
% En effet, les  pics étant très proches, il y a des phénomènes d'interférences qui font que les 2 réduisent leur importance
% On trouve alors d'autres pics, proches certes, mais pas tout à fait aux bons endroits.
% Quand on checke les emplacements en effet, on ne retrouve pas les valeurs de f0 = 0.25 et f = f0+1/N.

% On aura bon augmenter Nfft pour améliorer la précision, le problème est
% ailleurs et ne sera pas résolu. Il faut donc augmenter la résolution
% spectrale

%% HR Method on sound synthesised s - Question 3.2.5
n = 32;
K = 2; % Dimension de l'espace signal

% Method ESPRIT + LeastSquares
[delta_e, f_e] = ESPRIT(s, n, K);
[a_e, phi_e] = LeastSquares(s, delta_e, f_e);

% Comparision values esprit _e and values synthesised _s
fprintf(1,'Synthesised\t\tEstimated:\n');
fprintf(1,'Frequencies:\n');
fprintf(1,'%f\t%f\n\n',[f_s f_e]');
fprintf(1,'Amplitudes:\n');
fprintf(1,'%f\t%f\n\n',[a_s a_e]');
fprintf(1,'Decroissances:\n');
fprintf(1,'%f\t%f\n\n',[delta_s delta_e]');
fprintf(1,'Phases:\n');
fprintf(1,'%f\t%f\n\n',[phi_s phi_e]');

%% Method MUSIC - Question 3.2.6
MUSIC(s, n, K);

