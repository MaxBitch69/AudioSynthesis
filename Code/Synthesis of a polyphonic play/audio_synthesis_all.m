%% Script for synthesis of a violin sound
clear all;
close all;
clc;

%% Read signal
[x, Fs] = audioread('Violin/Vibrato_G3.aif');
x = x(:,1); %(x(:,1)+x(:,2))/2; % Stereo --> Mono
x = x(1:100000);

soundsc(x, Fs);
%N = length(x);

% Display original sound
figure();
plot(real(x));
title('Real sound');

%% HR Method on real sound x 
% Parameters
N = 1000; % Length of signal excerpted - Analysis window length
K = 100; % Signal space dimension
n = 300; % 256 % Full space dimension (signal + noise)
         % Noise space dimension = n-K
l = N-n+1; % N = n+l-1, l completes n to get N


%% Get signal properties
for k = 1:floor(length(x)/N)-1
    
    h = (k-1)*N+1
    
    % Signal extraction
    y = x(h:h+N);

    % Method ESPRIT + LeastSquares
    [delta_e, f_e] = ESPRIT(y, n, K);
    [a_e, phi_e] = LeastSquares(y, delta_e, f_e);

    f(:,k) = f_e;
    a(:,k) = a_e;
    phi(:,k) = phi_e;
    delta(:,k) = delta_e;
    
end

%% Re-synthesis of the signal
y = []';
for k = 1:size(phi, 2)
   
    k 
    
    length = 1.5*N; % durée plus longue pour mettre en évidence les résonances
    s = Synthesis(length, delta(:,k), 1.5*f(:,k), a(:, k), phi(:, k));

    y = [y; s];
    size(y)
end


%% Display results
figure();
plot(real(y));
title('Sound synthesised');

% Listen to result
soundsc(real(y), Fs);


%% Filtering signal
coeffs = 1/10*ones(1, 10);
y = filter(coeffs, 1, y);

% Listen to result
soundsc(real(y), Fs);

