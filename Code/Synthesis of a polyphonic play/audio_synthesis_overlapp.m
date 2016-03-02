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
n = 256; % Full space dimension (signal + noise)
         % Noise space dimension = n-K
l = N-n+1; % N = n+l-1, l completes n to get N

% Parameters overlapp
overlap = 4; % Overlap of 50%
I = floor(N/overlap); 
Nt = fix((length(x)-N)/I); % Nb of trames to compute on

%% Get signal properties - Analysis
for k = 1:Nt
    deb = (k-1)*I+1
    fin = deb+N-1;

    % Signal extraction
    y = x(deb:fin);
    
    % Method ESPRIT + LeastSquares
    [delta_e, f_e] = ESPRIT(y, n, K);
    [a_e, phi_e] = LeastSquares(y, delta_e, f_e);

    f(:,k) = f_e;
    a(:,k) = a_e;
    phi(:,k) = phi_e;
    delta(:,k) = delta_e;
end

%% Re-synthesis of the signal
stretch = 1; % Time stretching
pitch = 1.3; % Pitch shifting

Ns = floor(stretch*N); % durée plus longue pour mettre en évidence les résonances
Is = floor(I*stretch); % Overlap proportional to stretch
% Or Is = floor(Ns*stretch)

y = zeros(floor(stretch*length(x)), 1);

for k = 1:Nt   
    k 
    
    % Synthesis
    s = Synthesis(Ns, delta(:,k), pitch*f(:,k), a(:, k), phi(:, k));

    deb = (k-1)*Is+1;
    fin = deb+Ns-1;
    
    y(deb:fin) = y(deb:fin) + s;
end

y = y/overlap;

%% Display results
figure();
plot(real(y));
title('Sound synthesised');

% Listen to result
soundsc(real(y), Fs);


%% Filtering signal
coeffs = 1/N*ones(1, N);
y = filter(coeffs, 1, y);

% Listen to result
soundsc(real(y), Fs);

