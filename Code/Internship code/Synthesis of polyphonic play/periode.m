function P = periode(x, Fs, Pmin, Pmax, seuil)
% Period estimation on frames of 25ms  
% if no detection, periode = 10ms.Fs -> frequency = 100Hz 
% Si voise = 0, P est égal à 10ms.Fs

%% Read signal
x = x(:);
x = x - mean(x); % Get rid of mean
N = length(x);

%% Parameters
if nargin<5, seuil = 0.7; end; % Detect if there is a period
if nargin<4, Pmax = 1/20; end; % Frequency min 20Hz
if nargin<3, Pmin = 1/1500; end; % max frequency, 1500Hz

Nmin = 1 + ceil(Pmin*Fs);
Nmax = 1 + floor(Pmax*Fs);
Nmax = min([Nmax,N]);

%% Compute FFT
Nfft = 2^nextpow2(2*N-1);
X = fft(x,Nfft);
S = X .* conj(X) / N;
r = real(ifft(S));

%% Find period
[rmax,I] = max(r(Nmin:Nmax));
P = (I+Nmin-2);

%% Check if the period is acceptable
corr = (rmax/r(1)) * (N/(N-P));
if corr < seuil, P = round(10e-3*Fs); end;

end
