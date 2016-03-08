function locs = onsetDetection(x, Fs)
% Finds the onset of signal x
%   Finds the onsets of signal x using spectral difference (SD) method

filterLength = 500*10^-3; % 25ms average smoothing
filterPoints = floor(filterLength/75);

filt = 1/filterPoints*ones(1, filterPoints); % Average smoothing filter
x = filter(filt, 1, x); % Filtering

%% Parameters
N = length(x); % Signal length
Nfft = 1024; % Number of points for FFT

L_s = 250*10^-3; % Anaysis windows of 85 ms, here in seconds
Step_s = 75*10^-3; % step of 0.5ms, here in seconds 

L_n = floor(L_s*Fs); % In points
Step_n = floor(Step_s*Fs); % In points

w = hanning(L_n); % Hanning analysis window

%% STFT
Nt = floor((N-L_n)/Step_n); % Number of frames
f = zeros(Nt, 1); % Function of onsets
formerX = zeros(Nfft, 1); % FFT at previous loop

for k = 1:Nt
    deb = (k-1)*Step_n +1;
    fin = deb + L_n -1;
    tx = x(deb:fin).*w; % extract windowed signal
    X = fft(tx, Nfft); % Get FFT

    % Spectral difference method
    diffX = abs(X) - abs(formerX);
    h = ((diffX + abs(diffX))/2).^2;
    f(k) = sum(h);
    
    formerX = X;
end

%% Additionnal smoothing
filterLength = 500*10^-3; % 25ms average smoothing
filterPoints = floor(filterLength/Step_s);

filt = 1/filterPoints*ones(1, filterPoints); % Average smoothing filter
f = filter(filt, 1, f); % Filtering

figure();
plot(f)

%% Find local maxima
[pks, locs] = findpeaks(f)%,'MinPeakDistance',250/75)

end

