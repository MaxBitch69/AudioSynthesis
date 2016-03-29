function [pks, locs, signals, bins, peaks_bin] = onsetDetection(x, Fs)
% Finds the onset of signal x
%   Finds the onsets of signal x using spectral difference (SD) method

%% Parameters
N = length(x); % Signal length
Nfft = 1024; % Number of points for FFT

L_s = 0.085; %250*10^-3; % Anaysis windows of 85 ms, here in seconds
Step_s = 0.005; %75*10^-3; % step of 0.5ms, here in seconds 

L_n = floor(L_s*Fs); % In points
Step_n = floor(Step_s*Fs); % In points

w = hanning(L_n); % Hanning analysis window

%% STFT
Nt = floor((N-L_n)/Step_n); % Number of frames
f = zeros(Nt, 1); % Function of onsets
g = zeros(Nt, 1); % Function of onsets - 2nd method

former_bin = 0; % Bin corresponding to max of frequency
peaks_bin = [];
bins = [];

formerX = zeros(Nfft, 1); % FFT at previous loop

for k = 1:Nt
    deb = (k-1)*Step_n +1;
    fin = deb + L_n -1;
    tx = x(deb:fin).*w; % extract windowed signal
    X = fft(tx, Nfft); % Get FFT

    % Direct method - increase high frequencies's relevancy
    g(k) = 1/Nfft*sum((1:Nfft)'.^2.*abs(X).^2);

    % Spectral difference method
    diffX = abs(X) - abs(formerX);
    h = ((diffX + abs(diffX))/2).^2;
    f(k) = sum(h);

    % Method based on pitch detection
%     bin = binDetection(X);
%     if abs(bin-former_bin) > 4 % Frequency change
%         bins = [bins bin]; % Get bin
%         peaks_bin = [peaks_bin deb]; % Get window's beginning
%         former_bin = bin;
%     end
    [maxValue, indexMax] = max(abs(fft(tx-mean(tx))));
    if abs(indexMax-former_bin) > 50
        bins = [bins indexMax*Fs/L_n];
        peaks_bin = [peaks_bin deb]; % Get window's beginning
        former_bin = indexMax;
    end
    
    formerX = X;
end

figure();
plot(f);
title('Spectral Difference method - f');

figure();
plot(g);
title('Direct method - g');

mean(f)
median(f)
[pks, locs] = findpeaks(f, 'MINPEAKHEIGHT', 3*mean(f), 'MINPEAKDISTANCE', 25);%, Fs); %, 'MinPeakDistance', 3); 

%% Additionnal smoothing
filterLength = 0.05; %25; %L_s; %0.2; %10; %500*10^-3; % 25ms average smoothing
filterPoints = floor(filterLength/Step_s);

filt = 1/filterPoints*ones(1, filterPoints); % Average smoothing filter
f = filter(filt, 1, f); % Filtering
g = filter(filt, 1, g); % Filtering

figure();
plot(f);
title('Spectral Difference method after filtering');

figure();
plot(g);
title('Direct method after filtering');

figure();
plot(peaks_bin, bins);
title('Peaks and pitch shifting');

%% Find local maxima
mean(f)
%[pks, locs] = findpeaks(f, 'MINPEAKHEIGHT', 1.5*mean(f), 'MINPEAKDISTANCE', 10);

%[pks, locs] = findpeaks(f, Fs, 'MinPeakDistance', 3); %, 'MinPeakProminence', 1.3*mean(f));

f = find(f>1.1*mean(f)); % 10% superior to mean

X = {};%zeros(length(f)-1,1); % Vector of onsets
h = 1;

for k = 1:length(f)-1
    if f(k+1) - f(k) > 3
        X{h} = x(f(k):f(k+1));
        h = h+1;
    end
end

%% Get signals
signals = {};%cell(length(pks), 1);
for k=1:length(pks)-1
    debut = (locs(k)-1)*Step_n+1;
    fin = (locs(k+1)-1)*Step_n+1;
    signals{k} = x(debut: fin);
end

signals{length(pks)} = x((locs(length(pks))-1)*Step_n+1:end);

end
