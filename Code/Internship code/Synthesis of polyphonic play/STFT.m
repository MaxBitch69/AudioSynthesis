function Xtilde_m = STFT(x, Fs, Nw, overlap, display)
% Computes signal x's Short Term Fourier Transform
% And returns Xtilde_m back.
%
% Inputs:
%       - x: original signal
%       - Nw = analysis window's length. By default, Nfft = Nw, 
%       and the synthesised and analysis windows have the same length
%       - overlap in %: min of 0.75 = 75% advised
%       - display: 1 to see spectrogram, 0 by default
%
% Output:
%       - original STFT


%% Arguments
if nargin < 5
    display = 0;
end

%% Parameters
N = length(x); % Signal's duration
Nfft = Nw; % FFT points - Nw by default

overlap = 1-overlap; % Conversion 
I = floor(Nw*overlap); % Analysis hop size in points

Nt = floor((N-Nw)/I); % Trames/FFT number

Xtilde_m = zeros(Nfft, Nt); % STFT

%% Windowing
w = hanning(Nw); % Analysis window

%% STFT
for k=2:Nt;  % Loop on timeframes
    % Analysis
    deb = (k-1)*I +1; % frame's beginning - x(n+kI)
    fin = deb + Nw -1; % frame's end
    tx = x(deb:fin).*w; % Timeframe

    Xtilde_m(:,k) = fft(tx, Nfft); % FFT on the timeframe    
end

%% Display - Spectrogramm
if display
    figure();
    spectrogram(x, w, I, Nfft,  'yaxis');
    ylim([0 3000/Fs]); % Cut between 0 and 3000Hz
end

end

