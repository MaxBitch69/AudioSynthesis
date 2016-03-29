function y = timestretch(x, stretch)
% Time stretching using STFT of signal x, using stretching percentage
% "stretch".

%% Parameters
N = length(x); % Signal's duration
Nw = 1024; %round(L_s*Fs); % Analysis window's length in points
Nfft = 1024; % FFT points

overlap = 0.125; % overlap in %, here 82%
I = floor(Nw*overlap); % Hop size in points
R = floor(stretch*I); % Hop size for stretching

Nt = floor((N-Nw)/I); % Trames/FFT number
y = zeros(floor(N*stretch)+Nw, 1); % Synthesised signals

%% Windowing
w = hanning(Nw); % Analysis window
ws = w; % Synthesis window

% Check for perfect reconstruction - ie h = w.*ws == 1
h = w.*ws;
output = ola(h, R, 30); % Check reconstruction

% window's normalisation - w.*ws == 1
amp = max(output);
ws = ws./amp; 

%% STFT
% Initialisation
puls = 2*pi*I*(0:Nfft-1)'/Nfft; % Canals' pulsations
Xtilde_m = zeros(Nfft, Nt); % Matrix containing fft
Xtilde_m(:,1) = fft(x(1:Nw), Nfft); % 1st fft

% Parameters for time stretching
phase = angle(Xtilde_m(:,1));
former_phase = phase;

% STFT
for k=2:Nt;  % Loop on timeframes
    % Analysis
    deb = (k-1)*I +1; % Beginning - x(n+kI)
    fin = deb + Nw -1; % End
    tx = x(deb:fin).*w; % Timeframe
    X = fft(tx,Nfft); % FFT
    
    % Time stretching
    diff_phase = (angle(X) - former_phase) - puls;
    diff_phase = mod(diff_phase + pi,-2*pi) + pi;
    diff_phase = (diff_phase + puls) * stretch;

    phase = phase + diff_phase;

    Y = abs(X).*exp(1i*phase);

    former_phase = angle(X);
        
    % Synthèse
    deb = (k-1)*R +1; % début de trame - on écarte les instants de synthèse
    fin = deb + Nw -1; % fin de trame
    
    % Reconstruction
    ys = real(ifft(Y, 'symmetric')); % TFD inverse
    ys = ys.*ws; % pondération par la fenêtre de synthèse
    y(deb:fin)=y(deb:fin)+ys; % overlap add
end

y = y(1:floor(N*stretch)); % Resizing of vector

end

