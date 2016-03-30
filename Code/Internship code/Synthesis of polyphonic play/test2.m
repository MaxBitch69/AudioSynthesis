function [y, y_test, timeDifference] = test2(x, Fs, pitchs)
% Time stretching using STFT of signal x, using stretching percentage
% "stretch".
%
% Pitchs = list of new pitchs for pitch shifting, decimals

%% Parameters
N = length(x); % Signal's duration
Nw = 1*Fs;%0.5*Fs; % Excerpts of signals

overlap = 0.25; %0.25; % overlap in %, here 75%
I = floor(Nw*overlap); % Hop size in points

Nt = floor((N-Nw)/I); % Trames/FFT number

numberShiftings = length(pitchs); % Number of pitch-shiftings
y = zeros(N, 1); % Synthesised signal
y_test = zeros(N, numberShiftings);

% Metropolis-Hastings sampling
mu = 0; % Mean for normal distribution
sigma = 40; % standard deviation in ms

%% Windowing
w = hanning(Nw); % Analysis window
ws = w; % Synthesis window

% Check for perfect reconstruction - ie h = w.*ws == 1
h = w.*ws;
output = ola(h, I, 30); % Check reconstruction

% window's normalisation - w.*ws == 1
amp = max(output);
w = w./amp; 

%% Metropolis - Hastings - Compute Time Difference
timeDifference = zeros(numberShiftings, Nt); % Time Difference
for k = 1:numberShiftings
    timeDifference(k,:) = MetropolisHastings(mu, sigma, Nt);
end
    
%% Metro - Hastings for amplitude modulation
amplitudeModulation = zeros(numberShiftings, Nt); % Time Difference
for k = 1:numberShiftings
    amplitudeModulation(k,:) = abs(MetropolisHastings(1, 0.5, Nt)); % Amplitude
end

%% STFT
for k=2:Nt-1;  % Loop on timeframes
    % Analysis
    deb = (k-1)*I +1; % Beginning - x(n+kI)
    fin = deb + Nw -1; % End
    tx = x(deb:fin).*w; % Timeframe
    
    k
    
    % Treatment on signals
    for h = 1:numberShiftings
        
        % Pitch shifting
        [p, q] = rat(pitchs(h)); % Get fraction
        ttx = resample(tx, q, p); % New time base vector - p/q, p/q times smaller
        
        % Time Difference
        deb = deb + floor(timeDifference(h, k)*10^-3*Fs);
        fin = fin + floor(timeDifference(h, k+1)*10^-3*Fs);
        
        % Time Stretching
        ys = timestretch(ttx, (fin-deb+1)/length(ttx));
        ys = [ys; zeros(fin-deb+1-length(ys),1)]; % Correcting mistakes led by floor()
        
        % Amplitude modulation
        factorAmp = amplitudeModulation(h, k)/sum(amplitudeModulation(:,k)); % Normalisation
        ys = factorAmp*ys;
        
        % Synthèse
        % Reconstruction
        ws = hanning(fin-deb+1); % Synthesis window
        ys = ys.*ws; % pondération par la fenêtre de synthèse
        y_test(deb:fin, h) = y_test(deb:fin, h) + ys; % Each signal - y_test: stereo 
                                                      % if 2 pitchs: soundsc(y_test, Fs)
        y(deb:fin)=y(deb:fin)+ys; % overlap add - sum of signals
    end
end

end

