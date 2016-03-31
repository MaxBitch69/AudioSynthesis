function [y, y_test, pitchs] = test4(x, Fs, nbViolins)
% Time stretching using STFT of signal x, using stretching percentage
% "stretch".
%
% Pitchs = list of new pitchs for pitch shifting, decimals

%% Parameters
N = length(x); % Signal's duration
Nw = 2048; % Excerpts of signal
Nfft = 2048; % Precision of FFT

overlap = 0.25; % overlap in %, here 75%
I = floor(Nw*overlap); % Hop size in points

Nt = floor((N-Nw)/I); % Trames/FFT number

y = zeros(N, 1); % Synthesised signal
y_test = zeros(N, nbViolins); % Different channels

% Metropolis-Hastings sampling, for Time Difference
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

%% Compute Pitchs
pitchs = zeros(nbViolins, 1);
for k = 1:nbViolins
        % Pick a random number following normal distribution
        pitchs(k) = normrnd(1, 0.01); % 1% of pitch modification
end

%% Metropolis - Hastings - Compute Time Difference
timeDifference = zeros(nbViolins, Nt); % Time Difference
for k = 1:nbViolins
    timeDifference(k,:) = MetropolisHastings(mu, sigma, Nt);

    % Low-frequency sampling, ie smoothing
    filt = 1/20*ones(20,1); % simple smoother, corresponding
    timeDifference(k,:) = filter(filt, 1, timeDifference(k,:)); % Smooth on 1s
    
    % Display results
    figure();
    plot(timeDifference(k,:));
    str = sprintf('Time Difference for violin n°%d', k);
    title(str);
end

%% Metro - Hastings for amplitude modulation
amplitudeModulation = zeros(nbViolins, Nt); % Time Difference
for k = 1:nbViolins
    amplitudeModulation(k,:) = abs(MetropolisHastings(1, 0.4, Nt)); % Amplitude

    % Low-frequency sampling, ie smoothing
    filt = 1/20*ones(20,1); % simple smoother, corresponding to 1s
    amplitudeModulation(k,:) = filter(filt, 1, amplitudeModulation(k,:)); % Smooth on 1s
end

%% STFT
% Initialisation
puls = 2*pi*I*(0:Nfft-1)'/Nfft; % Canals' pulsations
Xtilde_m = zeros(Nfft, Nt); % Matrix containing fft
Xtilde_m(:,1) = fft(x(1:Nw), Nfft); % 1st fft

% Parameters for time stretching
phase = angle(Xtilde_m(:,1));
former_phase = phase;

for k=2:Nt-1;  % Loop on timeframes
    % Analysis
    deb = (k-1)*I +1; % Beginning - x(n+kI)
    fin = deb + Nw -1; % End
    tx = x(deb:fin).*w; % Timeframe
        
    % Treatment on signals
    for h = 1:nbViolins
        
        % Pitch shifting
        [p, q] = rat(pitchs(h)); % Get fraction
        ttx = resample(tx, q, p); % New time base vector - p/q, p/q times smaller

        % FFT
        X = fft(tx,Nfft); 

        % Time Difference
        deb = deb + floor(timeDifference(h, k)*10^-3*Fs);
        fin = fin + floor(timeDifference(h, k+1)*10^-3*Fs);    

        % Time stretching
        stretch = (fin-deb+1)/length(ttx);
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
    
            
        % Amplitude modulation
        factorAmp = amplitudeModulation(h, k)/sum(amplitudeModulation(:,k)); % Normalisation
        ys = factorAmp.*ys;
        
    
    y(deb:fin)=y(deb:fin)+ys; % overlap add
end

        
        % Synthèse
        % Reconstruction
        ws = hanning(fin-deb+1); % Synthesis window
        ys = ys.*ws; % pondération par la fenêtre de synthèse
        y_test(deb:fin, h) = y_test(deb:fin, h) + ys; % Each signal - y_test: stereo 
                                                      % if 2 pitchs: soundsc(y_test, Fs)
        y(deb:fin)=y(deb:fin)+ys; % overlap add - sum of signals
    end
end

%y = x + mean(abs(x))/mean(abs(y))*y; % Add original signal 

end

