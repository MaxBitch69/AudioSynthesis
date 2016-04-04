function [y, y_test] = test4(x, Fs, nbViolins)
% Time stretching using STFT of signal x, using stretching percentage
% "stretch".
%
% Pitchs = list of new pitchs for pitch shifting, decimals

close all

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
sigma = 45; % standard deviation in ms

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
        pitchs(k) = normrnd(1, 0.1); % 1% of pitch modification
end

%% Metropolis - Hastings - Compute Time Difference
timeDifference = zeros(nbViolins, Nt); % Time Difference
for k = 1:nbViolins
    timeDifference(k,:) = MetropolisHastings(mu, sigma, Nt);

    % Low-frequency sampling, ie smoothing
    filt = 1/10*hanning(200); % Hanning filter
    % filt = 1/40*ones(100,1); % simple smoother, corresponding
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
    filt = 1/10*hanning(30); % simple smoother, corresponding to 1s
    % filt = 1/20*ones(20,1); % simple smoother / Mean filter, corresponding to 1s
    amplitudeModulation(k,:) = filter(filt, 1, amplitudeModulation(k,:)); % Smooth on 1s
end

%% STFT
% Initialisation
puls = 2*pi*I*(0:Nfft-1)'/Nfft; % Canals' pulsations
Xtilde_m = zeros(Nfft, Nt); % Matrix containing fft
Xtilde_m(:,1) = fft(x(1:Nw), Nfft); % 1st fft

% Parameters for time stretching
phase = angle(Xtilde_m(:,1))*ones(1, nbViolins); % initialisation matrix
former_phase = phase;

for k=2:Nt-20  % Loop on timeframes
    % Display progression
    str = sprintf('Treatment progression: %.1f %%', 100*k/Nt);
    disp(str);
    
    %%% ANALYSIS
    deb = (k-1)*I +1; % Beginning - x(n+kI)
    fin = deb + Nw -1; % End
%    tx = x(deb:fin).*w; % Timeframe
        
    % Treatment on signals
    for h = 1:nbViolins

        % Time Difference
        deb1 = deb + floor(timeDifference(h, k)*10^-3*Fs);
        fin1 = deb1 + Nw -1; % fin de trame
        tx= x(deb1:fin1).*w;
        
        % Pitch shifting
        [p, q] = rat(pitchs(h)); % Get fraction
        ttx = resample(tx, q, p); % New time base vector - p/q, p/q times smaller

        % FFT
        X = fft(ttx,Nfft); 

        % Time stretching
        stretch = p/q;%(Nw+floor(timeDifference(h, k)*10^-3*Fs))/length(ttx);%p/q; %Nw/length(ttx); % p/q
        diff_phase = (angle(X) - former_phase(:,h)) - puls;
        diff_phase = mod(diff_phase + pi,-2*pi) + pi;
        diff_phase = (diff_phase + puls) * stretch;

        phase(:,h) = phase(:,h) + diff_phase;
        Y = abs(X).*exp(1i*phase(:,h));
        former_phase(:,h) = angle(X);
        
        %%% SYNTHESIS
%         % Time Difference
%         deb1 = deb + floor(timeDifference(h, k)*10^-3*Fs);
%         fin1 = deb1 + Nw -1; % fin de trame

        % Reconstruction
        ys = real(ifft(Y, 'symmetric')); % TFD inverse
        ys = ys.*ws; % pondération par la fenêtre de synthèse
        
        % Amplitude modulation
        factorAmp = amplitudeModulation(h, k)/sum(amplitudeModulation(:,k)); % Normalisation
        ys = factorAmp.*ys;
        
        y_test(deb:fin, h) = y_test(deb:fin, h) + ys; % Each signal - y_test: stereo 
                                                      % if 2 pitchs: soundsc(y_test, Fs)
        y(deb:fin)=y(deb:fin)+ys; % overlap add - sum of signals
    end
end
end