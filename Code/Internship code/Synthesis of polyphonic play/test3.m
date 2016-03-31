function [y, y_test] = test3(x, Fs, nbViolins)
% Time stretching using STFT of signal x, using stretching percentage
% "stretch".
%
% Pitchs = list of new pitchs for pitch shifting, decimals

%% Parameters
N = length(x); % Signal's duration
Nw = 0.5*Fs; %0.5*Fs; % Excerpts of signal, here 1s

overlap = 0.125; %0.25; % overlap in %, here 75%
I = floor(Nw*overlap); % Hop size in points

Nt = floor((N-Nw)/I); % Trames/FFT number

y = zeros(N, 1); % Synthesised signal
y_test = zeros(N, nbViolins);

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

%% Compute Pitchs
pitchs = zeros(nbViolins, Nt);
for k = 1:nbViolins
    for h = 1:Nt
        % Pick a random number following normal distribution
        pitchs(k, h) = normrnd(1, 0.05); % 10% modification
    end
end

%% Metropolis - Hastings - Compute Time Difference
timeDifference = zeros(nbViolins, Nt); % Time Difference
for k = 1:nbViolins
    timeDifference(k,:) = MetropolisHastings(mu, sigma, Nt);

    % Display results
    figure();
    plot(timeDifference(k,:));
    str = sprintf('Time Difference for violin n°%d', k);
    title(str);
end

%% Metro - Hastings for amplitude modulation
amplitudeModulation = zeros(nbViolins, Nt); % Time Difference
for k = 1:nbViolins
    amplitudeModulation(k,:) = abs(MetropolisHastings(1, 0.5, Nt)); % Amplitude
end

%% STFT
for k=2:Nt-1;  % Loop on timeframes
    % Analysis
    deb = (k-1)*I +1; % Beginning - x(n+kI)
    fin = deb + Nw -1; % End
    tx = x(deb:fin).*w; % Timeframe
    
    % Display
    str = sprintf('Treatment progression: %.1f %%', 100*k/Nt);
    disp(str);
    
    % Treatment on signals
    for h = 1:nbViolins
        
        % Pitch shifting
        [p, q] = rat(pitchs(h, k), 1e-4); % rat(pitchs(h, k)); % Get fraction
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

