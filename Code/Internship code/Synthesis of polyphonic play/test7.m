function [y, y_test] = test7(x, Fs, nbViolins)
% Time stretching using STFT of signal x, using stretching percentage
% "stretch".
%
% Work on each channel, independantly
% Try to do as in the article, but don't understand, where occurs the pitch
% shift ? here, sample and resample, so nothing has changed !
% Good result though
%
% Pitchs = list of new pitchs for pitch shifting, decimals

close all

%% Parameters
N = length(x); % Signal's duration
Nw = 2048; % Excerpts of signal
Nfft = 2048; % Precision of FFT

overlap = 0.25; % overlap in %, here 75%
I = floor(Nw*overlap); % Hop size in points

y = zeros(N, 1); % Synthesised signal
y_test = zeros(N, nbViolins); % Different channels

%% Windowing
w = hanning(Nw); % Analysis window
ws = w; % Synthesis window

% Check for perfect reconstruction - ie h = w.*ws == 1
h = w.*ws;
output = ola(h, I, 30); % Check reconstruction

% window's normalisation - w.*ws == 1
amp = max(output);
w = w./amp; 


% Display progression
strf = 'Algorithm progression:';

%% Work on every violin
for h = 1:nbViolins

    % Compute pitch with a random number following normal distribution
    pitch = normrnd(1, 0.1); % 1% of pitch modification

    % resample
    [p, q] = rat(pitch); % Get fraction
    signal = resample(x, q, p); % New time base vector - p/q, p/q times smaller

    N = length(signal);
    Nt = floor((N-Nw)/I); % Trames/FFT number

    %% Metro Hastings Sampling
    timeDifference = zeros(1, Nt); % Time Difference
    amplitudeModulation = zeros(1, Nt); % Amplitude modulation

    % Time Difference - Metropolis-Hastings sampling
    timeDifference = MetropolisHastings(0, 45, Nt); % mean = 0 ms
                                                    % standard deviation = 45 ms

    % Low-frequency sampling, ie smoothing
    filt = 1/30*hanning(2500); % Hanning filter - Violin
    %filt = 1/20*hanning(900); % Hanning filter - Trumpet
    timeDifference = filter(filt, 1, timeDifference); % Smooth on 1s

%     % Time Difference for test
%     timeDifference = randn(1)*200*ones(1,Nt);
    
    % Amplitude modulation - Metropolis-Hastings sampling
    amplitudeModulation = abs(MetropolisHastings(1, 0.4, Nt)); % mean = 1
                                                               % standard
                                                               % deviation
                                                               % = 40%
                                                              
    % Low-frequency sampling, ie smoothing
    filt = 1/30*hanning(600); % simple smoother, corresponding to 1s
    amplitudeModulation = filter(filt, 1, amplitudeModulation); % Smooth on 1s

    % Display results
    figure();
    plot(timeDifference);
    str = sprintf('Time Difference for violin n°%d', h);
    title(str);


    %% STFT    
    y1 = zeros(N, 1);
    
    % Initialisation
    puls = 2*pi*I*(0:Nfft-1)'/Nfft; % Canals' pulsations
    Xtilde_m = zeros(Nfft, Nt); % Matrix containing fft
    Xtilde_m(:,1) = fft(x(1:Nw), Nfft); % 1st fft

    % Parameters for time stretching
    phase = angle(Xtilde_m(:,1));
    former_phase = phase;

    for k=2:Nt-20  % Loop on timeframes
        % Display progression
        clc;
        str = sprintf('Violin n°%d, treatment progression: %.1f %%', h, 100*k/Nt);
        disp(strf);
        disp(str);

        %%% ANALYSIS
        % Time-base vector
        deb = (k-1)*I +1; % Beginning - x(n+kI)
        deb = deb + floor(timeDifference(k)*10^-3*Fs); % Time difference
        fin = deb + Nw -1;
        tx = signal(deb:fin).*w; % Timeframe

        % FFT
        X = fft(tx,Nfft); 
        Y = X;
%         % Time stretching
%         stretch = pitch;
%         diff_phase = (angle(X) - former_phase) - puls;
%         diff_phase = mod(diff_phase + pi,-2*pi) + pi;
%         diff_phase = (diff_phase + puls) * stretch;
% 
%         phase = phase + diff_phase;
%         Y = abs(X).*exp(1i*phase);
%         former_phase = angle(X);
                
        %%% SYNTHESIS
%         % Time stretching
%         R = floor(stretch*I);
%         deb = (k-1)*R+1;
%         fin = deb + Nw -1; % fin de trame

        % Reconstruction
        ys = real(ifft(Y, 'symmetric')); % TFD inverse
        ys = ys.*ws; % pondération par la fenêtre de synthèse
        
        % Amplitude modulation
        %factorAmp = amplitudeModulation(h, k)/sum(amplitudeModulation(:,k)); % Normalisation
        %ys = factorAmp.*ys;

        deb = (k-1)*I +1; % Beginning - x(n+kI)
        fin = deb + Nw -1;
        
        y1(deb:fin) = y1(deb:fin)+ys; % overlap add - sum of signals
    end
    
    y2 = resample(y1, p, q);
    y_test(:,h) = y2(1:length(y_test(:,h)));
    y = y + y_test(:,h); % Each signal - y_test: stereo 
                                                      % if 2 pitchs: soundsc(y_test, Fs)

    strf = sprintf('%s\n%s',strf, str);
end
end