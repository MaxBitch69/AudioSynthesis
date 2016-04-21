function [y, y_test] = final(x, Fs, nbViolins, phaselock)
% Pitch Time Amplitude (PTA) algorithm developped by J. Pätynen
% With different techniques of phase locking and transient dectection

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
[amp, P] = max(output); % P is when max is ok (ie 1)
w = w./amp;

%% Algorithm
% Display progression
strf = 'Algorithm progression:';

% Work on every violin
for h = 1:nbViolins
    
    % Compute pitch with a random number following normal distribution
    pitch = normrnd(1, 0.005); % 1% of pitch modification
    
    % resample
    [p, q] = rat(pitch); % Get fraction
    signal = resample(x, q, p); % New time base vector - p/q, p/q times smaller
    
    N = length(signal);
    Nt = floor((N-Nw)/I); % Trames/FFT number
    
    %% Metro Hastings Sampling
    timeDifference = zeros(1, Nt); % Time Difference
    amplitudeModulation = zeros(1, Nt); % Amplitude modulation
    
    
    %% Time Difference - Metropolis-Hastings sampling
    %timeDifference = MetropolisHastings(0, 45, Nt); % mean = 0 ms
    timeDifferenceTemp = MetropolisHastings(0, 45, 2*Nt); % mean = 0 ms
                                                          % standard deviation = 45 ms
    
    % Low-frequency sampling, ie smoothing with Hanning filter
    %filt = 1/2*hanning(1000); % For great variations
    %filt = 1/30*hanning(1500); % Normal use
    filt = 1/10*hanning(floor(Fs/(5*I))); % cutoff frequency 5Hz
    
    %filt = 1/20*hanning(900); % Hanning filter - Trumpet

    timeDifferenceTemp = filter(filt, 1, timeDifferenceTemp); % Smooth on 1s
    %timeDifference = filter(filt, 1, timeDifference); % Smooth on 1s

    % Add offset or start to 0
    %timeDifference = timeDifference + normrnd(0, 30);%30*randn(1);
    timeDifference = timeDifferenceTemp(Nt+1:end); % take last part
    
    
    %% Amplitude modulation - Metropolis-Hastings sampling
    amplitudeModulation = abs(MetropolisHastings(1, 0.1, Nt)); % mean = 1
                                                               % standard
                                                               % deviation
                                                               % = 10%
    
    % Low-frequency sampling, ie smoothing
    % 5hz low frequency in the paper. Here, 1 pt is I, ie floor(Nw*0.25) pts
    % We want 5 Hz, ie Fs/N (frequence coupure pour hanning(N)), so N =
    % Fs/5 pts --> Fs/(5*I)
    len = floor(Fs/(5*I));
    filt = 1/len*hanning(len); % simple smoother, corresponding to 1s
    amplitudeModulation = filter(filt, 1, amplitudeModulation); % Smooth on 1s
    
    % Display results
    figure();
    plot(timeDifference);
    str = sprintf('Time Difference for violin n°%d', h);
    ylabel('Time Difference in ms');
    xlabel('Number of frames');
    title(str);
    
    
    %% STFT
    
    % Stretch factor
    stretch = pitch;

    % Initialisation
    puls = 2*pi*I*(0:Nfft-1)'/Nfft; % Canals' pulsations
    Xtilde_m = zeros(Nfft, Nt); % Matrix containing STFT of original signal
    Ytilde_m = zeros(Nfft, Nt); % Matrix containing STFT of synthesised signal
    
    deb = 1 + floor(timeDifference(1)*10^-3*Fs); % first time frame
    if deb <0 % If < 0 start from the beginning, else to tha right frame
       deb = 1;
    end
    fin = deb + Nw -1;
    Xtilde_m(:,1) = fft(signal(deb:fin).*w, Nfft); % 1st fft
    Ytilde_m(:,1) = Xtilde_m(:,1);
    
    % Parameters for time stretching
    phase = angle(Xtilde_m(:,1));
    former_phase = phase;

    % Parameters for transient detection
    transientDetected = 0;
    transientCenter = 0;
    
    % Parameters for phase locking
    pks_before = 0;
    former_pks = [];
    former_minAmp = [];
    
    for k=2:Nt  % Loop on timeframes
        
        % Display progression
        clc;
        str = sprintf('Violin n°%d, treatment progression: %.1f %%', h, 100*k/Nt);
        disp(strf);
        disp(str);
        
        %%% ANALYSIS
        % Time-base vector
        deb = (k-1)*I +1; % Beginning - x(n+kI)
        deb = deb + floor(timeDifference(k)*10^-3*Fs); % Time difference
        if deb <0 % Start to 1 if negative
            deb = 1;
        end
        fin = deb + Nw -1;
        if fin > N % Stay at the end if too much advance
            fin = N;
            deb = fin - Nw +1;
        end
        tx = signal(deb:fin).*w; % Timeframe
        
        % FFT
        Xtilde_m(:,k) = fft(tx,Nfft);
        
        % Get peaks, and minimum amplitude locations between peaks
        [pks, minAmp] = peakDetection(Xtilde_m(:,k), 3);

        transientCenter = 0;
%         % Look for transients
%         COG = computeCOG(signal(deb:fin), w, Nfft); % COG in % of the analysis window
%         if transientDetected
%             if COG < 0.04 % Threshold for being in the center of transient
%                 transientCenter = 1;
%             end
%         else
%             if COG > 0.3 % transientDetected
%                 transientDetected = 1;
%             end
%         end
        
        %% Main algorithm 
        if transientCenter % Reinitialize the phase
            % Parameters for time stretching
            Ytilde_m(:,k) = Xtilde_m(:,k);
            phase = angle(Xtilde_m(:,k));
            former_phase = phase;

            % Parameters for transient detection
            transientDetected = 0;
            transientCenter = 0;
    
            % Parameters for phase locking
            pks_before = 0;
            former_pks = [];
            former_minAmp = [];
            
        else % Compute new phase - phase lock techniques
            diff_time = I + floor(timeDifference(k)*10^-3*Fs)-floor(timeDifference(k-1)*10^-3*Fs); % Time interval

            [Ytilde_m(:,k), phase, former_phase, pks_before, former_pks, former_minAmp] = phaseUpdate(Xtilde_m(:,k), diff_time, stretch, I, former_phase, phase, pks, minAmp, pks_before, former_pks, former_minAmp, phaselock);
        end
        
        %%% SYNTHESIS
        % Time stretching
        R = floor(stretch*I);
        deb = (k-1)*R+1;
        fin = deb + Nw -1; % fin de trame
        
        % Reconstruction
        ys = real(ifft(Ytilde_m(:,k), 'symmetric')); % TFD inverse
        ys = ys.*ws; % pondération par la fenêtre de synthèse
        
        % Amplitude modulation
        %factorAmp = amplitudeModulation(h, k)/sum(amplitudeModulation(:,k)); % Normalisation
        %ys = factorAmp.*ys;
        
        y_test(deb:fin, h) = y_test(deb:fin, h) + ys; % Each signal - y_test: stereo
        % if 2 pitchs: soundsc(y_test, Fs)
        y(deb:fin)=y(deb:fin)+ys; % overlap add - sum of signals
    end
    
    Dm = STFT_consistency(y_test(:,h), Ytilde_m, R, P);
    
    strf = sprintf('%s\n%s\nDm = %fdB', strf, str, Dm);
end

display(strf);

end