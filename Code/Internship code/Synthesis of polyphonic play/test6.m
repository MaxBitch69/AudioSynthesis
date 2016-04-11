function [y, y_test] = test6(x, Fs, nbViolins)
% Time stretching using STFT of signal x, using stretching percentage
% "stretch".
%
% Work on each channel, independantly
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
    pitch = normrnd(1, 0.005); % 1% of pitch modification

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
    filt = 1/60*hanning(1500); %1/30*hanning(1500); % Hanning filter - Violin
    %filt = 1/20*hanning(900); % Hanning filter - Trumpet
    timeDifference = filter(filt, 1, timeDifference); % Smooth on 1s

    % Add offset or start to 0
    timeDifference = timeDifference + normrnd(0, 45);%30*randn(1);
    
%     % Time Difference for test
%     timeDifference = randn(1)*200*ones(1,Nt);
    
    % Amplitude modulation - Metropolis-Hastings sampling
    amplitudeModulation = abs(MetropolisHastings(1, 0.4, Nt)); % mean = 1
                                                               % standard
                                                               % deviation
                                                               % = 40%
                                                              
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
    title(str);


    %% STFT    
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
        if deb <0
            deb = 1;
        end
        fin = deb + Nw -1;
        tx = signal(deb:fin).*w; % Timeframe

        % FFT
        X = fft(tx,Nfft); 

        % Time stretching
%         stretch = pitch;
%         diff_phase = (angle(X) - former_phase) - puls;
%         diff_phase = mod(diff_phase + pi,-2*pi) + pi;
%         %diff_phase = (diff_phase + puls) * stretch;
%         diff_phase = (diff_phase + puls) * stretch * I/(I+(timeDifference(k)-timeDifference(k-1))*10^-3*Fs);

        % With former method
        stretch = pitch;
        diff_phase = (angle(X) - former_phase); % Phase difference
        diff_time = I + floor(timeDifference(k)*10^-3*Fs)-floor(timeDifference(k-1)*10^-3*Fs); % Time interval
        diff_phase = diff_phase - 2*pi*diff_time*(0:Nfft-1)'/Nfft; % Remove analysis window phase
        diff_phase = mod(diff_phase + pi,-2*pi) + pi;
        freq_inst = diff_phase/diff_time+2*pi*(0:Nfft-1)'/Nfft; % Freq instant
        diff_phase = freq_inst*stretch*I;
        
        phase = phase + diff_phase; 
        Y = abs(X).*exp(1i*phase);
        former_phase = angle(X);
                
        %%% SYNTHESIS
        % Time stretching
        R = floor(stretch*I);
        deb = (k-1)*R+1;
        fin = deb + Nw -1; % fin de trame

        % Reconstruction
        ys = real(ifft(Y, 'symmetric')); % TFD inverse
        ys = ys.*ws; % pondération par la fenêtre de synthèse
        
        % Amplitude modulation
        %factorAmp = amplitudeModulation(h, k)/sum(amplitudeModulation(:,k)); % Normalisation
        %ys = factorAmp.*ys;
        
        y_test(deb:fin, h) = y_test(deb:fin, h) + ys; % Each signal - y_test: stereo 
                                                      % if 2 pitchs: soundsc(y_test, Fs)
        y(deb:fin)=y(deb:fin)+ys; % overlap add - sum of signals
    end
    
    strf = sprintf('%s\n%s',strf, str);
end

% %% Amplitude modulation, using Metropolis-Hastings
% N = length(x);
% ampMod = zeros(nbViolins, N);
% 
% for h = 1:nbViolins
%     ampMod(h,:) = MetropolisHastings(1, 0.1, N); % mean = 1
%                                                   % standard
%                                                   % deviation
%                                                   % = 10%
%                                                               
%     % Low-frequency sampling, ie smoothing
%     % 5hz low frequency in the paper. Here, 1 pt is I, ie floor(Nw*0.25) pts
%     % We want 5 Hz, ie Fs/N (frequence coupure pour hanning(N)), so N =
%     % Fs/5 pts --> Fs/(5*I)
%     len = floor(Fs/5);
%     filt = 1/len*hanning(len); % simple smoother, corresponding to 1s
%     ampMod(h,:) = filter(filt, 1, ampMod(h,:)); % Smooth on 1s
% end
% 
% % Scale to 1
% for k = 1:N
%     ampMod(:,k) = ampMod(:,k)/sum(ampMod(:,k));
% end
% 
% % Amplitude modulation
% for h = 1:nbViolins
%    y_test(:, h) = y_test(:, h).*ampMod(h,:)'; % Each signal - y_test: stereo 
%    y = y+y_test(:, h);                                                   % if 2 pitchs: soundsc(y_test, Fs)    
% end
end