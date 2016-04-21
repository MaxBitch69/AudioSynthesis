function [y, y_test] = PTAPhaseLock(x, Fs, nbViolins, phaselock)
% Pitch Time Amplitude (PTA) algorithm developped by J. Pätynen
% With different techniques of phase locking

close all

%% Parameters
N = length(x); % Signal's duration
Nw = 2048; % Excerpts of signal
Nfft = 2048; % Precision of FFT

overlap = 0.25; % overlap in %, here 75%
I = floor(Nw*overlap); % Hop size in points

y = zeros(N, 1); % Synthesised signal
y_test = zeros(N, nbViolins); % Different channels

% Phase lock technique
%phaselock = 'scaled'; %loose'; % Loose phase lock

%% Windowing
w = hanning(Nw); % Analysis window
ws = w; % Synthesis window

% Check for perfect reconstruction - ie h = w.*ws == 1
h = w.*ws;
output = ola(h, I, 30); % Check reconstruction

% window's normalisation - w.*ws == 1
[amp, P] = max(output); % P is when max is ok (ie 1)
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
    %timeDifference = MetropolisHastings(0, 45, Nt); % mean = 0 ms
    timeDifferenceTemp = MetropolisHastings(0, 45, 2*Nt); % mean = 0 ms
    
    % standard deviation = 45 ms
    
    % Low-frequency sampling, ie smoothing
    % To see well --> filt = 1/2*hanning(1000); %1/30*hanning(1500); % Hanning filter - Violin
    %filt = 1/2*hanning(1000);
    filt = 1/30*hanning(1500);
    %filt = 1/10*hanning(floor(Fs/(5*I))); % 5Hz
    
    %filt = 1/20*hanning(900); % Hanning filter - Trumpet

    timeDifferenceTemp = filter(filt, 1, timeDifferenceTemp); % Smooth on 1s
    %timeDifference = filter(filt, 1, timeDifference); % Smooth on 1s

    % Add offset or start to 0
    %timeDifference = timeDifference + normrnd(0, 30);%30*randn(1);
    timeDifference = timeDifferenceTemp(Nt+1:end); % take last part
    
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
    Xtilde_m = zeros(Nfft, Nt); % Matrix containing STFT of original signal
    Ytilde_m = zeros(Nfft, Nt); % Matrix containing STFT of synthesised signal
    
    
    % MODIF ICI
    %Xtilde_m(:,1) = fft(x(1:Nw), Nfft); % 1st fft
    deb = 1 + floor(timeDifference(1)*10^-3*Fs); % Time difference
    if deb <0
       deb = 1;
    end
    fin = deb + Nw -1;

    Xtilde_m(:,1) = fft(x(deb:fin), Nfft); % 1st fft
    Ytilde_m(:,1) = Xtilde_m(:,1);
    
    % Parameters for time stretching
    phase = angle(Xtilde_m(:,1));
    former_phase = phase;
    
    emptyPeaks = 1; % for scaled phase locking
    
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
        Xtilde_m(:,k) = fft(tx,Nfft);
        
        %% Phase lock techniques
        if strcmp(phaselock, '') % Normal phase vocoder approach
            % Time stretching
            stretch = pitch;
            diff_phase = (angle(Xtilde_m(:,k)) - former_phase); % Phase difference
            diff_time = I + floor(timeDifference(k)*10^-3*Fs)-floor(timeDifference(k-1)*10^-3*Fs); % Time interval
            diff_phase = diff_phase - 2*pi*diff_time*(0:Nfft-1)'/Nfft; % Remove analysis window phase
            diff_phase = mod(diff_phase + pi,-2*pi) + pi;
            freq_inst = diff_phase/diff_time+2*pi*(0:Nfft-1)'/Nfft; % Freq instant
            diff_phase = freq_inst*stretch*I;
            
            phase = phase + diff_phase;
            
            Ytilde_m(:,k) = abs(Xtilde_m(:,k)).*exp(1i*phase);
            
        elseif strcmp(phaselock, 'loose') % Loose phase lock
            stretch = pitch;
            diff_phase = (angle(Xtilde_m(:,k)) - former_phase); % Phase difference
            diff_time = I + floor(timeDifference(k)*10^-3*Fs)-floor(timeDifference(k-1)*10^-3*Fs); % Time interval
            diff_phase = diff_phase - 2*pi*diff_time*(0:Nfft-1)'/Nfft; % Remove analysis window phase
            diff_phase = mod(diff_phase + pi,-2*pi) + pi;
            freq_inst = diff_phase/diff_time+2*pi*(0:Nfft-1)'/Nfft; % Freq instant
            diff_phase = freq_inst*stretch*I;
            
            phase = phase + diff_phase;
            
            Ytilde_m(:,k) = abs(Xtilde_m(:,k)).*exp(1i*phase);
            
            Y = zeros(Nfft, 1);
            
            % Add synthesised FFT
            Y(1) = Ytilde_m(1,k) + Ytilde_m(2,k);
            Y(Nfft) = Ytilde_m(Nfft-1,k) + Ytilde_m(Nfft,k);
            Y(2:Nfft-1) = Ytilde_m(1:Nfft-2,k) + Ytilde_m(2:Nfft-1,k) + Ytilde_m(3:Nfft,k);
            
            % Get phase
            phaselocked = angle(Y);
            
            % Change phase for channel k
            Ytilde_m(:,k) = abs(Xtilde_m(:,k)).*exp(1i*phaselocked); %abs(Xtilde_m(:,k)).*exp(1i*phaselocked); % Vector of synthesised channels
            
            phase = phaselocked; % Update ? 
            
        elseif strcmp(phaselock, 'rigid') % Rigid phase lock

            % Find peaks
            [~, locs] = findpeaks(abs(Xtilde_m(:,k)),'MinPeakDistance',3); % Peaks spaced by 4
            
            if isempty(locs)
                stretch = pitch;
                diff_phase = (angle(Xtilde_m(:,k)) - former_phase); % Phase difference
                diff_time = I + floor(timeDifference(k)*10^-3*Fs)-floor(timeDifference(k-1)*10^-3*Fs); % Time interval
                diff_phase = diff_phase - 2*pi*diff_time*(0:Nfft-1)'/Nfft; % Remove analysis window phase
                diff_phase = mod(diff_phase + pi,-2*pi) + pi;
                freq_inst = diff_phase/diff_time+2*pi*(0:Nfft-1)'/Nfft; % Freq instant
                diff_phase = freq_inst*stretch*I;
                
                phase = phase + diff_phase;
                
                Ytilde_m(:,k) = abs(Xtilde_m(:,k)).*exp(1i*phase);
            else            
                % Time stretching for peaks
                stretch = pitch;
                diff_phase = (angle(Xtilde_m(locs,k)) - former_phase(locs)); % Phase difference
                diff_time = I + floor(timeDifference(k)*10^-3*Fs)-floor(timeDifference(k-1)*10^-3*Fs); % Time interval
                diff_phase = diff_phase - 2*pi*diff_time*(locs)/Nfft; % Remove analysis window phase
                diff_phase = mod(diff_phase + pi,-2*pi) + pi;
                freq_inst = diff_phase/diff_time+2*pi*(locs)/Nfft; % Freq instant
                diff_phase = freq_inst*stretch*I;
                
                phase(locs) = phase(locs) + diff_phase;
                
                % Teta and phasor Z
%                 beta = 2/3+pitch/3; % when = 1: unity scale
%                 teta = beta*(phase(locs) - angle(Xtilde_m(locs,k)));
                teta = phase(locs) - angle(Xtilde_m(locs,k));

                Z = exp(1i.*teta);
                
                min_inf = 1;
                l = 1;
                
                while min_inf < locs(length(locs)-1) % Last peak
                    % Define area for peak l
                    [~, min_sup] = min(abs(Xtilde_m(locs(l):locs(l+1),k)));
                    min_sup = min_sup + locs(l);
                    
                    % Rigid phase lock in the area
                    Ytilde_m(min_inf:min_sup,k) = Z(l).*Xtilde_m(min_inf:min_sup,k);

                    % Increment variables
                    min_inf = min_sup;
                    l = l+1;
                end
                min_sup = length(Xtilde_m(:,k));
                Ytilde_m(min_inf:min_sup,k) = Z(l).*Xtilde_m(min_inf:min_sup,k);
            
                phase = angle(Ytilde_m(:,k)); % Update phase for the other case
            end
        elseif strcmp(phaselock, 'scaled')
            % Find peaks
            [~, locs] = findpeaks(abs(Xtilde_m(:,k)),'MinPeakDistance',3); % Peaks spaced by 4
            
            if isempty(locs)
                %strf = sprintf('%s\nEMPTY ! at %.1f %%', strf, 100*k/Nt);
                stretch = pitch;
                diff_phase = (angle(Xtilde_m(:,k)) - former_phase); % Phase difference
                diff_time = I + floor(timeDifference(k)*10^-3*Fs)-floor(timeDifference(k-1)*10^-3*Fs); % Time interval
                diff_phase = diff_phase - 2*pi*diff_time*(0:Nfft-1)'/Nfft; % Remove analysis window phase
                diff_phase = mod(diff_phase + pi,-2*pi) + pi;
                freq_inst = diff_phase/diff_time+2*pi*(0:Nfft-1)'/Nfft; % Freq instant
                diff_phase = freq_inst*stretch*I;
                
                phase = phase + diff_phase;
                
                Ytilde_m(:,k) = abs(Xtilde_m(:,k)).*exp(1i*phase);
                emptyPeaks = 1;
            else          
                if emptyPeaks == 1 % No peaks before
                    % Time stretching for peaks
                    stretch = pitch;
                    diff_phase = (angle(Xtilde_m(locs,k)) - former_phase(locs)); % Phase difference
                    diff_time = I + floor(timeDifference(k)*10^-3*Fs)-floor(timeDifference(k-1)*10^-3*Fs); % Time interval
                    diff_phase = diff_phase - 2*pi*diff_time*(locs)/Nfft; % Remove analysis window phase
                    diff_phase = mod(diff_phase + pi,-2*pi) + pi;
                    freq_inst = diff_phase/diff_time+2*pi*(locs)/Nfft; % Freq instant
                    diff_phase = freq_inst*stretch*I;
                
                    phase(locs) = phase(locs) + diff_phase;
                
                    % Teta and phasor Z
                    teta = phase(locs) - angle(Xtilde_m(locs,k));
                    Z = exp(1i.*teta);
                
                    min_inf = 1;
                    l = 1;
                    min_sups = zeros(length(locs), 1); % Vector of min_sup
                    
                    while min_inf < locs(length(locs)-1) % Last peak
                        % Define area for peak l
                        [~, min_sup] = min(abs(Xtilde_m(locs(l):locs(l+1),k)));
                        min_sup = min_sup + locs(l);
                        min_sups(l) = min_sup;
                        
                        % Rigid phase lock in the area
                        Ytilde_m(min_inf:min_sup,k) = Z(l).*Xtilde_m(min_inf:min_sup,k);

                        % Increment variables
                        min_inf = min_sup;
                        l = l+1;
                    end
                    min_sup = length(Xtilde_m(:,k));
                    Ytilde_m(min_inf:min_sup,k) = Z(l).*Xtilde_m(min_inf:min_sup,k);
                    min_sups(l) = min_sup;
                    
                    phase = angle(Ytilde_m(:,k)); % Update phase for the other case
                    emptyPeaks = 0;
                    formerlocs = locs;
                else
                    % Time stretching for peaks
                    stretch = pitch;

                    % Find former locs
                    previouspeak = zeros(length(locs), 1);
                    for j = 1:length(locs)
                        [~, previouspeak_ind] = min(abs(locs(j)-min_sups)); % Check which is the closest sup
                        previouspeak(j) = formerlocs(previouspeak_ind); % Belongs to the same former peak
                    end
                    diff_phase = angle(Xtilde_m(locs,k)) - former_phase(previouspeak); % Phase difference
                    diff_time = I + floor(timeDifference(k)*10^-3*Fs)-floor(timeDifference(k-1)*10^-3*Fs); % Time interval
                    diff_phase = diff_phase - 2*pi*diff_time*(locs)/Nfft; % Remove analysis window phase
                    diff_phase = mod(diff_phase + pi,-2*pi) + pi;
                    freq_inst = diff_phase/diff_time+2*pi*(locs)/Nfft; % Freq instant
                    diff_phase = freq_inst*stretch*I;
                
                    phase(locs) = phase(previouspeak) + diff_phase;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BETA
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCALED
                    min_inf = 1;
                    l = 1;
                    beta = 2/3+pitch/3; % when = 1: unity scale                
                    min_sups = zeros(length(locs), 1); % Vector of min_sup
                    
                    while min_inf < locs(length(locs)-1) % Last peak
                        % Define area for peak l
                        [~, min_sup] = min(abs(Xtilde_m(locs(l):locs(l+1),k)));
                        min_sup = min_sup + locs(l);
                        min_sups(l) = min_sup;
                        
                        % Unwrap phases around peaks
                        diff_phase = angle(Xtilde_m(min_inf:min_sup, k)) - angle(Xtilde_m(locs(l), k)); % Phase difference
                        diff_phase = mod(diff_phase + pi,-2*pi) + pi;
                        
                        % Scaled phase lock in the area
                        phase(min_inf:min_sup) = phase(locs(l)) + beta*(diff_phase);
                        Ytilde_m(min_inf:min_sup,k) = abs(Xtilde_m(min_inf:min_sup,k)).*exp(1i*phase(min_inf:min_sup));

                        % Increment variables
                        min_inf = min_sup;
                        l = l+1;
                    end
                    min_sup = length(Xtilde_m(:,k));
                    % Unwrap phases around peaks
                    diff_phase = angle(Xtilde_m(min_inf:min_sup, k)) - angle(Xtilde_m(locs(l), k)); % Phase difference
                    diff_phase = mod(diff_phase + pi,-2*pi) + pi;
                        
                    % Scaled phase lock in the area
                    phase(min_inf:min_sup) = phase(locs(l)) + beta*(diff_phase);
                    Ytilde_m(min_inf:min_sup,k) = abs(Xtilde_m(min_inf:min_sup,k)).*exp(1i*phase(min_inf:min_sup));
                    
                    min_sups(l) = min_sup;
                    
                    phase = angle(Ytilde_m(:,k)); % Update phase for the other case
                    emptyPeaks = 0;
                    formerlocs = locs;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
%                     % Teta and phasor Z
%                     teta = phase(locs) - angle(Xtilde_m(locs,k));
%                     Z = exp(1i.*teta);
%                 
%                     min_inf = 1;
%                     l = 1;
%                     min_sups = zeros(length(locs), 1); % Vector of min_sup
%                     
%                     while min_inf < locs(length(locs)-1) % Last peak
%                         % Define area for peak l
%                         [~, min_sup] = min(abs(Xtilde_m(locs(l):locs(l+1),k)));
%                         min_sup = min_sup + locs(l);
%                         min_sups(l) = min_sup;
%                         
%                         % Rigid phase lock in the area
%                         Ytilde_m(min_inf:min_sup,k) = Z(l).*Xtilde_m(min_inf:min_sup,k);
% 
%                         % Increment variables
%                         min_inf = min_sup;
%                         l = l+1;
%                     end
%                     min_sup = length(Xtilde_m(:,k));
%                     Ytilde_m(min_inf:min_sup,k) = Z(l).*Xtilde_m(min_inf:min_sup,k);
%                     min_sups(l) = min_sup;
%                     
%                     phase = angle(Ytilde_m(:,k)); % Update phase for the other case
%                     emptyPeaks = 0;
%                     formerlocs = locs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                end

        end
        
        former_phase = angle(Xtilde_m(:,k));
        
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