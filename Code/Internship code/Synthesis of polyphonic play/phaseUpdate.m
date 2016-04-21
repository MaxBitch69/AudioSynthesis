function [Ytilde_m, phase, former_phase, pks_before, former_pks, former_minAmp] = phaseUpdate(Xtilde_m, diff_time, stretch, I, former_phase, phase, pks, minAmp, pks_before, former_pks, former_minAmp, phaselock);
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
Nfft = length(Xtilde_m);
Ytilde_m = zeros(Nfft, 1);

%% Normal phase stretch
if strcmp(phaselock, '')
    
    % Time stretching
    % Instantaneous phase difference
    diff_phase = (angle(Xtilde_m) - former_phase); % Phase difference
    diff_phase = diff_phase - 2*pi*diff_time*(0:Nfft-1)'/Nfft; % Remove analysis window phase
    diff_phase = mod(diff_phase + pi,-2*pi) + pi;
    
    % Instantaneous pulsation
    puls_inst = (diff_phase+2*pi*diff_time*(0:Nfft-1)'/Nfft) / diff_time; % Freq instant
    diff_phase = puls_inst*stretch*I;
    
    phase = phase + diff_phase;
    
    Ytilde_m = abs(Xtilde_m).*exp(1i*phase);

    former_phase = angle(Xtilde_m);
        
%% Loose phase lock
elseif strcmp(phaselock, 'loose') % Loose phase lock
    
    % Time stretching
    % Instantaneous phase difference
    diff_phase = (angle(Xtilde_m) - former_phase); % Phase difference
    diff_phase = diff_phase - 2*pi*diff_time*(0:Nfft-1)'/Nfft; % Remove analysis window phase
    diff_phase = mod(diff_phase + pi,-2*pi) + pi;
    
    % Instantaneous pulsation
    puls_inst = (diff_phase+2*pi*diff_time*(0:Nfft-1)'/Nfft) / diff_time; % Freq instant
    diff_phase = puls_inst*stretch*I;
    
    phase = phase + diff_phase;
    
    Ytilde_m = abs(Xtilde_m).*exp(1i*phase);
    
    
    % Loose phase locking
    Y = zeros(Nfft, 1);
    
    % Add synthesised FFT's
    Y(1) = Ytilde_m(1) + Ytilde_m(2);
    Y(Nfft) = Ytilde_m(Nfft-1) + Ytilde_m(Nfft);
    Y(2:Nfft-1) = Ytilde_m(1:Nfft-2) + Ytilde_m(2:Nfft-1) + Ytilde_m(3:Nfft);
    
    % Get phase
    phaseloose = angle(Y);
    
    % Change phase for each channel k
    Ytilde_m = abs(Xtilde_m).*exp(1i*phaseloose); % Vector of synthesised channels
    % phase = phaselocked; % Update ? --> NO
    
    former_phase = angle(Xtilde_m);
    
  
    

%% Rigid && scaled phase locks
elseif strcmp(phaselock(1:5), 'rigid') || strcmp(phaselock(1:6), 'scaled')

    if isempty(pks) % If no peaks, usual method
        % Time stretching
        % Instantaneous phase difference
        diff_phase = (angle(Xtilde_m) - former_phase); % Phase difference
        diff_phase = diff_phase - 2*pi*diff_time*(0:Nfft-1)'/Nfft; % Remove analysis window phase
        diff_phase = mod(diff_phase + pi,-2*pi) + pi;
        
        % Instantaneous pulsation
        puls_inst = (diff_phase+2*pi*diff_time*(0:Nfft-1)'/Nfft) / diff_time; % Freq instant
        diff_phase = puls_inst*stretch*I;
        
        phase = phase + diff_phase;
        
        Ytilde_m = abs(Xtilde_m).*exp(1i*phase);
        
        former_phase = angle(Xtilde_m);
        
        pks_before = 0;
        former_pks = [];
        
    else % If peaks
        if strcmp(phaselock(end-9:end), 'WithinPeak') || ~pks_before % update within the same peaks
    
            % Time stretching for peaks
            diff_phase = angle(Xtilde_m(pks)) - former_phase(pks); % Phase difference
            diff_phase = diff_phase - 2*pi*diff_time*(pks)/Nfft; % Remove analysis window phase
            diff_phase = mod(diff_phase + pi,-2*pi) + pi;
        
            puls_inst = diff_phase/diff_time+2*pi*(pks)/Nfft; % Instantaneous pulsation
            diff_phase = puls_inst*stretch*I;
        
            phase(pks) = phase(pks) + diff_phase;
  
            pks_before = 1;
            former_pks = pks;
            former_minAmp = minAmp;
            
        elseif strcmp(phaselock(end-8:end), 'InterPeak') % Can change peaks
            
            % Looking for associated peaks - previous_pks VS pks
            previous_pks = zeros(length(pks),1);
            for k = 1:length(pks)
                [~, ind] = max(former_minAmp>pks(k));
                previous_pks(k) = former_pks(ind-1);
            end
            
            % Time stretching for peaks
            diff_phase = angle(Xtilde_m(pks)) - former_phase(previous_pks); % Phase difference
            diff_phase = diff_phase - 2*pi*diff_time*(pks)/Nfft; % Remove analysis window phase
            diff_phase = mod(diff_phase + pi,-2*pi) + pi;
        
            puls_inst = diff_phase/diff_time+2*pi*(pks)/Nfft; % Instantaneous pulsation
            diff_phase = puls_inst*stretch*I;
        
            phase(pks) = phase(previous_pks) + diff_phase;
            
            pks_before = 1;
            former_pks = pks;
            former_minAmp = minAmp;
            
        else
            error('Wrong phase lock technique: possible choices: '''', ''loose'', ''rigidWithinPeak'', ''rigidInterPeak'', ''scaledWithinPeak'', ''scaledInterPeak''');
        
        end
    
    
        if strcmp(phaselock(1:5), 'rigid') % Rigid update

            % Teta and phasor Z
            teta = phase(pks) - angle(Xtilde_m(pks));
            Z = exp(1i.*teta);
        
            % Update areas aroud peaks
            Ytilde_m(1:minAmp(2)) = Z(1).*Xtilde_m(1:minAmp(2));
            for k = 2:length(minAmp)-1
                Ytilde_m(minAmp(k):minAmp(k+1)) = Z(k).*Xtilde_m(minAmp(k):minAmp(k+1));
            end
                
            phase = angle(Ytilde_m); % Update phase for the other case
            former_phase = angle(Xtilde_m); % Update former phase
        
        
        else % Scaled
            beta = 2/3+stretch/3; % when = 1: unity scale

            % Update areas aroud peaks
            % Unwrap phases around peaks
            diff_phase = angle(Xtilde_m(1:minAmp(2))) - angle(Xtilde_m(pks(1))); % Phase difference
            diff_phase = mod(diff_phase + pi,-2*pi) + pi;
                
            % Scaled phase lock in the area
            phase(1:minAmp(2)) = phase(pks(1)) + beta*(diff_phase);
                
            for k = 2:length(minAmp)-1
                % Unwrap phases around peaks
                diff_phase = angle(Xtilde_m(minAmp(k):minAmp(k+1))) - angle(Xtilde_m(pks(k))); % Phase difference
                diff_phase = mod(diff_phase + pi,-2*pi) + pi;
                
                % Scaled phase lock in the area
                phase(minAmp(k):minAmp(k+1)) = phase(pks(k)) + beta*(diff_phase);
            end
        
            Ytilde_m = abs(Xtilde_m).*exp(1i*phase);
            former_phase = angle(Xtilde_m);
        end
    end
    
    
%% Wrong entry
else
    error('Wrong phase lock technique: possible choices: '''', ''loose'', ''rigidWithinPeak'', ''rigidInterPeak'', ''scaledWithinPeak'', ''scaledInterPeak''');
end

end
