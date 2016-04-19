function [periodMarks, pitchs, notes] = analysePitch(x, Fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
N = length(x);
analysisLength = 10*10^-3*Fs; % analysis frame = 25ms

periodMarks(1) = 1; % Analysis marks, begin to 1st sample
pitchs(1) = 10*10^-3*Fs; % pitchs found
notes(1) = 1;

%% Analyse signal
n = 1;
q = 1; % For onsets

while (periodMarks(n) + 2.5*analysisLength) < N % Until end of file
    % Get period
    s = x(periodMarks(n):periodMarks(n) + floor(2.5*analysisLength)); % Excerpt to analyse
    P = periode(s, Fs); % Get period
    
    % Compute data
    n = n+1;
    periodMarks(n) = periodMarks(n-1) + P; % Analysis marks, on lobes
    pitchs(n) = Fs/P; % Get pitch
    if (pitchs(n)/pitchs(n-1)) > (2^(1/12)) % Note has changed
        q = q+1;
        notes(q) = periodMarks(n); % Get onset
    end
    
    % Display progression
%     clc;
    str = sprintf('Treatment progression: %.1f %%', 100*(periodMarks(n) + 2.5*analysisLength)/N);
    disp(str);
end
    
%% Display results
figure();
plot(periodMarks/Fs, pitchs);

end
