clear all;
close all;
clc;

%% Reading signal
[x, Fs] = audioread('../Audio/Violin/Vibrato_G3.aif');

% For stereo sounds
if size(x, 2) == 2
    x = (x(:,1)+x(:,2))/2; % For stereo sounds
end

%y = timestretch(x, 2);
y = pitchshift(x, 1.5);
z = pitchshift(x, 2);
w = pitchshift(x, 1.3);

dimension = min([length(y) length(w) length(z)]); 
u = y(1:dimension) + z(1:dimension) + w(1:dimension);

%locs = onsetDetection(x, Fs)