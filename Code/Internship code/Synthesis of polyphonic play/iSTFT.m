function y = iSTFT(Ytilde_m, Nw, overlap, stretch)
% Computes signal Ytilde_m's Inverse Short Term Fourier Transform
% and gets synthesised signal y
%
% Inputs:
%       - Ytilde_m: modified STFT
%       - Nw = analysis window's length. By default, Nfft = Nw, 
%       and the synthesised and analysis windows have the same length
%       - overlap: min of 75% advised
%       - stretch: stretch factor for the synthesis marks (Ytilde_m is
%       already modified)
%
% Output:
%       - synthesised signal


%% Parameters
Nt = size(Ytilde_m, 2); % Number of frames

overlap = 1-overlap;
I = floor(Nw*overlap); % Analysis hop size in points
R = floor(stretch*I); % Synthesis hop size

y = zeros(floor(Nt*R)+Nw, 1); % Synthesised signal

%% Windowing and reconstruction hypothesis
ws = hanning(Nw); % Synthesis window
w = ws; % Analysis window

% Check for perfect reconstruction - ie h = w.*ws == 1
h = w.*ws;
output = ola(h, I, 30); % Check reconstruction

% window's normalisation - w.*ws == 1
amp = max(output); 
ws = ws./amp;

%% iSTFT
for k=2:Nt;  % Loop on timeframes
    deb = (k-1)*R+1; % frame's begin - x(n+kR)
    fin = deb + Nw -1; % frame's end

    % Reconstruction
    ys = real(ifft(Ytilde_m(:,k), 'symmetric')); % iFFT
    ys = ys.*ws; % synthesis window
        
    y(deb:fin)=y(deb:fin)+ys; % overlap add - sum of signals
end

y = y(1:floor(Nt*R)); % Resizing vector

end

