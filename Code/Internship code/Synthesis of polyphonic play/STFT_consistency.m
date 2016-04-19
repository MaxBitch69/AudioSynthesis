function Dm = STFT_consistency(y, Ytilde_m, R, P)
% Computes the consistency of the synthesised STFT - shows the
% quality of the reconstructed signal
%
% Inputs:
%       - y: synthesised signal
%       - Ytilde_m: modified STFT
%       - overlap: overlap in %
%       - P: to exclude the P first and P last few frames to avoid taking 
%       into account errors due to missing overlapped segments in the 
%       resynthesis formula.
%
% Output:
%       - Dm = error between the two STFT in dB

%% Parameters
Nfft = size(Ytilde_m, 1); % Nfft, number of channels
Nt = size(Ytilde_m, 2); % Number of frames

%% y's STFT
Nw = Nfft;
w = hanning(Nw);
Y = zeros(Nfft, Nt);

for k=1:Nt
    deb = (k-1)*R +1; % frame's beginning - x(n+kI)
    fin = deb + Nw -1; % frame's end
    ty = y(deb:fin).*w; % Timeframe

    Y(:,k) = fft(ty, Nfft); % FFT on the timeframe    
end 

%% Compute Dm
num = 0;
denum = 0;
for k = 1+P:Nt-P
    vect = sum((abs(Y(:,k)) - abs(Ytilde_m(:,k))).^2);
    num = num + vect;
    
    denum = denum + sum(abs(Ytilde_m(:,k)).^2);
end

Dm = 20*log(num/denum);

end
