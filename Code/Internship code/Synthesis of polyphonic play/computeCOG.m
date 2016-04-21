function COG = computeCOG(x, w, Nfft)
% Function that computes t, based on time reassignement

%% Parameters
Nw = length(w); % Analysis window length and length of the execerpt

%% Th
mid = floor(Nw/2); % Start from the middle of the frame
T = [zeros(1, mid) (0:mid-1)]'; % Half ramp - Ones or 0 ?
Th = w.*T;

%% STFT's
X = fft(x.*w, Nfft);
Xt = fft(x.*Th, Nfft);

%% Reassigne time
% Phase derivate
diff_phase = real((Xt.*conj(X))./(abs(X).^2));

% Center Of Gravity - trapz does the integration on discrete values
cog_num = trapz((1:Nfft), -diff_phase.*(abs(X).^2));
cog_denum = trapz((1:Nfft), abs(X).^2);
COG = cog_num/cog_denum;

% Scale to window length
COG = COG/Nw;

end

