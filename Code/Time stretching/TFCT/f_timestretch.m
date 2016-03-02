function y_v = f_timestretch(x_v, stretch)

N = length(x_v);
N_win = 1024;
Nfft = N_win;

% Pas d'avancement
R_n = floor(N_win / 4);
% Pas d'avancement après timestretching (pour reconstitution)
R_stretch_n = floor(R_n * stretch);

Nt = floor((N - N_win) / R_n);

ws_v = hann(N_win);
% Normalisation de la fenêtre
ws_v = ws_v./max(ola(ws_v .* ws_v, R_stretch_n, 30));

Xtilde_m = f_tfct (x_v, Nfft, N_win, R_n);
Xtilde_stretched_m = zeros(size(Xtilde_m));

puls = 2 * pi * R_n * (0:Nfft-1)' / Nfft;

phase_v = angle(Xtilde_m(:,1)); % Phase de X
prev_phase_v = phase_v;

for frame_n = 2 : Nt
    frame_v = Xtilde_m(:,frame_n);
    
    diff_phase_v = (angle(frame_v) - prev_phase_v) - puls;
    diff_phase_v = princarg(diff_phase_v);
    diff_phase_v = (diff_phase_v + puls) * stretch;

    phase_v = phase_v + diff_phase_v;

    Xtilde_stretched_m(:,frame_n) = abs(frame_v).*exp(1i*phase_v);

    prev_phase_v = angle(frame_v);
end

y_v = f_itfct (Xtilde_stretched_m, N_win, R_stretch_n, ws_v);

end