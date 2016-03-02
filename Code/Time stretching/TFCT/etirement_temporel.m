% Analyse d'un signal à l'aide de la TFCT
% Transformée de Fourier à Court Terme
clear all;
close all;
clc;

%% Lecture du signal
[x,Fe] = audioread('signaux/signal40.wav');
x = x(1:10*Fe,1);
x = x(:); % ainsi x est un vect. colonne
% monovoie (voie gauche si stéréo)

%% Paramètres
affich = 0; % pour affichage du spectrogramme, 0 pour
% pour faire analyse/modif/synthèse sans affichage
% note: cf. spectrogram sous matlab
window = 0; % 0: fenêtre de Hann, else: fenêtre de Blackman-Harris
stretch = 2; % Facteur de streching

N = length(x); % longueur du signal
Nw = 1024; % longueur de la fenêtre
Nfft = 1024; % ordre de la TFD 

if window == 0
    % fenêtrage - Hann
    w = hanning(Nw); % définition de la fenetre d'analyse
    ws = w; % définition de la fenêtre de synthèse
    
    I = floor(Nw/8); % incrément sur les temps d'analyse,
    % appelé hop size, t_a=uR
    % Recouvrement de 88%
else
    w = blackmanharris(Nw); % définition de la fenetre d'analyse
    
    ws = rectwin(Nw); % définition de la fenêtre de synthèse: fenêtre
    % rectangulaire de même longueur
    
    I = floor(Nw/12); % incrément sur les temps d'analyse,
    % Blackman-Harris - Question 1.2.f
    % recouvrement: 92%
end

R = floor(stretch*I); % Etirement du signal: 2 fois plus long

Nt = fix((N-Nw)/I); % calcul du nombre de tfd à calculer
y = zeros(floor(N*stretch),1); % signal de synthèse

%% Condition de reconstruction parfaite
% Question 1.2.d et 1.2.f
h = w.*ws;
output = ola(h, I, 30); % recouvrement de 75% si Hann, 92% si Blackman

% Normalisation de la fenêtre
amp = max(output);
ws = ws./amp; % On nomralise la fenêtre de synthèse pour avoir une
% reconstruction = 1, peut aussi se faire sur la fenêtre
% d'analyse

if affich
    Xtilde = zeros(Nfft,Nt);
end

%% Initialisations
% Phases et pulsations
phase_X = zeros(Nfft, 1); % Phase de X
phase_Y = zeros(Nfft, 1); % Phase de Y
puls_instant = zeros(Nfft, Nt);

puls = 2*pi*(0:Nfft-1)'/Nfft; % Pulsations des canaux

%% Analyse-Synthèse par TFCT
for k=1:Nt;  % boucle sur les trames
    % Analyse
    deb = (k-1)*I +1; % début de trame - x(n+kR)
    fin = deb + Nw -1; % fin de trame
    tx = x(deb:fin).*w; % calcul de la trame
    X = fft(tx,Nfft); % tfd à l'instant b
    
    % Estimation des pulsations instantanées
    % Paramètres
    phase_ant_X = phase_X; % phase du tour d'avant/antérieure, 0 au premier tour
    phase_X = angle(X); % Phases pour tous les canaux du tour présent
    
    % Calcul des pulsations instantanées
    phase_instant = (phase_X-phase_ant_X);
    puls_instant = phase_instant./I + puls;
    
    % Synthèse
    deb = (k-1)*R +1; % début de trame - on écarte les instants de synthèse
    fin = deb + Nw -1; % fin de trame
    
    % Incrément de phase
    if k == 1 % premier tour
        phase_Y = phase_X;
    else
        phase_Y = phase_Y + (puls_instant - puls)*R;
    end
    
    % opérations de transformation - Etirement temporel

    Y = abs(X).*exp(1i*phase_Y); 

    % fin des opération de transformation
    
    % Reconstruction
    ys = real(ifft(Y)); % TFD inverse
    ys = ys(1:Nw).*ws; % pondération par la fenêtre de synthèse
    y(deb:fin)=y(deb:fin)+ys; % overlap add
end

%soundsc(y, Fe);

if affich
    figure();
    plot(y);
    title('Signal reconstruit');
    
    %% Affichage spectrogramme
    figure();
    freq = (0:Nfft/2)/Nfft*Fe; % échelle de fréquence
    b = ([0:Nt-1]*I+Nw/2)/Fe; % tps d'analyse (centre de la fenêtre)
    imagesc(b,freq,db(abs(Xtilde(1:L,:))))
    axis xy    
end