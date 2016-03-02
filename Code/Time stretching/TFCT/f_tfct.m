function Xtilde_m = f_tfct (x_v, Nfft, L_n, R)

if size(x_v,2) >= 2 
    x_v = (x_v(:,1) + x_v(:,2)) / 2;
end

N = length(x_v); % longueur du signal
w_v = hanning(L_n); % définition de la fenetre d'analyse

Nt = fix((N-L_n)/R); % calcul du nombre de tfd à calculer

Xtilde_m = zeros(Nfft,Nt);

for k=1:Nt;  % boucle sur les trames
   deb = (k-1)*R +1; % début de trame
   fin = deb + L_n -1; % fin de trame
   tx = x_v(deb:fin).*w_v; % calcul de la trame  
   Xtilde_m(:,k)=fft(tx,Nfft);
end

end