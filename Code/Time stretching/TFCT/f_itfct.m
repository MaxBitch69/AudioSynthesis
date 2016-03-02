function x_v = f_itfct (Xtilde_m, L_n, R, w_v)

Nt = size(Xtilde_m, 2);
N = Nt * R + L_n;
x_v = zeros(N,1);

for k=1:Nt;  % boucle sur les trames
   deb = (k-1)*R +1; % début de trame
   fin = deb + L_n -1; % fin de trame
   X_frame_v = Xtilde_m(:,k);
   x_temp_v = real(ifft(X_frame_v, 'symmetric'));
   x_temp_v = x_temp_v .* w_v;
   x_v(deb:fin) = x_v(deb:fin) + x_temp_v;
end

end