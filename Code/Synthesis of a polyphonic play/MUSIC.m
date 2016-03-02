function [ output_args ] = MUSIC(x, n, K)
% Méthode Haute Résolution MUSIC
%   Affiche le pseudo-spectre MUSIC du signal x

%% Parameters
N = length(x);
l = N-n+1;

%% Correlation matrix - Question 3.2.1
a = x(1:n);
b = x(n:N);

X = hankel(a, b);
Rxx = 1/l*(X*X');

%% Method MUSIC - Question 3.2.6
%% Noise space estimation
[U1, Lambda, U2] = svd(Rxx);
W_ortho = U1(:, K+1:n); % espace bruit

%% Pseudo-spectre
ind_f = 1;
ind_delta = 1;
f_v = 0:0.001:1; % intervalle de recherche des fréquences
delta_v = -0.1:0.001:0.1; % intervalle de recherche des deltas
for f = f_v
    for delta = delta_v
        v = exp((delta+2*1i*pi*f)*(0:n-1)');
        P(ind_delta, ind_f) = 1/(norm(W_ortho'*v))^2;
        ind_delta = ind_delta + 1;
    end
    ind_delta=1;
    ind_f = ind_f+1;
end

%% Display pseudo-spectre
figure();
surf(f_v, delta_v, log10(P));

end

