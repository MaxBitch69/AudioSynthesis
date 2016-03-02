function [a, phi] = LeastSquares(x, delta, f)
% Estime a et phi d'après x, delta et f avec la méthode 
% des Moindres Carrés

%% Question 3.2.4
%% Parameters
N = length(x);

%% Vandermonde matrix
T = (0:N-1)'; % vector column
A = delta' + 2*1i*pi*f'; % vector row

% Astuce de calcul pour la matrice de Vandermonde
v = T*A; % ln(V)
V = exp(v);

%% Computation of a and phi
alpha = pinv(V)*x;

a = abs(alpha);
phi = angle(alpha);

end

