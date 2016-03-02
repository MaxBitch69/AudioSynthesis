function [delta, f] = ESPRIT(x, n, K)
% Methode ESPRIT: calcule l'espace signal de dimension K du signal x
% et en extrait les fr�quences fk et les facteurs d'amortissement delta_k

%% Parameters
N = length(x);
l = N-n+1;

%% Correlation matrix - Question 3.2.1
a = x(1:n);
b = x(n:N);

X = hankel(a, b);
Rxx = 1/l*(X*X');

%% Signal space estimation - Question 3.2.2
[U1, Lambda, U2] = svd(Rxx); % On d�compose la matrice de corr�lation sur la base 
                             % des vecteurs propres/ valeurs propres.
                             % Lambda repr�sente la matrice des valeurs
                             % propres diagonale par ordre d�croissant
                             % et U1 et U2 les matrices des vecteurs
                             % propres associ�s aux valeurs propres
                             % Ici, la matrice �tant carr�e, on a
                             % clairement U1 = U2 (cas o� ils ne sont pas
                             % �gaux, quand il n'y pas de bruit)
                             
% on a U1 qui est, � peu de choses pr�s une base de la matrice de
% corr�lation. L'espace signal W est la partie de cette matrice de
% corr�lation de dimension K qui a les valeurs propres les plus
% importantes. Ainsi ici, on a K = 2, et svd trie les valeurs propres par
% ordre d�croissant. On prendra donc comme base de l'espace signal les 2
% premiers vecteurs de U1, qui sont associ�s aux 2 valeurs propres les plus
% importantes de la matrice de corr�lation.

%% Method ESPRIT
%% Question 3.2.3
% Estimation des fr�quences et des facteurs d'amortissement
W = U1(:, 1:K);
W_bas = W(1:n-1, :);
W_haut = W(2:n, :);

Phi = pinv(W_bas)*W_haut;
z = eig(Phi);

delta = log(abs(z));
f = 1/(2*pi)*angle(z);

end

