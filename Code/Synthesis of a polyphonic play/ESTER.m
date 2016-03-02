function bestK = ESTER(x, n, display)
% ESTER Method
% Finds the best K, so the best signal space's dimension
%
% Inputs:
%       - x: signal studied / windowed
%       - n: full space dimension (signal and noise spaces)
%       - display: 0 by default, if 1, plots inverse error between original
%       and synthesised signals
%
% Output:
%       - bestK: best signal space dimension

if nargin < 3
    display = 0;
end

%% Parameters
N = length(x); % Length of signal excerpted - Analysis window length

%% ESTER Method
range = 150;
error = zeros(range, 1);

for K = 1:range    
    % Computes synthesised signal
    % Method ESPRIT + LeastSquares
    [delta_e, f_e] = ESPRIT(x, n, K);
    [a_e, phi_e] = LeastSquares(x, delta_e, f_e);

    % Re-synthesis of the signal
    s = Synthesis(N, delta_e, f_e, a_e, phi_e);

    % Computes error
    error(K, 1) = 1/norm(x-s);
end

if display
    figure();
    plot(error);
    title('ESTER Method - Error between synthesised and original signal')
end

[~, bestK] = max(error);

end

