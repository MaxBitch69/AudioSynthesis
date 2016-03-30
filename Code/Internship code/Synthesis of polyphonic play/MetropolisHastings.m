function metro = MetropolisHastings(mu, sigma, N)
% Uses the Metropolis Hastings sampling to generate Time Difference (TD)
% Here, using a normal distribution of mean mu and standard deviation sigma
% N is the length of the list of samplings wanted

%% Metropolis-Hastings sampling
x = normrnd(mu, sigma); % Pick a random number following normal distribution
metro = [x]; % List containing the samples

for k = 1:N-1
   %y = x + normrnd(mu, 0.6); % Pick a random number following normal distribution
    y = normrnd(mu, sigma); % Pick a random number following normal distribution

    alpha = min( exp(-1/2*(y-mu).^2) / exp(-1/2*(x-mu).^2), 1); % Probability to go from x to y

%     if exp(-1/2*(x-mu).^2) > 0
%         alpha = min( exp(-1/2*(y-mu).^2) / exp(-1/2*(x-mu).^2), 1); % Probability to go from x to y
%     else
%         alpha = 1;
%     end
    u = normrnd(0, 1); % value to compare to alpha

    if u < alpha
        x = y;
    end
    
    metro = [metro x];
end

end
