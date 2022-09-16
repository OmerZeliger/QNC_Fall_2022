%% simulate gaussian distribution
% doi: 10.1097/j.pain.0000000000001973

% approx statistics
mu = 30;
sigma = 7;
N = 1000;

samples = normrnd(mu, sigma, N, 1);

histogram(samples);
title({'Simulated labeled neuron size data', ...
    ['mean: ' num2str(mean(samples)) ', SD: ' num2str(std(samples))]});