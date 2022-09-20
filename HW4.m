% Compute confidence/credible intervals based on the four methods above
% for simulated data sampled from a population that is Gaussian distributed
% with mean ? = 10 and standard deviation ? = 2,
% for n=5, 10, 20, 40, 80, 160, 1000 at a 95% confidence level.


%% Simulate the data
mu = 10;
sigma = 2;
sampleSizes = [5, 10, 20, 40, 80, 160, 1000];
confidenceLevel = 0.95;
for n = sampleSizes
    sample = normrnd(mu, sigma, n, 1);
    
    %% 1. The simple, analytic approach with large n and/or known standard deviation.
    % use known standard deviation, aka sigma
    avg = mean(sample);
    sem = sigma / sqrt(n);
    
    % use z-score calculator
    zscore = norminv((1 - confidenceLevel) / 2);
    confidenceInterval = abs(sem * zscore);
    disp(['sample size = ' num2str(n) ', z-tested confidence interval is ' ...
        num2str(avg - confidenceInterval) ' - ' num2str(avg + confidenceInterval)]);
    
    %% 2. The simple, analytic approach with small n and unknown population standard deviation
    % use calculated standard deviation
    avg = mean(sample);
    sem = std(sample) / sqrt(n);
    
    % use t-score calculator
    tscore = tinv((1 - confidenceLevel) / 2 ,n-1);
    confidenceInterval = abs(sem * tscore);
    disp(['sample size = ' num2str(n) ', t-tested confidence interval is ' ...
        num2str(avg - confidenceInterval) ' - ' num2str(avg + confidenceInterval)]);
    
    %% 3. Bootstrapped confidence intervals
    % resample with replacement 1000 times, save means
    numReps = 1000;
    bootstrappedMeans = zeros(numReps, 1);
    for rep = 1:numReps
        bootstrappedMeans(rep) = mean(sample(randi(n, [n 1])));
    end
    
    % get bootstrapped distribution statistics
    avg = mean(bootstrappedMeans);
    sd = std(bootstrappedMeans);
    
    % use z-score calculator
    zscore = norminv((1 - confidenceLevel) / 2);
    confidenceInterval = abs(sd * zscore);
    disp(['sample size = ' num2str(n) ', bootstrapped confidence interval is ' ...
        num2str(avg - confidenceInterval) ' - ' num2str(avg + confidenceInterval)]);
    
    %% 4. Bayesian credible intervals
    % flat prior, so bayesian credible interval is identical to confidence
    % interval for gaussian dirtsibution (?)
    
    % use sample mean and sem
    avg = mean(sample);
    sem = std(sample) / sqrt(n);
    
    % use z-score calculator
    zscore = norminv((1 - confidenceLevel) / 2);
    confidenceInterval = abs(sem * zscore);
    disp(['sample size = ' num2str(n) ', credible interval is ' ...
        num2str(avg - confidenceInterval) ' - ' num2str(avg + confidenceInterval)]);
    
end


