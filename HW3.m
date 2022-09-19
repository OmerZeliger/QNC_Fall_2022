%% Binomial distribution
%
%
% Distribution of the number of successful outcomes for a given number 
%  of Bernoulli trials in a given experiment, defined by two parameters:
%  p = probability of success on each trial (constant across trials)
%  n = number of trials
%
%  This is also a DISCRETE distribution because by definition it only 
%  takes on specific, discrete values (in this case non-negative integers 
%  corresponding to the number of successes in n tries; thus, values 0:n).

% Choose some values for the parameters n and p
p = 0.7;
n = 1000;

% Generate random picks. Note that this is as if we did the Bernoulli
% trials as above, but then just counted the successes
outcome = binornd(n,p);

% Print out the results
fprintf('%d successes out of %d trials, simulated p = %.2f, empirical p = %.2f\n', ...
   n, outcome, p, outcome/n);

% The full probability distribution describes the probabilty of obtaining
%  each possible number of successes (k), given n and p. If we set n=10,
%  the the possible values of k are 0, 1, ..., 10. Now we use binornd
%  to simulate many different picks to get a full distribution
p = 0.7;
n = 10;       % number of "trials" per "experiment"
NumExperiments = 1000;     % number of "experiments"
outcomes = binornd(n,p,NumExperiments,1);

% Make histogram of all possible outcomes. We want bins centered on whole
% numbers so we offset the edges
edges = -0.5:10.5;
counts = histcounts(outcomes, edges);

% Show a bar plot of the simulated bionimal distribution
clf;
xs = edges(1:end-1)+diff(edges)/2;
bar(xs, counts);
title(sprintf('Histogram of binomial distribution, n=%d, p=%.2f, N=%d simulations', ...
   n,p,N));
xlabel(sprintf('Number of successes in %d tries', n));
ylabel('Count');

% Normalize it to make it a pdf. Here counts (the x-axis of the histogram)
%  is a DISCRETE variable, so we just have to add up the values
bar(xs, counts./sum(counts));
title(sprintf('PDF of binomial distribution, n=%d, p=%.2f, N=%d simulations', ...
   n,p,N));
xlabel(sprintf('Number of successes in %d tries', n));
ylabel('Probability');

% Compare it to the real binomial distribution, which we find using binopdf
Y = binopdf(xs,n,p);
hold on;
plot(xs,Y,'ro-', 'LineWidth', 2, 'MarkerSize', 10);
legend('Simulated', 'Theoretical')

% Fun animation showing different values of p and N. Note that, of course,
% the more picks you make of the random variable (hihger N), the more
% closely the outcomes (the "samples") match the true distribution.
for p = 0.1:0.1:0.9
   for N = round(logspace(1,5,10))
      
      % Get the pdf
      Y = binopdf(xs,n,p);
      
      % Get the random picks
      outcomes = binornd(n,p,N,1);
      counts = histcounts(outcomes, edges);
      
      % Set up the plot
      cla reset; hold on;
      
      % Show both
      bar(xs, counts./sum(counts));
      plot(xs,Y,'ro-', 'LineWidth', 2, 'MarkerSize', 10);
      
      % Labels, etc
      title(sprintf('p=%.1f, n=10, N=%d', p, N))
      xlabel(sprintf('Number of successes in %d tries', n));
      ylabel('Probability');
      legend('Simulated', 'Theoretical')
      axis([0 10 0 0.45])
      
      % Wait
      pause(0.15);
   end
end

% The cumulative distribution function is just the proability of obtaining
% an outcome that is equal to OR LESS THAN a particular value.
p = 0.7;
n = 10;        % number of "trials" 
N = 10000;     % number of "simulations"
outcomes = binornd(n,p,N,1);

% Make histogram of all possible outcomes. We want bins centered on whole
% numbers so we offset the edges
edges = -0.5:10.5;
counts = histcounts(outcomes, edges);

% Show a bar plot of the cdf (the cumulative sum of the normalized counts)
% Clear the figure
clf;
% Set up the x axis based on the edges
xs = edges(1:end-1)+diff(edges)/2;
% Compute the cumulative sum of the counts normalized by the total counts
%  (so it is a probability function and not just a histogram -- note that
%  the final value in the cdf should equal 1 because every value should be
%  equal to or less than that value).
scounts = cumsum(counts./sum(counts));
bar(xs, scounts);

% Compare it to the real binomial cumulative distribution, 
%  which we find using binocdf
Y = binocdf(xs,n,p);
hold on;
plot(xs,Y,'ro-', 'LineWidth', 2, 'MarkerSize', 10);

% Labels, etc
title(sprintf('Cumulative binomial distribution, p=%.1f, n=%d, N=%d', p, n, N))
xlabel(sprintf('Number of successes in %d tries', n));
ylabel('Cumulative probability');
legend('Simulated', 'Theoretical')


%% Exercise 1
% Assume that there are 10 quanta available in a nerve terminal, and for a 
% given release event each is released with a probability of 0.2. 
% For one such event, what is the probability that 0, 1, 2, 3, 4, 5, 6, 7, 
% 8, 9, or 10 quanta will be released?

p = 0.2;
n = 10;
N = 10000;

realDistribution = binopdf(0:n, n, p);
simulatedResults = binornd(n, p, N, 1);
simulatedDistribution = histcounts(simulatedResults,  -0.5:n+.5) / N;

hold on;
bar(0:n, simulatedDistribution);
plot(0:n, realDistribution, '-ro', 'LineWidth', 3);
hold off;

%% Exercise 2

totalQuanta = 14;
quantaReleased = 8;
p = 0:0.1:1;
probResult = binopdf(quantaReleased, totalQuanta, p);
plot(p, probResult, '-ro', 'LineWidth', 3);
[~,maxIndex] = max(probResult);
title(['Probability with max likelihood = ' num2str(p(maxIndex))]);

%% Exercise 3
totalQuanta = 14;
quantaReleased = [5 8];
p = 0:0.1:1;

probResult = zeros(length(quantaReleased), length(p));
for trialIdx = 1:length(quantaReleased)
    probResult(trialIdx, :) = binopdf(quantaReleased(trialIdx), totalQuanta, p);
end

% likelihood
subplot(1, 2, 1);
likelihood = prod(probResult);
plot(p, likelihood, '-ro', 'LineWidth', 3);
[~,maxIndex] = max(likelihood);
title(['Probability with max likelihood = ' num2str(p(maxIndex))]);

% log likelihood
subplot(1, 2, 2);
logLikelihood = sum(log(probResult));
plot(p, logLikelihood, '-ro', 'LineWidth', 3);
[~,maxIndex] = max(logLikelihood);
title(['Probability with max log likelihood = ' num2str(p(maxIndex))]);


%% Exercise 4
totalQuanta = 14;
quantaReleased = 0:totalQuanta;
counts = [0 0 3 7 10 19 26 16 16 5 5 0 0 0 0]; % guessed a reasonable number for 3 quanta released

p = 0:0.01:1;

probResult = zeros(length(quantaReleased), length(p));
for trialIdx = 1:length(quantaReleased)
    probResult(trialIdx, :) = binopdf(quantaReleased(trialIdx), totalQuanta, p);
end

%%TODO: fit to find p-hat
% likelihood
subplot(1, 2, 1);
likelihood = prod(power(probResult, counts'), 1);
plot(p, likelihood, '-ro', 'LineWidth', 1, 'MarkerSize', 3);
[~,maxIndex] = max(likelihood);
title(['Probability with max likelihood = ' num2str(p(maxIndex))]);

% log likelihood
subplot(1, 2, 2);
logLikelihood = counts * log(probResult);
plot(p, logLikelihood, '-ro', 'LineWidth', 1, 'MarkerSize', 3);
[~,maxIndex] = max(logLikelihood);
title(['Probability with max log likelihood = ' num2str(p(maxIndex))]);

%% Exercise 5
totalQuanta = 14;
quantaReleased = 7;
counts = 1;
originalProbability = 0.3;
p = 0:0.01:1;

probResult = zeros(length(quantaReleased), length(p));
for trialIdx = 1:length(quantaReleased)
    probResult(trialIdx, :) = binopdf(quantaReleased(trialIdx), totalQuanta, p);
end

% likelihood
likelihood = prod(power(probResult, counts'), 1);
[~,maxIndex] = max(likelihood);
disp(['Probability with max likelihood = ' num2str(p(maxIndex))]);
disp(['Likelihood of result given null hypothesis = ' num2str(likelihood(p == originalProbability))]);

% we can't conclude that changing had an effect because the likelihood is
% above the thresholf of 0.05

%% Bonus
%%TODO: the bonus