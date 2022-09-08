%% 1
% Exercise #1: 
%   If someone gets a positive test, is it "statistically significant" at 
%   the p<0.05 level? Why or why not?

% No - the probability of getting a positive test given the null hypothesis 
% (no infection) is 0.05, which is not less than 0.05.


%% 2
% Exercise #2: 
%   What is the probability that if someone gets a positive test, 
%   that person is infected?

% calculate (infected & positive test) / (all positive test)
falsePositiveRate = 0.05;
infectionRates = 0:0.1:1;
chanceTruePositive = infectionRates ./ (infectionRates + (falsePositiveRate * (1 - infectionRates)));

disp([infectionRates; chanceTruePositive]);
