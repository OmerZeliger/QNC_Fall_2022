wing = [10.4 10.8 11.1 10.2 10.3 10.2 10.7 10.5 10.8 11.2 10.6 11.4];
tail = [7.4 7.6 7.9 7.2 7.4 7.1 7.4 7.2 7.8 7.7 7.8 8.3];

%% Exercise 1
scatter(wing, tail);

%% Exercise 2
meanWing = mean(wing);
meanTail = mean(tail);
manualR = sum((wing - meanWing) .* (tail - meanTail)) / ...
    sqrt(sum((wing - meanWing) .^ 2) * sum((tail - meanTail) .^ 2));

buildInR = corrcoef(wing, tail);

%% Exercise 3
standardError = sqrt((1 - r^2) / (length(wing) - 2));

z = 0.5 * log((1 + r) / (1 - r));
sd = sqrt(1 / (length(wing) - 3));
zSpaceConfInt = [z - (sd * norminv(1 - 0.025)) z + (sd * norminv(1 - 0.025))];
rConfidenceInterval = (exp(2*zSpaceConfInt) - 1) ./ (exp(2*zSpaceConfInt) + 1);

%% Exercise 4
tStat = r / standardError;
% reject the null hypothesis

%% Exercise 5
zr = 0.5 * log((1 + 0.75) / (1 - 0.75));
gamma = (zr - z) * sqrt(1 / (length(wing) - 3));
% do not reject thenull hypothesis

%% Exercise 6
theoreticalR = 0.5;
sampleSize = 1:1000;
tStat = theoreticalR ./ sqrt((1 - theoreticalR^2) ./ (sampleSize - 2));
necessarySampleSize = sampleSize(find(tStat >= norminv(1 - 0.025), 1));