age = [3 4 5 6 7 8 9 11 12 14 15 16 17];
wing = [1.4 1.5 2.2 2.4 3.1 3.2 3.2 3.9 4.1 4.7 4.5 5.2 5.0];

%% Exercise 1
scatter(age, wing);

%% Exercise 2
reg = fitlm(age, wing);
plot(reg);

%% Exercise 3

%% Exercise 4
confInt = tinv(0.975, length(age) - 2) * ...
    sqrt(reg.MSE / sum((age - mean(age)) .^ 2));

%% Exercise 5
rSquared = reg.Rsquared.Ordinary;

%% Exercise 6
r = sum((age - mean(age)) .* (wing - mean(wing))) / ...
    (sqrt(sum((age - mean(age)) .^ 2)) * sqrt(sum((wing - mean(wing)) .^ 2)));

%% Exercise 7
noise = normrnd(0, 1, size(wing));
noisyWing = wing + noise;

reg = fitlm(age, noisyWing);
plot(reg);
rSquared = reg.Rsquared.Ordinary;