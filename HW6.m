reps = 1000;
differenceBetweenMeans = 0:5;
pVals = zeros(reps, length(differenceBetweenMeans));

mu = 0;
sigma = 4;
N = 50;

for diffIdx = 1:length(differenceBetweenMeans)
    X1 = normrnd(mu, sigma, N, reps);
    X2 = normrnd(mu + differenceBetweenMeans(diffIdx), sigma, N, reps);
    Sp = sqrt((var(X1) + var(X2)) ./ 2);
    tU = (mean(X1) - mean(X2)) ./ (Sp .* sqrt(2 ./ N));
    pVals(:, diffIdx) = 2 .* (1-tcdf(abs(tU), 2 * N - 2))';
end

disp(['number of significant values without correction: ' num2str(sum(pVals < 0.05))]);
disp(['number of significant values with Bonferroni correction: ' ...
    num2str(sum(pVals < (0.05 / reps)))]);
disp(['number of significant values with Benjamini?Hochberg procedure: ' ...
    num2str(sum(pVals < (0.05 ./ repmat((1:reps)', 1, length(differenceBetweenMeans)))))]);