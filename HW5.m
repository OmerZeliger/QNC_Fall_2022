% simulating null hypothesis
numTrials = 100; % assume 100 trials/session
numSimulations = 10000;
pupilData = normrnd(0, 1, [numTrials, numSimulations]);
LCData = normrnd(1, 1, [numTrials numSimulations]);

NHCorrelations = corr(pupilData, LCData);
NHCorrelations = NHCorrelations(logical(eye(numSimulations)));
figure 1; histogram(NHCorrelations);

% calculate correlation coefficients necessary for an experiment with the given power
power = 0.8;
NHParameters = [mean(NHCorrelations) std(NHCorrelations)];
effectSizes = (0:0.01:1)';
alternativeMeans = NHParameters(2) * effectSizes;
samplesRequired = sampsizepwr('z', NHParameters, alternativeMeans, power);

figure 2; plot(effectSizes, samplesRequired);
ylabel(['samples needed for power of ' num2str(power)]);
xlabel('effect size');