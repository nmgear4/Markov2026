P = [ ...
    0    1/2  1/2  0    0    0    0    0;
    0    0    1/2  1/2  0    0    0    0;
    1/2  0    0    0    1/2  0    0    0;
    1/3  0    1/3  0    1/3  0    0    0;
    0    1/3  0    0    0    1/3  1/3  0;
    0    0    0    0    0    1    0    0;
    0    0    0    0    0    0    0    1;
    0    0    0    0    0    0    1    0];

d = 0.85;
n = 8;

G = d*P + (1-d)/n * ones(n);

A = G' - eye(n);
A(n,:) = ones(1,n);

rhs = zeros(n,1);
rhs(n) = 1;

pi = (A\rhs)';

R = 100;
T = 10^5;

rng(1);

state = ones(R,1);

CDF = cumsum(G,2);

counts = zeros(R,n);

recordTimes = [10^2 10^3 10^4 10^5];

errors = zeros(size(recordTimes));

recordNumber = 1;

for t = 1:T

    U = rand(R,1);

    currentCDF = CDF(state,:);

    nextState = 1 + sum(U > currentCDF,2);

    state = nextState;

    counts = counts + full(sparse( ...
        (1:R)',state,1,R,n));

    if recordNumber <= length(recordTimes) && ...
            t == recordTimes(recordNumber)

        hatPi = counts/t;

        maxError = max(abs(hatPi-pi),[],2);

        errors(recordNumber) = sqrt(mean(maxError.^2));

        recordNumber = recordNumber + 1;
    end

end

hatPi = counts/T;


fprintf('\nExact PageRank:\n');
disp(pi);

fprintf('\nSurfer 1 estimate at T = 10^5:\n');
disp(hatPi(1,:));

fprintf('\nMaximum error for surfer 1:\n');
disp(max(abs(hatPi(1,:)-pi)));

figure;

bar(1:n,[hatPi(1,:)' pi']);

xlabel('Page');
ylabel('Probability');

title('Surfer 1 estimate vs. PageRank');

legend('Surfer 1','Exact PageRank');

grid on;

figure;

loglog(recordTimes,errors,'o-','LineWidth',1.5);

hold on;

fit = polyfit(log10(recordTimes), ...
              log10(errors),1);

slope = fit(1);

fittedErrors = 10.^polyval(fit,log10(recordTimes));

loglog(recordTimes,fittedErrors,'--','LineWidth',1.5);

xlabel('T');
ylabel('RMS of max_i |hat{\pi}_i - \pi_i|');

title('Monte Carlo PageRank error');

legend('Simulation', ...
       sprintf('Fitted slope = %.3f',slope), ...
       'Location','southwest');

grid on;

fprintf('\n----------------------------------\n');
fprintf('Monte Carlo results\n');
fprintf('----------------------------------\n');

fprintf('T          RMS error\n');

for k = 1:length(recordTimes)
    fprintf('%-8d   %.6e\n', ...
        recordTimes(k),errors(k));
end

fprintf('\nFitted log-log slope = %.4f\n',slope);
